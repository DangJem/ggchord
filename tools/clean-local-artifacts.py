#!/usr/bin/env python3
"""Inventory disposable local artifacts; delete only with --clean.

Run from any directory. Package sources, tracked files, original inputs, audit
images, frozen docs and Git history are never cleanup targets.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def tracked_files():
    raw = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT)
    return [name for name in raw.decode().split('\0') if name]


def source_hashes(names):
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in names if (ROOT / name).is_file()}


def usage(path):
    files = [path] if path.is_file() else (
        Path(parent) / name for parent, _, names in os.walk(path, followlinks=False)
        for name in names)
    count = logical = allocated = 0
    for file in files:
        info = file.lstat()
        count += 1
        logical += info.st_size
        allocated += info.st_blocks * 512
    return dict(files=count, logical_bytes=logical, allocated_bytes=allocated)


def candidates():
    paths = [ROOT / '.site', ROOT / '.codex-tmp/ocr-libs']
    paths += [p for p in ROOT.glob('*.Rcheck') if (p / '00check.log').is_file()]
    disposable = ('check', 'library', 'download', 'download-check',
                  'reissue-download', 'devices', 'devices-final')
    for parent in (ROOT / '.codex-tmp').glob('*'):
        if re.fullmatch(r'(release|cran)-v\d+\.\d+\.\d+', parent.name):
            if parent.is_symlink():
                raise SystemExit('Refusing symlink artifact parent: ' + str(parent))
            paths += [parent / name for name in disposable]
    return sorted({path for path in paths if path.exists()})


def preserve_records(path):
    records = []
    if path.name == '.site':
        records += [path / name for name in ('source-snapshot.json', 'build-manifest.json',
                    'render-processes.json', 'build.log', 'en-build.log', 'zh-build.log')]
    checks = [path] if path.name.endswith('.Rcheck') else list(path.glob('*.Rcheck'))
    for check in checks:
        records += [check / name for name in ('00check.log', '00install.out')]
        records += list(check.glob('*.Rout')) + list(check.glob('tests/*.Rout*'))
    saved = []
    for file in records:
        if not file.is_file() or file.is_symlink():
            continue
        digest = hashlib.sha256(file.read_bytes()).hexdigest()
        relative = file.relative_to(ROOT)
        destination = ROOT / '.codex-tmp/cleanup-records' / relative
        if destination.is_symlink() or any(part.is_symlink() for part in destination.parents):
            raise SystemExit('Refusing symlink evidence destination: ' + str(destination))
        if destination.exists() and destination.read_bytes() != file.read_bytes():
            destination = destination.with_name(destination.name + '.' + digest[:12])
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(file, destination)
        saved.append(dict(file=str(destination.relative_to(ROOT)), sha256=digest,
                          bytes=destination.stat().st_size))
    return saved


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--clean', action='store_true', help='Remove listed artifacts after safety checks')
    parser.add_argument('--report', type=Path, help='Write a new JSON record; never overwrite an existing file')
    args = parser.parse_args()
    if args.report and args.report.exists():
        raise SystemExit('Report already exists; choose a new path')
    tracked = tracked_files()
    paths = candidates()
    for path in paths:
        relative = path.relative_to(ROOT).as_posix()
        if not path.is_dir():
            raise SystemExit('Refusing artifact that is not a directory: ' + relative)
        if path.is_symlink() or any(part.is_symlink() for part in path.parents if part != ROOT):
            raise SystemExit('Refusing symlink artifact: ' + relative)
        if any(name == relative or name.startswith(relative + '/') for name in tracked):
            raise SystemExit('Refusing artifact containing tracked files: ' + relative)
        ignored = subprocess.run(['git', 'check-ignore', '-q', '--', relative], cwd=ROOT)
        if ignored.returncode != 0:
            raise SystemExit('Refusing artifact not covered by .gitignore: ' + relative)
    before = source_hashes(tracked)
    refs = subprocess.check_output(['git', 'show-ref'], cwd=ROOT)
    rows = [dict(path=str(path.relative_to(ROOT)), **usage(path)) for path in paths]
    saved = []
    if args.clean:
        for path in paths:
            saved += preserve_records(path)
            shutil.rmtree(path)
        if source_hashes(tracked) != before:
            raise SystemExit('Tracked contents changed during cleanup; inspect before continuing')
        if subprocess.check_output(['git', 'show-ref'], cwd=ROOT) != refs:
            raise SystemExit('Git refs changed during cleanup; inspect before continuing')
    result = dict(mode='cleaned' if args.clean else 'inventory only', artifacts=rows,
                  candidate_allocated_bytes=sum(row['allocated_bytes'] for row in rows),
                  preserved_records=saved, tracked_file_count=len(before),
                  tracked_contents_sha256=hashlib.sha256(json.dumps(before, sort_keys=True).encode()).hexdigest(),
                  git_refs_sha256=hashlib.sha256(refs).hexdigest(),
                  tracked_contents_preserved=True, git_refs_preserved=True,
                  scope='Only known ignored build/cache copies; no Git history rewrite, package checks or image regeneration')
    if args.report:
        with args.report.open('x') as handle:
            json.dump(result, handle, indent=2)
            handle.write('\n')
    for row in rows:
        print(f"{row['allocated_bytes'] / 2**20:8.2f} MiB  {row['path']}")
    print(f"{result['mode']}: {result['candidate_allocated_bytes'] / 2**20:.2f} MiB")
    print(f"Preserved tracked contents ({len(before)} files) and Git refs; kept {len(saved)} small evidence files.")


if __name__ == '__main__':
    main()
