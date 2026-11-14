#!/usr/bin/env python3
"""Apply two exact post-build corrections without reusing stale article inputs."""
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = '4b059e3'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
site = ROOT/'.site'
snapshot_path = site/'source-snapshot.json'
previous_sha = sha(snapshot_path)
snapshot = json.loads(snapshot_path.read_text())
assert snapshot['status'] == 'verified unchanged across rendering and finalization'
changes = {
    'tools/finalize-website.R': (
        "    id <- if(title %in% c('绘图模块', 'Plotting modules')) 'site-visualizations'",
        "    id <- if(title %in% c('开始使用', 'Get started')) 'site-start'\n      else if(title %in% c('绘图模块', 'Plotting modules')) 'site-visualizations'"),
    'tools/build-website.py': (
        "    for p in (ROOT/'man/figures').iterdir():\n        if p.is_file() and (p.stem in CAPTIONS or p.stem=='logo'):shutil.copy2(p,assets/p.name)",
        "    # Copy the declared inventory explicitly; directory enumeration can omit an asset.\n    for key in [*CAPTIONS, 'logo']:\n        p=ROOT/'man/figures'/(key+'.png')\n        shutil.copy2(p,assets/p.name)")}
actual = {name for name, digest in snapshot['render_inputs'].items() if sha(ROOT/name) != digest}
assert actual == set(changes), actual
for name, (before, after) in changes.items():
    original = subprocess.check_output(['git', 'show', BASE+':'+name], cwd=ROOT)
    assert hashlib.sha256(original).hexdigest() == snapshot['render_inputs'][name]
    text = original.decode()
    assert text.count(before) == 1
    assert (ROOT/name).read_text() == text.replace(before, after, 1), name
for name, digest in snapshot['archive_inputs'].items():
    assert sha(ROOT/name) == digest, name

html_changes = {}
for lang, old in [('en', 'dropdown-get-started'), ('zh', 'dropdown--')]:
    for path in (site/lang).rglob('*.html'):
        text = path.read_text()
        pattern = r'((?:id|aria-labelledby)=")'+re.escape(old)+r'(")'
        fresh, count = re.subn(pattern, lambda m: m[1]+'site-start'+m[2], text)
        if count:
            assert count == 2, path
            assert re.sub(r'((?:id|aria-labelledby)=")site-start(")',
                          lambda m: m[1]+old+m[2], fresh) == text
            html_changes[str(path.relative_to(site))] = dict(before=sha(path))
            path.write_text(fresh)
            html_changes[str(path.relative_to(site))]['after'] = sha(path)
assert html_changes

copied = {}
captions = json.loads((ROOT/'website/assets/captions.json').read_text())
for key in [*captions, 'logo']:
    source = ROOT/'man/figures'/(key+'.png')
    target = site/'assets'/source.name
    if not target.exists() or sha(source) != sha(target):
        copied[source.name] = sha(source)
        shutil.copy2(source, target)
    assert sha(source) == sha(target), key
assert set(copied) == {'showcasecollinear.png'}, copied

record = dict(previous_snapshot_sha256=previous_sha,
              source_files={name: sha(ROOT/name) for name in changes},
              html_files=len(html_changes), html_sha256=html_changes, copied_assets=copied,
              scope='Exact new start-menu id/ARIA normalization and missing PNG copy only; all article bodies, executable examples and archive inputs retained. No native content rerender or browser acceptance.')
(ROOT/'tools/site-audit/manual-static-correction.json').write_text(json.dumps(record, indent=2)+'\n')
for name in changes:
    snapshot['render_inputs'][name] = sha(ROOT/name)
snapshot.setdefault('derived_refreshes', []).append({k: v for k, v in record.items() if k != 'html_sha256'})
snapshot_path.write_text(json.dumps(snapshot, indent=2)+'\n')
manifest_path = site/'build-manifest.json'
manifest = json.loads(manifest_path.read_text())
manifest['source_sha256'].update(record['source_files'])
manifest['derived_refreshes'] = snapshot['derived_refreshes']
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
print('PASS: exact start-menu identifiers normalized in', len(html_changes), 'HTML files; one missing PNG copied')
