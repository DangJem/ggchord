"""Verify the clipboard text extracted from independently runnable examples."""
import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


class CopiedBlocks(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.stack = []
        self.blocks = []
        self.parts = None

    def handle_starttag(self, tag, attrs):
        if tag in ('area','base','br','col','embed','hr','img','input','link','meta','param','source','track','wbr'):
            return
        kind = dict(attrs).get('data-example-kind') or (self.stack[-1][1] if self.stack else None)
        self.stack.append((tag, kind))
        if kind and tag == 'code' and any(t == 'pre' for t, _ in self.stack):
            self.parts = []

    def handle_endtag(self, tag):
        if self.parts is not None and tag == 'code':
            self.blocks.append((self.stack[-1][1], ''.join(self.parts)))
            self.parts = None
        for index in range(len(self.stack)-1, -1, -1):
            if self.stack[index][0] == tag:
                del self.stack[index:]
                break

    def handle_data(self, text):
        if self.parts is not None:
            self.parts.append(text)


def verify(out, workers, library=None):
    locations = defaultdict(list)
    for page in sorted(out.glob('*/*.html')):
        if page.parent.name not in ('en', 'zh'):
            continue
        reader = CopiedBlocks()
        reader.feed(page.read_text())
        for i, (kind, code) in enumerate(reader.blocks, 1):
            locations[(kind, code)].append(f'{page.relative_to(out)}#{i}')
    if not any(kind == 'plot' for kind, _ in locations):
        raise ValueError('No independently runnable plots found in the generated pages')
    with tempfile.TemporaryDirectory(prefix='ggchord-copy-') as tmp:
        folder = Path(tmp)
        if library is None:
            library = folder / 'library'
            library.mkdir()
            subprocess.run(['R','CMD','INSTALL','--no-docs','--no-html','-l',str(library),str(ROOT)],
                           cwd=folder, check=True, stdout=subprocess.DEVNULL)
        env = dict(os.environ, R_LIBS=os.pathsep.join(filter(None, [str(library), os.environ.get('R_LIBS')])))

        def run(item):
            (kind, code), pages = item
            with tempfile.TemporaryDirectory(dir=folder) as work:
                proc = subprocess.run(['Rscript','--vanilla',str(ROOT/'tools/website/run-copy.R'),kind],
                                      input=code, text=True, cwd=work, env=env,
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if proc.returncode:
                print('FAILED', pages, proc.stdout[-800:], flush=True)
            return dict(pages=pages, kind=kind, sha256=hashlib.sha256(code.encode()).hexdigest(),
                        success=proc.returncode == 0, **({'error':proc.stdout[-800:]} if proc.returncode else {}))

        with ThreadPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(run, locations.items()))
    errors = [r for r in results if not r['success']]
    (out/'copy-examples.json').write_text(json.dumps(results, ensure_ascii=False, indent=2)+'\n')
    if errors:
        raise SystemExit(f'{len(errors)} of {len(results)} copied R blocks failed')
    print(f'{sum(len(r["pages"]) for r in results)} copy controls verified in {len(results)} independent R sessions.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'.site')
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--library', type=Path)
    args = parser.parse_args()
    verify(args.output.resolve(), args.workers, args.library)
