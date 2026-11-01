#!/usr/bin/env python3
"""Scoped source/site checks. Does not build the package or execute R examples."""
import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
from documentation_lib import ROOT, WEB, PAGES, CAPTIONS, expand

class Document(HTMLParser):
    def __init__(self, path):
        super().__init__(convert_charrefs=True)
        self.ids = []
        self.links = []
        self.images = []
        self.lang = None
        self.feed(path.read_text())
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if 'id' in a: self.ids.append(a['id'])
        if tag == 'html': self.lang = a.get('lang')
        if tag in ['a', 'link', 'script', 'img']:
            for key in ['href', 'src']:
                if a.get(key): self.links.append(a[key])
        if tag == 'img': self.images.append(a)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-only', action='store_true')
    args = parser.parse_args()
    errors = []
    def require(ok, message):
        if not ok: errors.append(message)
    slugs = [p['slug'] for p in PAGES]
    require(len(slugs) == len(set(slugs)), 'Duplicate page slug')
    def blocks(text):
        return re.findall(r'```[^\n]*\n(.*?)\n```', text, re.S)
    for p in PAGES:
        texts = [(WEB/lang/(p['slug']+'.md')).read_text() for lang in ['en','zh']]
        for lang, text in zip(['en','zh'],texts):
            require('中文' not in text, f'{lang}/{p["slug"]}: inconsistent language label')
            require(bool(p['title'][lang]), f'{p["slug"]}: missing localized title')
        require(re.findall(r'\{#([^}]+)\}', texts[0]) == re.findall(r'\{#([^}]+)\}', texts[1]), f'{p["slug"]}: unmatched anchors')
        require(re.findall(r'\{\{[^}]+\}\}', texts[0]) == re.findall(r'\{\{[^}]+\}\}', texts[1]), f'{p["slug"]}: unmatched shared code/asset markers')
        require(blocks(expand(texts[0], 'en')) == blocks(expand(texts[1], 'zh')), f'{p["slug"]}: mismatched executable code')
    for lang in ['en','zh']:
        require({f.stem for f in (WEB/lang).glob('*.md')} == set(slugs), f'{lang}: orphan or missing source page')
    subprocess.run(['python3','tools/sync-tutorials.py','--check'],cwd=ROOT,check=True)
    with (WEB/'assets/plots.csv').open() as f:
        assets = list(csv.DictReader(f))
    require({a['id'] for a in assets} == set(CAPTIONS), 'Figure manifest/captions mismatch')
    for a in assets:
        path = ROOT/a['file']
        require(path.is_file(),f'Missing figure {path}')
        if path.is_file(): require(hashlib.md5(path.read_bytes()).hexdigest() == a['md5'],f'Figure hash mismatch: {path}')
        require(a['code_md5'] == hashlib.md5((ROOT/'tools/documentation-examples.R').read_bytes()).hexdigest(), 'Shared plot code changed: regenerate documentation figures')
        require(float(a['width'])>0 and float(a['height'])>0 and a['units']=='in' and a['device']=='png', f'Invalid figure dimensions: {a["id"]}')
    aliases = set()
    rd_topics = set()
    for rd in (ROOT/'man').glob('*.Rd'):
        rd_topics.add(rd.stem)
        aliases.update(re.findall(r'\\alias\{([^}]+)\}',rd.read_text()))
    exports = set(re.findall(r'^export\(([^)]+)\)',(ROOT/'NAMESPACE').read_text(),re.M))
    require(not (exports-aliases), f'Exports without Rd: {sorted(exports-aliases)}')
    for lang in ['en','zh']:
        targets = {p['slug']: set(re.findall(r'\{#([^}]+)\}', (WEB/lang/(p['slug']+'.md')).read_text())) for p in PAGES}
        for p in PAGES:
            for url in re.findall(r'\]\(([^)\s]+)\)', expand((WEB/lang/(p['slug']+'.md')).read_text(),lang)):
                split = urlsplit(url)
                if split.scheme or split.netloc: continue
                if split.path.startswith('reference/'):
                    require(Path(split.path).stem in rd_topics, f'{lang}/{p["slug"]}: unknown help target {url}')
                elif split.path.endswith('.html'):
                    slug = Path(split.path).stem
                    require(slug in targets,f'{lang}/{p["slug"]}: unknown topic {url}')
                    if split.fragment and slug in targets:
                        require(split.fragment in targets[slug],f'{lang}/{p["slug"]}: unknown source anchor {url}')
                elif split.path.startswith('../assets/'):
                    name = Path(split.path).name
                    origin = ROOT/'man/figures'/name if name.endswith(('.png','.svg')) else WEB/'assets'/name
                    if name=='NEWS.md': origin=ROOT/'NEWS.md'
                    if name=='documentation-examples.R': origin=ROOT/'tools/documentation-examples.R'
                    require(origin.exists(), f'{lang}/{p["slug"]}: unknown asset {url}')
    icon = json.loads((WEB/'assets/icon.json').read_text())
    for record in icon['assets']+[icon['favicon']]:
        path=ROOT/record['file']
        require(path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256'],f'Icon hash mismatch: {path}')
    config = (ROOT/'_pkgdown.yml').read_text()
    require(all(re.search(r'^      - '+re.escape(t)+r'$',config,re.M) for t in rd_topics), 'Fallback pkgdown reference omits Rd topics')
    # README/vignette resources and links are source paths, not the frozen site.
    for path in [ROOT/'README.md',ROOT/'README-Hans.md',ROOT/'MIGRATION.md',ROOT/'MIGRATION-Hans.md',*(ROOT/'vignettes').glob('*.Rmd')]:
        text = path.read_text()
        ids = set(re.findall(r'\{#([^}]+)\}',text))
        for url in re.findall(r'\]\(([^)\s]+)\)',text):
            split = urlsplit(url)
            if split.scheme or split.netloc: continue
            if split.path:
                require((path.parent/unquote(split.path)).exists(), f'{path.name}: missing source link {url}')
            elif split.fragment:
                require(unquote(split.fragment) in ids,f'{path.name}: missing chapter {url}')
    if not args.source_only:
        out = ROOT/'.site'
        require(out.is_dir(), 'Build .site first')
        docs = {p.resolve():Document(p) for p in out.rglob('*.html')}
        for path,doc in docs.items():
            rel = str(path.relative_to(out))
            require(not [i for i,n in Counter(doc.ids).items() if n>1],f'{rel}: duplicate element IDs')
            if rel.startswith(('en/','zh/')):
                require(doc.lang == ('zh-Hans' if rel.startswith('zh/') else 'en'),f'{rel}: language metadata')
                require('content' in doc.ids,f'{rel}: missing skip target')
            for img in doc.images: require('alt' in img,f'{rel}: image missing alt')
            for url in doc.links:
                split = urlsplit(url)
                if split.scheme or split.netloc: continue
                target = (path.parent/unquote(split.path)).resolve() if split.path else path
                if target.is_dir(): target=target/'index.html'
                require(target.is_relative_to(out),f'{rel}: link escapes site {url}')
                require(target.exists(),f'{rel}: broken local link {url}')
                if split.fragment and target in docs:
                    require(unquote(split.fragment) in docs[target].ids,f'{rel}: broken anchor {url}')
        for lang in ['en','zh']:
            expected = {f'{s}.html' for s in slugs} | {f'reference/{s}.html' for s in rd_topics}
            actual = {str(p.relative_to(out/lang)) for p in (out/lang).rglob('*.html')}
            require(actual==expected,f'{lang}: generated page coverage')
            index = json.loads((out/lang/'search.json').read_text())
            require({p['url'] for p in index} == expected,f'{lang}: incomplete search coverage')
            require(any('track' in p['text'] for p in index),f'{lang}: empty search content')
        require(any('标签' in p['text'] for p in json.loads((out/'zh/search.json').read_text())), 'Han search content absent')
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f'PASS: {len(PAGES)} language pairs, identical shared code, assets, {len(exports)} exports / {len(rd_topics)} Rd topics'+(' and generated links/anchors/search.' if not args.source_only else '.'))

if __name__ == '__main__':
    main()
