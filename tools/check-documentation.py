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
from documentation_lib import ROOT, WEB, PAGES, CAPTIONS, RECIPES, expand, snippets, localize_r_comments
import importlib.util

fingerprint_spec = importlib.util.spec_from_file_location('plot_fingerprints', ROOT/'tools/documentation-fingerprints.py')
fingerprint_module = importlib.util.module_from_spec(fingerprint_spec)
fingerprint_spec.loader.exec_module(fingerprint_module)

class Document(HTMLParser):
    def __init__(self, path):
        super().__init__(convert_charrefs=True)
        self.ids = []
        self.links = []
        self.images = []
        self.codeblocks = []
        self._code_parts = None
        self.main_text = []
        self._in_main = False
        self.lang = None
        self.classes = set()
        self.headings = []
        self._heading = None
        self._heading_parts = []
        self.feed(path.read_text())
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        self.classes.update(a.get('class','').split())
        if tag in ['h1','h2','h3']:
            self._heading = tag
            self._heading_parts = []
        if tag == 'main': self._in_main = True
        if tag == 'pre': self._code_parts = []
        if 'id' in a: self.ids.append(a['id'])
        if tag == 'html': self.lang = a.get('lang')
        if tag in ['a', 'link', 'script', 'img']:
            for key in ['href', 'src']:
                if a.get(key): self.links.append(a[key])
        if tag == 'img': self.images.append(a)
    def handle_data(self, data):
        if self._heading is not None: self._heading_parts.append(data)
        if self._in_main: self.main_text.append(data)
        if self._code_parts is not None: self._code_parts.append(data)
    def handle_endtag(self, tag):
        if tag == self._heading:
            self.headings.append(''.join(self._heading_parts).strip())
            self._heading = None
        if tag == 'main': self._in_main = False
        if tag == 'pre' and self._code_parts is not None:
            self.codeblocks.append(''.join(self._code_parts))
            self._code_parts = None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-only', action='store_true')
    args = parser.parse_args()
    errors = []
    def require(ok, message):
        if not ok: errors.append(message)
    slugs = [p['slug'] for p in PAGES]
    require(len(slugs) == len(set(slugs)), 'Duplicate page slug')
    navigation = json.loads((WEB/'navigation.json').read_text())
    entries = [slug for group in navigation for slug in group['pages']]
    require(set(entries) == set(slugs) - {'index'} and len(entries) == len(set(entries)),
            'Site guide must cover every current learning page once and exclude the home page')
    require({'chord','collinear','plasmid'} <= {group['id'] for group in navigation},
            'Missing visualization task module')
    def blocks(text):
        return re.findall(r'```[^\n]*\n(.*?)\n```', text, re.S)
    def executable_text(block):
        return '\n'.join(line for line in block.replace('-zh.R', '.R').splitlines()
                          if not re.match(r'^\s*#', line))
    def executable_blocks(text):
        return [executable_text(block) for block in blocks(text)]
    for p in PAGES:
        texts = [(WEB/lang/(p['slug']+'.md')).read_text() for lang in ['en','zh']]
        for lang, text in zip(['en','zh'],texts):
            require('中文' not in text, f'{lang}/{p["slug"]}: inconsistent language label')
            require(bool(p['title'][lang]), f'{p["slug"]}: missing localized title')
            require('{{' not in expand(text,lang),f'{lang}/{p["slug"]}: unexpanded source marker')
            require(all(re.search(r'\{#[^}]+\}$',heading) for heading in re.findall(r'^#{2,4} .+$',text,re.M)),
                f'{lang}/{p["slug"]}: headings need explicit paired anchors')
        require(re.findall(r'\{#([^}]+)\}', texts[0]) == re.findall(r'\{#([^}]+)\}', texts[1]), f'{p["slug"]}: unmatched anchors')
        require(re.findall(r'\{\{[^}]+\}\}', texts[0]) == re.findall(r'\{\{[^}]+\}\}', texts[1]), f'{p["slug"]}: unmatched shared code/asset markers')
        require(executable_blocks(expand(texts[0], 'en')) == executable_blocks(expand(texts[1], 'zh')), f'{p["slug"]}: mismatched executable code')
    for lang in ['en','zh']:
        require({f.stem for f in (WEB/lang).glob('*.md')} == set(slugs), f'{lang}: orphan or missing source page')
    for filename in ['README.md', 'README-Hans.md']:
        require(not re.search(r'SnapGene|Geneious', (ROOT/filename).read_text(), re.I),
                filename + ': commercial product name in design references')
    comparisons = json.loads((WEB/'parameter-comparisons.json').read_text())
    for row in comparisons:
        require(all(key in snippets() for key in row['snippets']), row['id'] + ': missing comparison dependency')
        require(row['page'] in slugs, row['id'] + ': unknown comparison lesson')
        if row['baseline']:
            require(row['baseline'] in CAPTIONS, row['id'] + ': unknown baseline')
    subprocess.run(['python3','tools/generate-task-index.py','--check'],cwd=ROOT,check=True)
    subprocess.run(['python3','tools/generate-parameter-atlas.py','--check'],cwd=ROOT,check=True)
    subprocess.run(['python3','tools/generate-api-curriculum.py','--check'],cwd=ROOT,check=True)
    reviews=json.loads((WEB/'teaching-reviews.json').read_text())
    inventory=json.loads((ROOT/'tools/teaching-audit/api-signatures.json').read_text())
    require(set(reviews)=={e['name'] for e in inventory}, 'Every export needs an explicit curriculum decision')
    for lang in ['en','zh']:
        compatibility=(WEB/lang/'compatibility.md').read_text()
        documented_geoms=set(re.findall(r'^\| \[\x60(geom_\w+)\(\)\x60\]',compatibility,re.M))
        require(documented_geoms=={e['name'] for e in inventory if e['name'].startswith('geom_')},lang+': compatibility table must retain each individual geom row')
    contracts=json.loads((ROOT/'tools/teaching-audit/parameter-contracts.json').read_text())
    require({(r['api'],r['parameter']) for r in contracts} ==
            {(e['name'],p) for e in inventory for p in (e['parameters'] or {})}, 'Incomplete parameter contracts')
    runtime_record=json.loads((ROOT/'tools/teaching-audit/lesson-execution.json').read_text())
    require(runtime_record.get('shared_code_sha256')==hashlib.sha256((ROOT/'tools/documentation-examples.R').read_bytes()).hexdigest(),'Lesson execution uses stale shared R source')
    for path, digest in runtime_record.get('engine_sha256',{}).items():
        require(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'Lesson execution uses stale engine: '+path)
    execution=runtime_record['pages']
    expanded_for_validation={p['slug']:expand((WEB/'en'/(p['slug']+'.md')).read_text(),'en') for p in PAGES}
    all_code='\n'.join('\n'.join(re.findall(r'```r\n(.*?)```',text,re.S)) for text in expanded_for_validation.values())
    for entry in inventory:
        name=entry['name']; review=reviews[name]
        require(all(review.get(k) for k in ['purpose','reason','page','anchor','modes','scope']), name+': incomplete semantic review')
        if entry['owner']=='ggchord' and not review.get('alias_of'):
            require(bool(re.search(r'(?<![\w.])'+re.escape(name)+r'\s*\(',all_code)),name+': missing canonical lesson call')
            runtime_pages=[slug for slug, result in execution.items() if result['result']=='success' and re.search(r'(?<![\w.])'+re.escape(name)+r'\s*\(',expanded_for_validation[slug])]
            require(bool(runtime_pages),name+': no executed lesson call')
    # Occurrence and execution are prerequisites, never automatic depth upgrades.
    subprocess.run(['python3','tools/reference-translations.py','--check'],cwd=ROOT,check=True)
    subprocess.run(['python3','tools/sync-tutorials.py','--check'],cwd=ROOT,check=True)
    with (WEB/'assets/plots.csv').open() as f:
        assets = list(csv.DictReader(f))
    require({a['id'] for a in assets} == set(CAPTIONS), 'Figure manifest/captions mismatch')
    by_id = {a['id']: a for a in assets}
    for row in comparisons:
        asset = by_id.get(row['id'], {})
        require(all(asset.get(key) == row[key] for key in ['object', 'coordinate']),
                row['id'] + ': comparison object or mode differs from rendered asset')
        require(all(float(asset.get(key, 0)) == float(row[key]) for key in ['width', 'height']),
                row['id'] + ': comparison dimensions differ from rendered asset')
    for group in navigation:
        expected = {'chord': 'chord', 'collinear': 'collinear', 'plasmid': 'circular'}.get(group['id'])
        if not expected: continue
        for slug in group['pages']:
            source = (WEB/'en'/(slug+'.md')).read_text()
            for figure in re.findall(r'\{\{figure:([a-z]+)\}\}', source):
                require(by_id.get(figure, {}).get('coordinate') == expected,
                        f'{slug}: figure {figure} must use {expected} coordinates')
    for page in PAGES:
        source = (WEB/'en'/(page['slug']+'.md')).read_text()
        for keys in re.findall(r'\{\{code:([a-z,]+)\}\}', source):
            if keys in {'showcasechord', 'showcasecollinear', 'showcaseplasmid'}: continue
            lines = sum(len(snippets()[key].splitlines()) for key in keys.split(','))
            require(lines <= 22, f'{page["slug"]}: simplify long teaching block {keys} ({lines} lines)')
    current_fingerprints = fingerprint_module.fingerprints()
    for a in assets:
        path = ROOT/a['file']
        require(path.is_file(),f'Missing figure {path}')
        if path.is_file(): require(hashlib.md5(path.read_bytes()).hexdigest() == a['md5'],f'Figure hash mismatch: {path}')
        require(a.get('dependency_md5') == current_fingerprints['assets'].get(a['id']),
                f'Plot dependencies changed: regenerate {a["id"]}')
        require(float(a['width'])>0 and float(a['height'])>0 and a['units']=='in' and a['device']=='png', f'Invalid figure dimensions: {a["id"]}')
    aliases = set()
    rd_topics = set()
    for rd in (ROOT/'man').glob('*.Rd'):
        rd_topics.add(rd.stem)
        aliases.update(re.findall(r'\\alias\{([^}]+)\}',rd.read_text()))
    exports = set(re.findall(r'^export\(([^)]+)\)',(ROOT/'NAMESPACE').read_text(),re.M))
    require(not (exports-aliases), f'Exports without Rd: {sorted(exports-aliases)}')
    for lang in ['en','zh']:
        targets = {p['slug']: {a or b for a,b in re.findall(r'\{#([^}]+)\}|\bid="([^"]+)"', (WEB/lang/(p['slug']+'.md')).read_text())} for p in PAGES}
        for p in PAGES:
            for url in re.findall(r'\]\(([^)\s]+)\)', expand((WEB/lang/(p['slug']+'.md')).read_text(),lang,figures='../assets')):
                split = urlsplit(url)
                if split.scheme or split.netloc: continue
                if split.path.startswith('reference/'):
                    require(split.path == 'reference/index.html' or Path(split.path).stem in rd_topics, f'{lang}/{p["slug"]}: unknown help target {url}')
                elif split.path.endswith('.html'):
                    slug = Path(split.path).stem
                    require(slug in targets,f'{lang}/{p["slug"]}: unknown topic {url}')
                    if split.fragment and slug in targets:
                        require(split.fragment in targets[slug],f'{lang}/{p["slug"]}: unknown source anchor {url}')
                elif split.path.startswith('../assets/'):
                    name = Path(split.path).name
                    origin = ROOT/'man/figures'/name if name.endswith(('.png','.svg')) else WEB/'assets'/name
                    if name in ['NEWS.md', 'MIGRATION.md', 'MIGRATION-Hans.md']: origin=ROOT/name
                    if name=='ASSET-SOURCES.md': origin=WEB/'ASSET-SOURCES.md'
                    if name=='documentation-examples.R': origin=ROOT/'tools/documentation-examples.R'
                    if split.path.startswith('../assets/source/tools/teaching-audit/'): origin=ROOT/'tools/teaching-audit'/name
                    if name in {f'prepare-{key}{suffix}.R' for key in RECIPES for suffix in ['', '-zh']}: origin=WEB/'showcase-recipes.json'
                    if name in ['tutorial-inputs.zip','input-hashes.json']: origin=ROOT/'examples'
                    require(origin.exists(), f'{lang}/{p["slug"]}: unknown asset {url}')
    icon = json.loads((WEB/'assets/icon.json').read_text())
    for record in icon['assets']+[icon['favicon']]:
        path=ROOT/record['file']
        require(path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256'],f'Icon hash mismatch: {path}')
    config = (ROOT/'_pkgdown.yml').read_text()
    names = [re.search(r'\\name\{([^}]+)\}',p.read_text())[1] for p in (ROOT/'man').glob('*.Rd')]
    require(all('`'+name+'`' in config for name in names), 'Native pkgdown reference omits Rd topics')
    coverage=json.loads((WEB/'coverage.json').read_text())
    code=snippets()
    require(len({f['id'] for f in coverage})==len(coverage),'Duplicate capability family')
    results=json.loads((WEB/'assets/example-results.json').read_text())
    for family in coverage:
        require(family['result'] is None or family['result'] in results,f'{family["id"]}: missing result record')
        for lang in ['en','zh']:
            prose='\n'.join((WEB/lang/(slug+'.md')).read_text() for slug in family['articles'])
            article_code = {key for marker in re.findall(r'\{\{(?:code|prepare):([a-z,]+)\}\}', prose) for key in marker.split(',')}
            require(all(key in code and key in article_code for key in family['snippets']),f'{family["id"]}: missing shared example in {lang}')
            require(all(key in CAPTIONS and '{{figure:'+key+'}}' in prose for key in family['figures']),f'{family["id"]}: missing output in {lang}')
    # README/vignette resources and links are source paths, not the frozen site.
    for path in [ROOT/'README.md',ROOT/'README-Hans.md',ROOT/'MIGRATION.md',ROOT/'MIGRATION-Hans.md',*(ROOT/'vignettes').glob('*.Rmd')]:
        text = path.read_text()
        ids = set(re.findall(r'\{#([^}]+)\}',text)) | set(re.findall(r'\bid="([^"]+)"',text))
        for url in re.findall(r'\]\(([^)\s]+)\)',text):
            split = urlsplit(url)
            if split.scheme or split.netloc: continue
            if split.path:
                require((path.parent/unquote(split.path)).exists(), f'{path.name}: missing source link {url}')
            elif split.fragment:
                require(unquote(split.fragment) in ids,f'{path.name}: missing chapter {url}')
    if not args.source_only:
        for key, chunks in RECIPES.items():
            for lang in ['en', 'zh']:
                filename = f'prepare-{key}' + ('-zh' if lang == 'zh' else '') + '.R'
                script = ROOT / '.site' / 'assets' / filename
                expected = '\n\n'.join(localize_r_comments(code[name], lang) for name in chunks) + '\n'
                require(script.is_file() and script.read_text() == expected,
                        f'Preparation script differs from shared source: {key} / {lang}')
        out = ROOT/'.site'
        require(out.is_dir(), 'Build .site first')
        docs = {p.resolve():Document(p) for p in out.rglob('*.html')}
        for asset in assets:
            generated=out/'assets'/Path(asset['file']).name
            require(generated.is_file() and hashlib.md5(generated.read_bytes()).hexdigest()==asset['md5'],f'Generated figure differs from source: {asset["id"]}')
        normalize=lambda text: re.sub(r'\s+',' ',text.strip())
        for lang in ['en','zh']:
            for page in PAGES:
                if page['slug']=='reference': continue
                path=(out/lang/(page['slug']+'.html')).resolve()
                expected=blocks(expand((WEB/lang/(page['slug']+'.md')).read_text(),lang,figures='../assets'))
                require(path in docs and list(map(normalize,docs[path].codeblocks))==list(map(normalize,expected)),f'{lang}/{page["slug"]}: rendered code differs from shared source')
        for path,doc in docs.items():
            rel = str(path.relative_to(out))
            require(not [i for i,n in Counter(doc.ids).items() if n>1],f'{rel}: duplicate element IDs')
            if rel.startswith(('en/','zh/')):
                require(doc.lang == ('zh-Hans' if rel.startswith('zh/') else 'en'),f'{rel}: language metadata')
                if Path(rel).name!='reference.html':
                    require('main' in doc.ids,f'{rel}: missing skip target')
                    require('site-language' in doc.ids and 'site-search' in doc.ids,f'{rel}: missing language/search navigation')
                    require({'chapter-nav','chapter-drawer','language-control'} <= doc.classes,
                            f'{rel}: missing site guide or language dropdown')
                    if Path(rel).parent == Path('zh/reference'):
                        for heading in doc.headings:
                            if heading in aliases: continue
                            require(not re.search('[A-Za-z]{3}',heading) or bool(re.search('[\u3400-\u9fff]',heading)),
                                    f'{rel}: untranslated reference heading {heading}')
            for img in doc.images:
                require('alt' in img,f'{rel}: image missing alt')
                if Path(urlsplit(img.get('src','')).path).stem in CAPTIONS:
                    require(int(img.get('width','0'))>0 and int(img.get('height','0'))>0,
                        f'{rel}: figure must reserve its intrinsic dimensions')
            for url in doc.links:
                require('rdrr.io/pkg/ggchord/man/' not in url,f'{rel}: current API autolink escapes local reference {url}')
                split = urlsplit(url)
                if split.scheme or split.netloc: continue
                target = (path.parent/unquote(split.path)).resolve() if split.path else path
                if target.is_dir(): target=target/'index.html'
                require(target.is_relative_to(out),f'{rel}: link escapes site {url}')
                require(target.exists(),f'{rel}: broken local link {url}')
                if split.fragment and target in docs:
                    require(unquote(split.fragment) in docs[target].ids,f'{rel}: broken anchor {url}')
        for lang in ['en','zh']:
            expected = {f'{s}.html' for s in slugs} | {f'reference/{s}.html' for s in rd_topics} | {'reference/index.html','authors.html','LICENSE.html','LICENSE-text.html','404.html'}
            actual = {str(p.relative_to(out/lang)) for p in (out/lang).rglob('*.html')}
            require(actual==expected,f'{lang}: generated page coverage')
            index = json.loads((out/lang/'search.json').read_text())
            require({p['url'] for p in index} == expected-{'reference.html'},f'{lang}: incomplete search coverage')
            require(any('track' in p['text'] for p in index),f'{lang}: empty search content')
        for path,doc in docs.items():
            if path.is_relative_to(out/'en') and path.name!='reference.html':
                partner=(out/'zh'/path.relative_to(out/'en')).resolve()
                require(partner in docs and set(doc.ids)==set(docs[partner].ids),f'{path.relative_to(out)}: counterpart anchors differ')
                if path.parent.name=='reference' and path.name!='index.html' and partner in docs:
                    require([normalize(executable_text(code)) for code in doc.codeblocks]==[normalize(executable_text(code)) for code in docs[partner].codeblocks],f'{path.name}: reference usage/examples differ by language')
                    require(len(re.findall('[\u3400-\u9fff]', ''.join(docs[partner].main_text)))>=10,f'{path.name}: translated reference body is missing')
                    require(bool(doc.codeblocks) or 'usage' not in doc.ids,f'{path.name}: missing usage code')
        require(any('标签' in p['text'] for p in json.loads((out/'zh/search.json').read_text())), 'Han search content absent')
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f'PASS: {len(PAGES)} language pairs, equivalent executable code, assets, {len(exports)} exports / {len(rd_topics)} Rd topics'+(' and generated links/anchors/search.' if not args.source_only else '.'))

if __name__ == '__main__':
    main()
