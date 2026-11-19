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
from urllib.parse import unquote, urljoin, urlsplit
from documentation_lib import ROOT, WEB, PAGES, CAPTIONS, RECIPES, expand, snippets, localize_r_comments
from website.source_links import check_source_links, check_repository_urls, asset_sources
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
        self.home_actions = []
        self.figure_code_targets = []
        self._heading = None
        self._heading_parts = []
        self.navbar_links = []
        self.navbar_text = []
        self._navbar_depth = 0
        self._navbar_link = False
        self.feed(path.read_text())
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == 'div' and (a.get('id') == 'navbar' or self._navbar_depth):
            self._navbar_depth += 1
        if self._navbar_depth and tag == 'a' and a.get('href') and 'nav-link' in a.get('class','').split():
            self.navbar_links.append(a['href'])
            self._navbar_link = True
        if tag == 'a' and {'module-entry','home-source'} & set(a.get('class','').split()):
            self.home_actions.append(a.get('href'))
        if tag == 'a' and a.get('data-code-url'):
            self.figure_code_targets.append(a['data-code-url'])
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
        if self._navbar_link and data.strip(): self.navbar_text.append(data.strip())
        if self._heading is not None: self._heading_parts.append(data)
        if self._in_main: self.main_text.append(data)
        if self._code_parts is not None: self._code_parts.append(data)
    def handle_endtag(self, tag):
        if tag == 'a':
            self._navbar_link = False
        if tag == 'div' and self._navbar_depth: self._navbar_depth -= 1
        if tag == self._heading:
            self.headings.append(''.join(self._heading_parts).strip())
            self._heading = None
        if tag == 'main': self._in_main = False
        if tag == 'pre' and self._code_parts is not None:
            self.codeblocks.append(''.join(self._code_parts))
            self._code_parts = None

class RedirectLinks(HTMLParser):
    """Collect the links a reader can follow when scripts are disabled."""
    def __init__(self, markup):
        super().__init__()
        self.links = {}
        self.refresh = False
        self.active = None
        self.feed(markup)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'a' and attrs.get('href'):
            self.active = (attrs['href'], attrs.get('hreflang'))
            self.links[self.active] = ''
        if tag == 'meta' and attrs.get('http-equiv', '').lower() == 'refresh':
            self.refresh = True

    def handle_data(self, text):
        if self.active:
            self.links[self.active] += text

    def handle_endtag(self, tag):
        if tag == 'a':
            self.active = None


def check_retired_reference_routes(out, docs, require):
    replacements = {
        'geom_gene': 'geom_feature',
        'geom_feature_plasmid': 'geom_feature',
        'geom_feature_label_repel': 'geom_feature_label',
        'position_plasmid': 'position_offset',
        'position_track': 'track_band',
    }
    scripts = {}
    cases = {}
    for old, new in replacements.items():
        for prefix in ('', 'en/', 'zh/'):
            entry = prefix + 'reference/' + old + '.html'
            path = out / entry
            require(path.is_file(), 'Missing retired reference entry: ' + entry)
            if not path.is_file():
                continue
            markup = path.read_text()
            fallback = RedirectLinks(markup)
            require(not fallback.refresh, 'Automatic refresh prevents manual selection: ' + entry)
            actual = {(urljoin(entry, href), tag): text.strip()
                      for (href, tag), text in fallback.links.items()}
            expected_links = {}
            for lang in ([prefix.rstrip('/')] if prefix else ['en', 'zh']):
                tag = 'en' if lang == 'en' else 'zh-Hans'
                expected_links[(lang + '/migration.html#removed-layers', tag)] = (
                    'API migration' if lang == 'en' else '接口迁移说明')
                expected_links[(lang + '/reference/' + new + '.html#ref-examples', tag)] = (
                    'Replacement example' if lang == 'en' else '替代用法示例')
            require(actual == expected_links, 'Incorrect manual migration links or labels: ' + entry)
            script = re.search(r'<script\b[^>]*>(.*?)</script>', markup, re.S)
            require(bool(script), 'Missing redirect script: ' + entry)
            if not script:
                continue
            scripts[entry] = script[1]
            language = prefix or 'en/'
            for fragment, destination in (
                ('', 'migration.html#removed-layers'),
                ('#ref-examples', 'reference/' + new + '.html#ref-examples'),
            ):
                target = language + destination
                page, anchor = target.split('#')
                document = docs.get((out / page).resolve())
                require(document is not None and anchor in document.ids,
                        'Missing migration destination: ' + target)
                cases[entry + fragment] = target
    # Execute only the generated redirect scripts. Capture location.replace and
    # resolve its argument as a URL relative to the page that issued it.
    runner = """
const fs = require('node:fs');
const vm = require('node:vm');
const { scripts, cases } = JSON.parse(fs.readFileSync(0, 'utf8'));
const base = 'https://docs.invalid/ggchord/';
const errors = [];
for (const [entry, expected] of Object.entries(cases)) {
  let current = new URL(entry, base);
  try {
    for (let hop = 0; hop < 3; hop++) {
      const script = scripts[current.href.slice(base.length).split('#')[0]];
      if (!script) break;
      const location = { hash: current.hash, replace: dest => { current = new URL(dest, current); } };
      vm.runInNewContext(script, { location }, { timeout: 500 });
    }
    if (current.href !== new URL(expected, base).href) {
      errors.push(entry + ' -> ' + current.href + ' (expected ' + expected + ')');
    }
  } catch (error) {
    errors.push(entry + ': ' + error.message);
  }
}
process.stdout.write(JSON.stringify(errors));
"""
    errors = json.loads(subprocess.check_output(['node', '-e', runner],
                        input=json.dumps(dict(scripts=scripts, cases=cases)), text=True))
    for error in errors:
        require(False, 'Retired reference redirect: ' + error)


def lesson_code_snapshot():
    result={}
    for path in sorted((WEB/'en').glob('*.md')):
        text=path.read_text()
        for name,keys in RECIPES.items():
            text=text.replace('{{recipe:'+name+'}}','{{prepare:'+','.join(keys)+'}}')
        code=re.findall(r'\{\{(?:code|prepare|output|console):[a-z,]+\}\}|```r\n[\s\S]*?```',text)
        result[str(path.relative_to(ROOT))]=hashlib.sha256(json.dumps(code).encode()).hexdigest()
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-only', action='store_true')
    parser.add_argument('--lesson-code-snapshot',action='store_true')
    args = parser.parse_args()
    if args.lesson_code_snapshot:
        print(json.dumps(lesson_code_snapshot(),sort_keys=True)); return
    subprocess.run(['python3','tools/sync-website-languages.py','--check'],cwd=ROOT,check=True)
    subprocess.run(['python3','tools/check-repository-links.py'],cwd=ROOT,check=True)
    subprocess.run(['python3', str(ROOT / 'tools/check-teaching-content.py')], cwd=ROOT, check=True)
    errors = []
    def require(ok, message):
        if not ok: errors.append(message)
    slugs = [p['slug'] for p in PAGES]
    require(len(slugs) == len(set(slugs)), 'Duplicate page slug')
    navigation = json.loads((WEB/'navigation.json').read_text())
    topics = json.loads((WEB/'manual-sections.json').read_text())
    page_by_slug = {page['slug']: page for page in PAGES}
    primary = [group for group in navigation if 'navbar' in group]
    require([group['navbar'] for group in primary] == ['start','plots','data','gallery','reference'],
            'The five document destinations must appear in the expected header order')
    entries = [slug for group in navigation for slug in group['pages']]
    require(set(entries) == {p['slug'] for p in PAGES if p.get('kind')!='alias'} - {'index'} and len(entries) == len(set(entries)),
            'The topic groups must cover every current article once')
    for group in navigation:
        slug = group['id']
        require(group['landing'] in group['pages'] and
                all(page_by_slug[key]['module'] == slug for key in group['pages']),
                f'{slug}: page categories must match the topic group')
        require(group.get('parent') in {None, *(g['id'] for g in primary)},
                f'{slug}: invalid parent group')
        links = [key for section in topics.get(slug, []) for key in section['pages']]
        require(set(group['pages']) - {group['landing']} <= set(links) and
                set(links) <= set(entries),
                f'{slug}: related topics must cover its pages and link to current articles')
    require(all('prerequisites' not in p and 'next_steps' not in p for p in PAGES),
            'Topic navigation should not prescribe an article order')
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
            errors.extend(check_source_links(WEB/lang/(p['slug']+'.md')))
            errors.extend(check_repository_urls(text))
            require('中文' not in text, f'{lang}/{p["slug"]}: inconsistent language label')
            require(bool(p['title'][lang]), f'{p["slug"]}: missing localized title')
            require('{{' not in expand(text,lang),f'{lang}/{p["slug"]}: unexpanded source marker')
            require(all(re.search(r'\{#[^}]+\}$',heading) for heading in re.findall(r'^#{2,4} .+$',text,re.M)),
                f'{lang}/{p["slug"]}: headings need explicit paired anchors')
        require(re.findall(r'\{#([^}]+)\}', texts[0]) == re.findall(r'\{#([^}]+)\}', texts[1]), f'{p["slug"]}: unmatched anchors')
        require(re.findall(r'\{\{[^}]+\}\}', texts[0]) == re.findall(r'\{\{[^}]+\}\}', texts[1]), f'{p["slug"]}: unmatched shared code/asset markers')
        require(executable_blocks(expand(texts[0], 'en')) == executable_blocks(expand(texts[1], 'zh')), f'{p["slug"]}: mismatched executable code')
    for lang in ['en','zh']:
        actual = {f.stem for f in (WEB/lang).glob('*.md')}
        expected = set(slugs)
        require(actual == expected,
                f'{lang}: orphan source pages: {sorted(actual - expected)}; '
                f'missing source pages: {sorted(expected - actual)}')
    for filename in ['README.md', 'README-Hans.md']:
        errors.extend(check_repository_urls((ROOT/filename).read_text()))
        require(not re.search(r'SnapGene|Geneious', (ROOT/filename).read_text(), re.I),
                filename + ': commercial product name in design references')
    # Only the current release is maintained; historical notes describe their own tags.
    release_version = re.search(r'^Version:\s*(\S+)', (ROOT/'DESCRIPTION').read_text(), re.M)[1]
    release_notes = ROOT/'tools/release-notes'/('v'+release_version+'.md')
    if release_notes.exists():
        errors.extend(check_repository_urls(release_notes.read_text()))
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
    # Public task entries and native reference examples serve different purposes.
    # The rendered clipboard runner checks self-contained plotting and data blocks.
    for entry in inventory:
        name = entry['name']
        review = reviews[name]
        require(all(review.get(key) for key in ('purpose', 'page', 'anchor', 'modes', 'scope')),
                name + ': incomplete function entry')
        page = WEB/'zh'/(review['page']+'.md')
        require(page.is_file() and review['anchor'] in page.read_text(), name + ': unresolved task entry')
    subprocess.run(['python3','tools/check-site-design.py'],cwd=ROOT,check=True)
    subprocess.run(['python3','tools/reference-translations.py','--check'],cwd=ROOT,check=True)
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
    for page in PAGES:
        source = (WEB/'en'/(page['slug']+'.md')).read_text()
        for keys in re.findall(r'\{\{code:([a-z,]+)\}\}', source):
            if keys in {'showcasechord', 'showcasecollinear', 'showcaseplasmid'}: continue
            require(all(key in snippets() and snippets()[key].strip() for key in keys.split(',')),
                    f'{page["slug"]}: empty or missing shared teaching block {keys}')
    subprocess.run(['python3','tools/site-audit/verify-existing-assets.py'],cwd=ROOT,check=True)
    verification=json.loads((ROOT/'tools/site-audit/existing-asset-equivalence.json').read_text())
    current_fingerprints={'assets':verification['expected_dependency_md5']}
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
                    if name in asset_sources(): origin=(WEB/'en'/asset_sources()[name]).resolve()
                    if name=='documentation-examples.R': origin=ROOT/'tools/documentation-examples.R'
                    if split.path.startswith('../assets/source/tools/teaching-audit/'): origin=ROOT/'tools/teaching-audit'/name
                    if name in {f'prepare-{key}{suffix}.R' for key in RECIPES for suffix in ['', '-zh']}: origin=WEB/'showcase-recipes.json'
                    if name in ['tutorial-inputs.zip','input-hashes.json']: origin=ROOT/'examples'
                    require(origin.exists(), f'{lang}/{p["slug"]}: unknown asset {url}')
                elif split.path:
                    require(False, f'{lang}/{p["slug"]}: source link has no site destination {url}')
    icon = json.loads((WEB/'assets/icon.json').read_text())
    for record in icon['assets']+[icon['favicon']]:
        path=ROOT/record['file']
        require(path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256'],f'Icon hash mismatch: {path}')
    config = (ROOT/'_pkgdown.yml').read_text()
    indexed_names = set(re.findall(r"`([^`]+)`", config))
    for rd in (ROOT/'man').glob('*.Rd'):
        names = set(re.findall(r'\\alias\{([^}]+)\}', rd.read_text()))
        require(bool(names & indexed_names), 'Native reference omits ' + rd.stem)
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
            for function in family.get('calls', []):
                require(function+'(' in prose, f'{family["id"]}: missing {function} example in {lang}')
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
        extra_aliases=set(json.loads((WEB/'page-migrations.json').read_text()))-{p['slug']+'.html' for p in PAGES}
        require(out.is_dir(), 'Build .site first')
        docs = {p.resolve():Document(p) for p in out.rglob('*.html')}
        for asset in assets:
            generated=out/'assets'/Path(asset['file']).name
            require(generated.is_file() and hashlib.md5(generated.read_bytes()).hexdigest()==asset['md5'],f'Generated figure differs from source: {asset["id"]}')
        normalize=lambda text: re.sub(r'\s+',' ',text.strip())
        for lang in ['en','zh']:
            for page in PAGES:
                if page['slug']=='reference' or page.get('kind')=='alias': continue
                path=(out/lang/(page['slug']+'.html')).resolve()
                expected=blocks(expand((WEB/lang/(page['slug']+'.md')).read_text(),lang,figures='../assets'))
                require(path in docs and list(map(normalize,docs[path].codeblocks))==list(map(normalize,expected)),f'{lang}/{page["slug"]}: rendered code differs from shared source')
        for path,doc in docs.items():
            rel = str(path.relative_to(out))
            require(not [i for i,n in Counter(doc.ids).items() if n>1],f'{rel}: duplicate element IDs')
            if rel.startswith(('en/','zh/')):
                require(doc.lang == ('zh-Hans' if rel.startswith('zh/') else 'en'),f'{rel}: language metadata')
                alias_paths={p['slug']+'.html' for p in PAGES if p.get('kind')=='alias'}
                if Path(rel).name!='reference.html' and str(Path(rel).relative_to(Path(rel).parts[0])) not in alias_paths|extra_aliases:
                    require('main' in doc.ids,f'{rel}: missing skip target')
                    require('site-language' in doc.ids and 'site-search' in doc.ids,f'{rel}: missing language/search navigation')
                    require('language-control' in doc.classes, f'{rel}: missing language selector')
                    lang, local = rel.split('/',1)
                    key = 'reference' if local.startswith('reference/') else Path(local).stem
                    module = page_by_slug.get(key, {}).get('module', 'information')
                    topic_nav = {'chapter-nav','chapter-drawer'}
                    require((topic_nav <= doc.classes) == bool(topics.get(module)),
                            f'{rel}: related topics must match the page category')
                    require(not {'section-number','lesson-pagination','lesson-stage'} & doc.classes,
                            f'{rel}: article numbering or sequential navigation returned')
                    titles = [group['title'][lang] for group in primary]
                    require([label for label in doc.navbar_text if label in titles] == titles,
                            f'{rel}: missing or misordered primary navigation')
                    hrefs = [urlsplit(url) for url in doc.navbar_links]
                    targets = [str((path.parent / unquote(url.path)).resolve().relative_to(out/lang))
                               for url in hrefs if not url.scheme and not url.netloc and url.path]
                    expected_nav = [page_by_slug[group['landing']]['canonical_path'] for group in primary]
                    require(targets == expected_nav, f'{rel}: primary navigation must use direct destinations')
                    if local == 'index.html':
                        home = json.loads((WEB/'home.json').read_text())['modules']
                        figures = {Path(img.get('src','')).stem for img in doc.images}
                        expected_figures = {key for group in home for key in group['figures']}
                        require(len(home) == 4 and len(expected_figures) == 7 and expected_figures <= figures,
                                f'{rel}: missing homepage group or preview figure')
                        require(doc.home_actions == [item for group in home for item in (group['entry'],group['case'])],
                                f'{rel}: missing try-it or complete-example link')
                        full_figures = {item['id']: item['page'] for item in
                                        json.loads((WEB/'manual-figures.json').read_text())}
                        previews = {item['id']: full_figures[item['full']] for item in
                                    json.loads((WEB/'home-previews.json').read_text())}
                        require(doc.figure_code_targets == [previews[key] for group in home for key in group['figures']],
                                f'{rel}: each preview must lead to the corresponding plotting example')
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
            expected = {f'{s}.html' for s in slugs} | {f'reference/{s}.html' for s in rd_topics} | {'reference/index.html','authors.html','LICENSE.html','LICENSE-text.html','404.html'} | extra_aliases
            actual = {str(p.relative_to(out/lang)) for p in (out/lang).rglob('*.html')}
            require(actual==expected,f'{lang}: generated page coverage')
            index = json.loads((out/lang/'search.json').read_text())
            require({p['url'].split('#')[0] for p in index} == expected-{'reference.html'}-extra_aliases-{p['slug']+'.html' for p in PAGES if p.get('kind')=='alias'},f'{lang}: incomplete canonical search coverage')
            require(any('track' in p['text'] for p in index),f'{lang}: empty search content')
            require({'guide','function','parameter','task','faq'} <= {p.get('kind') for p in index},f'{lang}: missing structured search categories')
            for result in index:
                url=urlsplit(result['url']);target=(out/lang/unquote(url.path)).resolve()
                require(target in docs,f'{lang}: search target missing: '+result['url'])
                if url.fragment: require(unquote(url.fragment) in docs[target].ids,f'{lang}: search section missing: '+result['url'])
        from website.copy_examples import CopiedBlocks
        from website.gallery import GROUPS
        downloads = json.loads((out/'assets/examples.json').read_text())
        gathered = {}
        for lang in ('en', 'zh'):
            gallery = docs[(out/lang/'gallery.html').resolve()]
            require(all('gallery-'+key in gallery.ids for key in GROUPS),
                    lang + ': missing gallery plot family')
            require(all(row['id'] in gallery.ids for row in assets),
                    lang + ': gallery lost an existing figure anchor')
            for page in PAGES:
                url = out/lang/(page['slug']+'.html')
                if not url.exists():
                    continue
                reader = CopiedBlocks()
                reader.feed(url.read_text())
                for index, (_, block) in enumerate(reader.blocks, 1):
                    gathered[(lang, page['slug'], index)] = block
        require(set(gathered) == {(r['language'],r['page'],r['index']) for r in downloads},
                'Every marked example must have a matching R download')
        for example in downloads:
            script = out/example['file']
            require(script.is_file() and hashlib.sha256(script.read_bytes()).hexdigest() == example['sha256'],
                    'Example script hash mismatch: '+example['file'])
            code = gathered.get((example['language'],example['page'],example['index']))
            require(script.is_file() and script.read_text().rstrip('\n') == code,
                    'Downloaded R differs from the code copied on the page: '+example['file'])

        for path,doc in docs.items():
            if path.is_relative_to(out/'en') and path.name!='reference.html':
                partner=(out/'zh'/path.relative_to(out/'en')).resolve()
                require(partner in docs and set(doc.ids)==set(docs[partner].ids),f'{path.relative_to(out)}: counterpart anchors differ')
                if path.parent.name=='reference' and path.stem in rd_topics and partner in docs:
                    require([normalize(executable_text(code)) for code in doc.codeblocks]==[normalize(executable_text(code)) for code in docs[partner].codeblocks],f'{path.name}: reference usage/examples differ by language')
                    require(len(re.findall('[\u3400-\u9fff]', ''.join(docs[partner].main_text)))>=10,f'{path.name}: translated reference body is missing')
                    require(bool(doc.codeblocks) or 'usage' not in doc.ids,f'{path.name}: missing usage code')
        check_retired_reference_routes(out, docs, require)
        require(any('标签' in p['text'] for p in json.loads((out/'zh/search.json').read_text())), 'Han search content absent')
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f'PASS: {len(PAGES)} language pairs, equivalent executable code, assets, {len(exports)} exports / {len(rd_topics)} Rd topics'+(' and generated links/anchors/search.' if not args.source_only else '.'))

if __name__ == '__main__':
    main()
