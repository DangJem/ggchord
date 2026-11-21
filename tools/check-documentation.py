#!/usr/bin/env python3
"""Check bilingual handbook sources and, by default, the local rendered website."""
import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit
from website.reference_manual import load as load_api_reference

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / 'website'
EXAMPLE = re.compile(r'\{\{example:([a-z0-9-]+)\}\}')
ANCHOR = re.compile(r'\{#([^}]+)\}|\bid="([^"]+)"')
LINK = re.compile(r'\]\(([^)\s]+)\)|\b(?:href|src)="([^"]+)"')


def read(path):
    return json.loads(path.read_text())


def source_check():
    subprocess.run([sys.executable, 'tools/sync-website-languages.py', '--check'], cwd=ROOT, check=True)
    structure, _ = load_api_reference()
    pages = read(WEB/'pages.json')
    nav = read(WEB/'navigation.json')
    examples = read(WEB/'handbook-examples.json')
    policy = read(WEB/'content-policy.json')
    errors = []

    def require(ok, message):
        if not ok:
            errors.append(message)

    require((policy.get('source_language'), policy.get('generated_language')) == ('en', 'zh'),
            'The handbook requires an English source and reviewed Han translation.')
    by_slug = {page['slug']: page for page in pages}
    require(len(by_slug) == len(pages), 'Duplicate page slugs')
    canonical = [page for page in pages if page['kind'] != 'alias']
    paths = [page['canonical_path'] for page in canonical]
    require(len(paths) == len(set(paths)), 'Duplicate canonical page paths')
    groups = {group['id'] for group in nav}
    require(len(groups) == len(nav), 'Duplicate navigation groups')
    listed = [key for group in nav for key in group['pages']]
    require(set(listed) | {'index'} == {p['slug'] for p in canonical}
            and len(listed) == len(set(listed)), 'Navigation must cover every handbook article exactly once')
    for group in nav:
        require(group['landing'] in group['pages'] and group.get('parent') in {None, *groups},
                'Unknown navigation landing or parent: ' + group['id'])
        for slug in group['pages']:
            require(by_slug.get(slug, {}).get('module') == group['id'],
                    'Incorrect module for navigation entry: ' + slug)
    sources = {row['source'] for row in policy['markdown']}
    active = {p['slug'] for p in canonical}
    require({f'en/{slug}.md' for slug in active} <= sources, 'Missing canonical Markdown translation entry')
    used = set()
    for page in canonical:
        slug = page['slug']
        texts = [(WEB/lang/(slug+'.md')).read_text() for lang in ('en','zh')]
        markers = [EXAMPLE.findall(t) for t in texts]
        require(markers[0] == markers[1], slug+': example identities or order differ by language')
        used.update(markers[0])
        require(ANCHOR.findall(texts[0]) == ANCHOR.findall(texts[1]), slug+': heading anchors differ')
        for lang, text in zip(('en', 'zh'), texts):
            require(bool(page['title'][lang].strip()) and '中文' not in text, lang+'/'+slug+': localized title or language label')
            path = WEB/lang/(slug+'.md')
            # Source links resolve in the checkout before conversion to website paths.
            prose = re.sub(r'^ *```[^\n]*\n.*?^ *``` *$', '', text, flags=re.M|re.S)
            for matches in LINK.findall(prose):
                url = next(v for v in matches if v)
                parts = urlsplit(url)
                if parts.scheme or parts.netloc:
                    continue
                dest = (path.parent/unquote(parts.path)).resolve() if parts.path else path
                require(dest.is_relative_to(ROOT) and dest.exists(), f'{lang}/{slug}: missing source link {url}')
                if parts.fragment and dest.is_file() and dest.suffix == '.md':
                    ids = {a or b for a,b in ANCHOR.findall(dest.read_text())}
                    require(unquote(parts.fragment) in ids, f'{lang}/{slug}: missing source anchor {url}')
    require(used == set(examples), f'Unknown example markers {sorted(used-set(examples))}, unused examples {sorted(set(examples)-used)}')
    files = {p.stem for p in (WEB/'examples').glob('*.R')}
    require(files == set(examples), 'Example scripts and handbook-examples.json differ: '+str(sorted(files ^ set(examples))))
    for name, spec in examples.items():
        require(bool(re.fullmatch(r'[a-z0-9-]+', name)) and spec['kind'] in ('plot','data'), name+': invalid ID or kind')
        require((WEB/'examples'/(name+'.R')).is_file(), name+': missing script')
        if spec['kind'] == 'plot':
            require(all(spec.get(dim, 0) > 0 for dim in ('width','height')), name+': missing plot dimensions')
    topics = {p.stem: re.findall(r'\\alias\{([^}]+)\}', p.read_text()) for p in (ROOT/'man').glob('*.Rd')}
    aliases = {name for names in topics.values() for name in names}
    exports = set(re.findall(r'^export\(([^)]+)\)', (ROOT/'NAMESPACE').read_text(), re.M))
    require(exports <= aliases, 'Missing Rd for '+str(sorted(exports-aliases)))
    flexible = read(WEB/'flexible-parameters.json')
    for api, params in flexible.items():
        require(api in exports and api in topics, 'Unknown flexible-parameter function: ' + api)
        if api in topics:
            terms = re.findall(r'\\item\{([^}]+)\}', (ROOT/'man'/(api+'.Rd')).read_text())
            documented = {p.strip() for term in terms for p in term.split(',')}
            require(set(params) <= documented, api+': a badge names an undocumented argument')
    require('flex-radius-error' in EXAMPLE.findall((WEB/'en/parameters.md').read_text()) and
            examples.get('flex-radius-error',{}).get('kind') == 'data',
            'The flexible-parameter guide lacks its independently runnable invalid input')
    migrations = read(WEB/'page-migrations.json')
    require(not any('reference/'+name+'.html' in migrations for name in structure['topics']),
            'A current API topic has a retired-page redirect')
    if errors:
        raise SystemExit('\n'.join(dict.fromkeys(errors)))
    print(f'PASS: {len(canonical)} bilingual entries, {len(examples)} example sources, {len(exports)} exports, '
          f'{len(topics)} native help topics and {sum(map(len,flexible.values()))} flexible arguments.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-only', action='store_true')
    parser.add_argument('--reference-only', action='store_true', help='Verify the paired API manual sources')
    parser.add_argument('--lesson-code-snapshot', action='store_true', help='Print registered example fingerprints')
    args = parser.parse_args()
    if args.reference_only:
        structure, _ = load_api_reference()
        print(f'PASS: {len(structure["topics"])} paired API entries')
        return
    if args.lesson_code_snapshot:
        print(json.dumps({name: hashlib.sha256((WEB/'examples'/(name+'.R')).read_bytes()).hexdigest()
                          for name in read(WEB/'handbook-examples.json')}, sort_keys=True))
        return
    source_check()
    if not args.source_only:
        subprocess.run([sys.executable, 'tools/site-audit/check-built.py'], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
