"""Repository-readable Markdown links, converted only at site expansion time."""
import csv
import json
import re
import subprocess
from functools import lru_cache
from pathlib import Path
from urllib.parse import unquote, urlsplit

from . import ROOT
WEB = ROOT / 'website'
INPUT_ARCHIVE = 'https://github.com/DangJem/ggchord/archive/refs/heads/main.zip'
LINK = re.compile(r'(?<=\]\()[^\s)]+(?=\))')
HTML_LINK = re.compile(r'''(\b(?:href|src)=)(['"])(.*?)\2''')


def rewrite_links(text, convert):
    text = LINK.sub(lambda match: convert(match[0]), text)
    return HTML_LINK.sub(lambda match: match[1]+match[2]+convert(match[3])+match[2], text)


@lru_cache(maxsize=1)
def tracked_paths():
    return set(subprocess.check_output(['git','ls-files','-z'], cwd=ROOT, text=True).split('\0'))


@lru_cache(maxsize=1)
def asset_sources():
    rows = list(csv.DictReader((WEB/'assets/plots.csv').open()))
    for name in ('course-figures.json', 'manual-figures.json', 'home-previews.json'):
        rows += json.loads((WEB/name).read_text())
    result = {r['id']+'.png': '../../'+r['file'] for r in rows}
    result['logo.svg'] = '../../man/figures/logo.svg'
    result.update({n: '../../'+n for n in ('NEWS.md','MIGRATION.md','MIGRATION-Hans.md')})
    result['ASSET-SOURCES.md'] = '../ASSET-SOURCES.md'
    result['documentation-examples.R'] = '../../tools/documentation-examples.R'
    for name in ('generate_website_data.R', 'generate_flagship_data.R',
                 'website-data-sources.tsv', 'flagship-data-sources.tsv'):
        result[name] = '../../data-raw/' + name
    return result


def source_url(url):
    part = urlsplit(url)
    if part.scheme or part.netloc or not part.path:
        return url
    path = unquote(part.path)
    fragment = '#'+part.fragment if part.fragment else ''
    if path == 'reference/index.html':
        return '../../man/'
    if path.startswith('reference/') and path.endswith('.html'):
        return '../../man/'+Path(path).stem+'.Rd'
    if path.endswith('.html'):
        return path[:-5]+'.md'+fragment
    if path == '../assets/tutorial-inputs.zip':
        return INPUT_ARCHIVE
    if path.startswith('../assets/source/tools/teaching-audit/'):
        return '../../tools/teaching-audit/'+Path(path).name+fragment
    if path.startswith('../assets/'):
        name = path[len('../assets/'):]
        if name in asset_sources():
            return asset_sources()[name]+fragment
    return url


def source_document(text):
    """Preserve code; add real GitHub anchors without changing semantic IDs."""
    parts = re.split(r'(^```[^\n]*\n.*?^```\s*$)', text, flags=re.M|re.S)
    existing = set(re.findall(r'\bid="([^"]+)"', text))
    for i in range(0,len(parts),2):
        parts[i] = rewrite_links(parts[i], source_url)
        def anchor(match):
            key=match[2]
            if key in existing:return match[0]
            existing.add(key)
            return '<a id="'+key+'"></a>\n\n'+match[0]
        parts[i] = re.sub(r'^(#{1,6} .+?) \{#([^}]+)\}$',anchor,parts[i],flags=re.M)
    return ''.join(parts)


def site_document(text, figures='../assets'):
    """Resolve repository paths for HTML/tutorial output, before adding macros."""
    # Explicit source anchors coexist with Pandoc IDs only in the Markdown.
    text = re.sub(r'<a id="([^"]+)"></a>\n\n(?=#{1,6} [^\n]*\{#\1\})','\n\n',text)
    assets = asset_sources()
    inverse = {}
    # Prefer the source basename where both declared output aliases exist.
    for name,path in assets.items():
        if path not in inverse or name==Path(path).name:inverse[path]=name
    def convert(url):
        if url == '../../':
            return figures+'/ggchord-source.zip'
        if url==INPUT_ARCHIVE:return figures+'/tutorial-inputs.zip'
        part=urlsplit(url)
        if part.scheme or part.netloc or not part.path:return url
        path=unquote(part.path);fragment='#'+part.fragment if part.fragment else ''
        if path in inverse:return figures+'/'+inverse[path]+fragment
        if path=='../../man/':return 'reference/index.html'
        if path.startswith('../../man/') and path.endswith('.Rd'):
            return 'reference/'+Path(path).stem+'.html'+fragment
        if path.startswith('../../tools/teaching-audit/'):
            return figures+'/source/tools/teaching-audit/'+Path(path).name+fragment
        if path.endswith('.md') and not path.startswith('../'):
            return path[:-3]+'.html'+fragment
        return url
    parts=re.split(r'(^```[^\n]*\n.*?^```\s*$)',text,flags=re.M|re.S)
    for i in range(0,len(parts),2):parts[i]=rewrite_links(parts[i],convert)
    return ''.join(parts)


def check_source_links(path):
    """Return broken Markdown/HTML targets and missing explicit source anchors."""
    text=path.read_text()
    text=re.sub(r'^```[^\n]*\n.*?^```\s*$','',text,flags=re.M|re.S)
    urls=LINK.findall(text)+[m[2] for m in HTML_LINK.findall(text)]
    errors=[]
    for url in urls:
        part=urlsplit(url)
        if part.scheme or part.netloc:continue
        target=(path.parent/unquote(part.path)).resolve() if part.path else path
        if not target.exists():
            errors.append(f'{path.relative_to(ROOT)}: missing source target {url}')
            continue
        if target.is_file() and target.is_relative_to(ROOT) and str(target.relative_to(ROOT)) not in tracked_paths():
            errors.append(f'{path.relative_to(ROOT)}: source target is not tracked {url}')
        if part.fragment and target.suffix=='.md':
            anchors=set(re.findall(r'\bid="([^"]+)"',target.read_text()))
            if unquote(part.fragment) not in anchors:
                errors.append(f'{path.relative_to(ROOT)}: missing source anchor {url}')
    return errors


def check_repository_urls(text):
    """Check repository URLs against explicit tags, never an ambiguous short ref."""
    errors = []
    base = 'https://github.com/DangJem/ggchord/'
    def revision(ref):
        if re.fullmatch(r'v\d+\.\d+\.\d+', ref):
            return 'refs/tags/'+ref
        if re.fullmatch(r'[0-9a-f]{7,40}', ref):
            return ref
        return 'refs/heads/'+ref
    def exists(spec):
        candidates = [spec]
        if spec.startswith('refs/heads/'):
            candidates.append(spec.replace('refs/heads/', 'refs/remotes/origin/', 1))
        return any(subprocess.run(['git','cat-file','-e',candidate], cwd=ROOT,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0
                   for candidate in candidates)
    for url in set(LINK.findall(text)):
        if not url.startswith(base):continue
        path = unquote(urlsplit(url).path[len('/DangJem/ggchord/'):])
        if match := re.fullmatch(r'(?:blob|tree)/([^/]+)/(.*)', path):
            ref, file = match.groups()
            valid = (ROOT/file).exists() if ref == 'main' else exists(revision(ref)+':'+file)
            if not valid:errors.append('Missing repository target: '+url)
        elif path.startswith('compare/'):
            for ref in path.removeprefix('compare/').split('...'):
                if not exists(revision(ref)+'^{commit}'):
                    errors.append('Missing comparison revision: '+url)
        elif path.startswith('archive/refs/heads/') and path.endswith('.zip'):
            if not exists('refs/heads/'+path[len('archive/refs/heads/'):-4]+'^{commit}'):
                errors.append('Missing archive branch: '+url)
    return errors
