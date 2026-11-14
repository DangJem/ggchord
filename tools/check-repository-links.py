#!/usr/bin/env python3
"""Check tracked documentation targets/anchors; list external URLs for explicit online review."""
import argparse
import importlib.util
import json
import re
import subprocess
import unicodedata
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
from site_source_links import ROOT, check_repository_urls


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls = []
    def handle_starttag(self, tag, attrs):
        self.urls.extend(value for key,value in attrs if key in ('href','src') and value)


def without_code(text):
    return re.sub(r'^\s*(```|~~~)[^\n]*\n.*?^\s*\1\s*$', '', text, flags=re.M|re.S)


def urls(path):
    text = path.read_text()
    extension = path.suffix.lower()
    found = []
    if extension in ('.md','.rmd'):
        text = without_code(text)
        found += [a or b for a,b in re.findall(
            r'\]\(([^\s)]+)(?:\s+["\x27][^\n]*?["\x27])?\)|^\s*\[[^\]]+\]:\s*(\S+)', text, re.M)]
    if extension in ('.md','.rmd','.html'):
        parser = Links(); parser.feed(text); found += parser.urls
    if extension == '.rd':
        found += re.findall(r'\\(?:url|href)\{([^}]+)\}', text)
    return [value.strip('<>').replace('\\%','%') for value in found]


def anchors(path):
    text = path.read_text()
    result = set(re.findall(r'\b(?:id|name)=["\x27]([^"\x27]+)', text))
    result.update(re.findall(r'\{#([^}]+)\}', text))
    if path.suffix.lower() in ('.md','.rmd'):
        seen = {}
        for heading in re.findall(r'^#{1,6}\s+(.+?)\s*#*$', without_code(text), re.M):
            heading = re.sub(r'<[^>]*>', '', heading)
            heading = re.sub(r'!?\[([^]]*)\]\([^)]*\)', r'\1', heading).replace('`','').lower()
            slug = ''.join(c for c in heading if unicodedata.category(c)[0] not in 'PS' or c in '-_').replace(' ','-')
            number = seen.get(slug,0); seen[slug] = number+1
            result.add(slug+('-'+str(number) if number else ''))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    tracked = set(subprocess.check_output(['git','ls-files','-z'], cwd=ROOT, text=True).split('\0'))
    pages = sorted(ROOT/name for name in tracked if Path(name).suffix.lower() in ('.md','.rmd','.html','.rd'))
    # This shared template is evaluated at the README output location in both
    # languages. Literal {{lang}} is not an independent navigation target.
    spec = importlib.util.spec_from_file_location('tutorials', ROOT/'tools/sync-tutorials.py')
    tutorials = importlib.util.module_from_spec(spec); spec.loader.exec_module(tutorials)
    errors = []
    for lang,name in [('en','README.md'),('zh','README-Hans.md')]:
        if tutorials.readme(lang) != (ROOT/name).read_text():errors.append('Stale generated README: '+name)
    cache = {}; external = {}; count = 0
    for path in pages:
        if path == ROOT/'website/readme-template.md':continue
        links = urls(path); count += len(links)
        for url in links:
            part = urlsplit(url)
            if part.scheme or part.netloc:
                if part.scheme in ('http','https'):
                    external.setdefault(url,[]).append(str(path.relative_to(ROOT)))
                    errors.extend(check_repository_urls('[target]('+url+')'))
                continue
            target = (path.parent/unquote(part.path)).resolve() if part.path else path
            if not target.is_relative_to(ROOT):
                errors.append(f'{path.relative_to(ROOT)}: target escapes repository: {url}')
                continue
            relative = str(target.relative_to(ROOT))
            if not target.exists() or (target.is_file() and relative not in tracked):
                errors.append(f'{path.relative_to(ROOT)}: missing/untracked target: {url}')
            elif part.fragment and target.suffix.lower() in ('.md','.rmd','.html'):
                if target not in cache:cache[target] = anchors(target)
                if unquote(part.fragment) not in cache[target]:
                    errors.append(f'{path.relative_to(ROOT)}: missing anchor: {url}')
    report = dict(files=len(pages),link_occurrences=count,errors=sorted(set(errors)),
                  template='README template validated through both generated outputs',
                  external_urls={u:sorted(set(v)) for u,v in sorted(external.items())},
                  external_scope='Listed for separate online review; HTTP access is not inferred from a static pass.')
    if args.report:args.report.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    if errors:raise SystemExit('\n'.join(sorted(set(errors))))
    print(f'PASS: {len(pages)} tracked documents, {count} links, local targets/anchors and explicit repository revisions; {len(external)} external URLs listed separately.')


if __name__ == '__main__':main()
