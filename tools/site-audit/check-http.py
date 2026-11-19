#!/usr/bin/env python3
"""Check registered pages, historical addresses and downloads on a running local preview."""
import argparse
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlsplit
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-url', default='http://127.0.0.1:8765/')
    parser.add_argument('--site', type=Path, default=ROOT/'.site')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    base = args.base_url.rstrip('/')+'/'
    if urlsplit(base).hostname not in {'127.0.0.1','localhost','::1'}:
        parser.error('A loopback preview URL is required')
    site = args.site.resolve()
    urls = {p.relative_to(site).as_posix() for p in site.rglob('*.html')}
    urls.add('assets/ggchord-source.zip')
    urls.update(read(ROOT/'website/historical-pages.json'))
    urls.update(record['file'] for record in read(site/'assets/examples.json'))
    urls.update('assets/handbook/'+key+'.png' for key in read(site/'assets/handbook/rendered.json'))
    errors = []
    for path in sorted(urls):
        try:
            with urlopen(base+quote(path,safe='/'),timeout=10) as response:
                if response.status!=200:
                    errors.append(f'{path}: HTTP {response.status}')
        except (OSError,URLError) as error:
            errors.append(f'{path}: {error}')
    try:
        with urlopen(base+'__ggchord_missing_route_for_audit__.html',timeout=10) as response:
            status=response.status
    except HTTPError as error:
        status=error.code
    if status!=404:
        errors.append(f'Unknown route returned {status}')
    with urlopen(base+'404.html',timeout=10) as response:
        body=response.read().decode()
    if not all(lang+'/index.html' in body for lang in ('en','zh')):
        errors.append('Missing bilingual recovery links on 404.html')
    if args.report:
        args.report.write_text(json.dumps(dict(checked_urls=len(urls),missing_page_status=status,errors=errors),indent=2)+'\n')
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f'PASS: {len(urls)} local HTTP URLs, expected 404 and bilingual recovery links.')


if __name__=='__main__':
    main()
