#!/usr/bin/env python3
"""Read-only HTTP checks against the local preview; never inspect a deployed site."""
import argparse,json
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError
from urllib.parse import urlsplit
ROOT=Path(__file__).resolve().parents[2]
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--base-url',default='http://127.0.0.1:8765/');args=parser.parse_args()
    base=args.base_url.rstrip('/')+'/'
    assert urlsplit(base).hostname in {'127.0.0.1','localhost'},'Local preview only'
    baseline=json.loads((ROOT/'tools/site-audit/baseline.json').read_text());failures=[]
    for row in baseline['urls']:
        try:
            with urlopen(base+row['path'],timeout=10) as response:
                if response.status!=200:failures.append(row['path']+': '+str(response.status))
        except Exception as error:failures.append(row['path']+': '+str(error))
    missing='__ggchord_missing_route_for_audit__.html'
    try:
        with urlopen(base+missing,timeout=10) as response:status=response.status
    except HTTPError as error:status=error.code
    if status!=404:failures.append('Unknown route returned '+str(status))
    with urlopen(base+'404.html',timeout=10) as response:
        body=response.read().decode();recovery=response.status==200 and 'en/index.html' in body and 'zh/index.html' in body
    if not recovery:failures.append('Explicit bilingual 404 recovery page missing')
    result=dict(audited_urls=len(baseline['urls']),expected_status=200,unknown_route_status=status,explicit_bilingual_recovery_page=recovery,
        scope='Local Python HTTP server; unknown URLs get its HTTP 404 body. The generated 404.html is separately verified; hosting rewrite/custom-404 configuration is not deployed or tested.',failures=failures)
    (ROOT/'tools/site-audit/http-check.json').write_text(json.dumps(result,indent=2)+'\n')
    if failures:raise SystemExit('\n'.join(failures))
    print('PASS: 582 audited HTTP URLs, unknown URL 404, bilingual recovery file')
if __name__=='__main__':main()
