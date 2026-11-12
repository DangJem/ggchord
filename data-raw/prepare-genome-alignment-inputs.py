#!/usr/bin/env python3
"""Verify unchanged pinned examples; --download explicitly retrieves originals."""
import argparse, hashlib, json, urllib.request
from pathlib import Path
root = Path(__file__).resolve().parents[1]
folder = root/'examples/genome-alignments'
parser = argparse.ArgumentParser()
parser.add_argument('--download', action='store_true')
args = parser.parse_args()
for item in json.loads((folder/'sources.json').read_text()):
    if item.get('status') == 'historical-only':
        print('Historical original retained; excluded from active regeneration:', item['file'])
        continue
    path = folder/item['file']
    if args.download:
        with urllib.request.urlopen(item['url']) as response:
            raw = response.read()
        if hashlib.sha256(raw).hexdigest() != item['sha256']:
            raise SystemExit('Downloaded bytes changed: '+item['file'])
        path.write_bytes(raw)
    if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
        raise SystemExit('Source bytes changed: '+item['file'])
    print('Verified unchanged input:', item['file'])
