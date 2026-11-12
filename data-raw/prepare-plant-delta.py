#!/usr/bin/env python3
"""Copy one complete observed plant sequence-pair block; never run an aligner."""
import argparse
import gzip
import hashlib
import json
import urllib.request
from pathlib import Path

root = Path(__file__).resolve().parents[1]
folder = root / 'examples/genome-alignments'
metadata = json.loads((folder / 'arabidopsis-source.json').read_text())
parser = argparse.ArgumentParser()
parser.add_argument('--cache', type=Path, help='Existing upstream compressed file and LICENSE')
args = parser.parse_args()
sources = {}
for key in ('compressed', 'license'):
    item = metadata[key]
    raw = (args.cache / item['upstream_name']).read_bytes() if args.cache else urllib.request.urlopen(item['url']).read()
    if hashlib.sha256(raw).hexdigest() != item['sha256']:
        raise SystemExit('Changed upstream bytes: ' + item['upstream_name'])
    sources[key] = raw
decoded = gzip.decompress(sources['compressed'])
if hashlib.sha256(decoded).hexdigest() != metadata['decompressed_sha256']:
    raise SystemExit('Changed decompressed source')
lines = decoded.splitlines(keepends=True)
headers = [i for i, line in enumerate(lines) if line.startswith(b'>')]
if headers[0] != 2 or lines[headers[0]].decode().strip() != metadata['selected_header']:
    raise SystemExit('Changed first complete sequence-pair block')
selected = b''.join(lines[:headers[1]])
if hashlib.sha256(selected).hexdigest() != metadata['product_sha256']:
    raise SystemExit('Changed selected prefix')
(folder / metadata['file']).write_bytes(selected)
(folder / 'ASSEMBLY-EXAMPLE-LICENSE.txt').write_bytes(sources['license'])
print('Verified and copied unchanged first sequence-pair block:', len(selected), 'bytes')
