#!/usr/bin/env python3
"""Retrieve two pinned upstream originals for the optional graph prototype.

Usage: python3 data-raw/prepare-graph-prototype.py /tmp/chosen-cache
No analysis, assembly or graph construction is performed. Runtime does not
download inputs. The caller explicitly chooses this development-only retrieval.
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

cache = Path(sys.argv[1])
cache.mkdir(parents=True, exist_ok=True)
base = 'https://raw.githubusercontent.com/pangenome/odgi/3f59d596cff3bbe9ddff3862e3212812254496e5/'
inputs = [
    ('test/LPA.gfa', 'LPA.gfa', '9017b433f35b604bdcafd9339318f1649585bceda4ccf263f1ee1f16adf0cdf3'),
    ('LICENSE', 'LICENSE', 'e8e719c7e069ef90fb30751232f628ed0e98293150049a72a0bba544a4c6fdd4'),
]
downloaded = []
for relative, filename, expected in inputs:
    with urllib.request.urlopen(base + relative) as response:
        raw = response.read()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise SystemExit('Pinned source bytes changed: ' + relative)
    downloaded.append((filename, raw))
for filename, raw in downloaded:
    target = cache / filename
    if target.exists() and target.read_bytes() != raw:
        raise SystemExit('Refusing to replace different cache file: ' + str(target))
    target.write_bytes(raw)
    print(target, 'original bytes verified')
