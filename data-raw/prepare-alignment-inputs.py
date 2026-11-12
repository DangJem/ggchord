#!/usr/bin/env python3
"""Re-encode upstream CLUSTAL records as FASTA; no alignment is computed.
Run with a directory containing the three unchanged upstream files named below.
The recorded hashes are checked before extraction. Only formatting changes.
"""
import hashlib, json, re, sys
from pathlib import Path
root = Path(__file__).resolve().parents[1]
dest = root / 'examples/alignment'
source = Path(sys.argv[1]) if len(sys.argv) > 1 else dest
manifest = json.loads((dest/'sources.json').read_text())
for item in manifest:
    raw = (source/item['file']).read_bytes()
    if hashlib.sha256(raw).hexdigest() != item['sha256']:
        raise SystemExit('Input hash differs: '+item['file'])
    if 'output' not in item: continue
    rows = {}
    for line in raw.decode().splitlines():
        if not line or line[0].isspace() or line.startswith('CLUSTAL'): continue
        fields = line.split()
        if len(fields)<2 or not re.fullmatch(r'[A-Za-z.*-]+', fields[1]):
            raise SystemExit('Unrecognized source row: '+line)
        rows[fields[0]] = rows.get(fields[0], '') + fields[1]
    if len({len(s) for s in rows.values()}) != 1: raise SystemExit('Unequal alignment lengths')
    (dest/item['output']).write_text(''.join('>'+k+'\n'+s+'\n' for k,s in rows.items()))
    print(item['output'], len(rows), len(next(iter(rows.values()))))
