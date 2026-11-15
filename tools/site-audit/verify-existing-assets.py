#!/usr/bin/env python3
"""Verify current website PNGs against their actual generator dependency receipt.

Run tools/documentation-plots.R after source changes. Historical render and
scoped equivalence records remain unchanged and are not current certificates.
This checks normal website PNGs only, never a package/release device matrix.
"""
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('fingerprints', ROOT/'tools/documentation-fingerprints.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def main():
    fingerprints = module.fingerprints()
    rows = list(csv.DictReader((ROOT/'website/assets/plots.csv').open()))
    assert len(rows) == len({r['id'] for r in rows}), 'Duplicate image IDs'
    assert {r['id'] for r in rows} == set(fingerprints['assets']), 'Image inventory changed'
    for row in rows:
        name = row['id']
        assert row['dependency_md5'] == fingerprints['assets'][name], name + ': stale render dependencies; rerun tools/documentation-plots.R'
        assert row['code_md5'] == fingerprints['code'], name + ': stale generator receipt'
        path = ROOT/row['file']
        assert hashlib.md5(path.read_bytes()).hexdigest() == row['md5'], name + ': PNG bytes changed'
    result = dict(status='passed current PNG dependency and byte checks',
                  generator='tools/documentation-plots.R', figure_count=len(rows),
                  expected_dependency_md5=fingerprints['assets'],
                  figures={r['id']:r['md5'] for r in rows},
                  scope='Normal website PNG generation only; historical receipts are preserved separately')
    (ROOT/'tools/site-audit/existing-asset-equivalence.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f'PASS: {len(rows)} website PNGs match current source dependencies and recorded bytes')


if __name__ == '__main__':
    main()
