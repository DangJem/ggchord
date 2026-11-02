#!/usr/bin/env python3
"""Conservative plot dependencies for reviewed, independent showcase chunks.

Generation-file hashes remain historical provenance. Dependency hashes establish
whether the current engine and executable preparation are still equivalent.
All chunks are dependencies except the two explicitly independent plot families.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def fingerprints():
    source = (ROOT / 'tools/documentation-examples.R').read_text()
    markers = list(re.finditer(r'^# ---- ([a-z]+) ----\n', source, re.M))
    chunks = {m[1]: source[m.end():markers[i + 1].start()].removesuffix('\n')
              for i, m in enumerate(markers[:-1])}
    # These chunks create only their own plots/tables. No other plot consumes
    # their outputs. Shared data, themes, imports and all R code stay included.
    comparison = next(k for k, v in chunks.items() if 'label_plots <-' in v)
    engine = '\n'.join(str(p.relative_to(ROOT)) + ':' +
                       hashlib.md5(p.read_bytes()).hexdigest()
                       for p in sorted((ROOT / 'R').glob('*.R')))
    captions = json.loads((ROOT / 'website/assets/captions.json').read_text())
    result = {}
    for key in captions:
        excluded = set()
        if key != 'showcasecollinear':
            excluded.add('showcasecollinear')
        if key not in {'labelauto', 'labelradial', 'labelarc'}:
            excluded.add(comparison)
        payload = engine + '\n' + '\n'.join(
            name + ':\n' + value for name, value in chunks.items() if name not in excluded)
        result[key] = hashlib.md5(payload.encode()).hexdigest()
    return {'code': hashlib.md5(source.encode()).hexdigest(), 'assets': result}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stamp-current', action='store_true')
    args = parser.parse_args()
    result = fingerprints()
    if args.stamp_current:
        path = ROOT / 'website/assets/plots.csv'
        with path.open() as handle:
            reader = csv.DictReader(handle)
            fields, rows = reader.fieldnames, list(reader)
        # This migration is permitted only immediately after a complete current
        # export. Never infer equivalence from old PNGs or silently refresh code hashes.
        for row in rows:
            assert row['code_md5'] == result['code'], 'Full current export required'
            assert hashlib.md5((ROOT / row['file']).read_bytes()).hexdigest() == row['md5']
            row['dependency_md5'] = result['assets'][row['id']]
        if 'dependency_md5' not in fields:
            fields.append('dependency_md5')
        with path.open('w') as handle:
            writer = csv.DictWriter(handle, fields, lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        print(f'Stamped reviewed dependencies for {len(rows)} current PNG exports.')
    else:
        print(json.dumps(result))


if __name__ == '__main__':
    main()
