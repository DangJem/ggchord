#!/usr/bin/env python3
"""Reproduce pinned plant relation subsets; copy complete original byte lines.

Usage: python3 data-raw/prepare-feature-relations.py [--cache DIRECTORY]
The optional cache contains local_name entries from the manifest. Otherwise
fetch only the pinned URLs into a temporary directory. No analysis is run.
The four-column tool position tables are not standard GFF3 or complete source
length dictionaries. Coordinate interpretation must be declared by the caller.
"""
import argparse
import hashlib
import json
import tempfile
import urllib.request
from pathlib import Path


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path)
    args = parser.parse_args()
    folder = Path(__file__).resolve().parents[1] / 'inst/extdata'
    manifest = json.loads((folder / 'feature-relation-sources.json').read_text())
    inputs = {}
    with tempfile.TemporaryDirectory(prefix='ggchord-relations-') as temporary:
        cache = args.cache or Path(temporary)
        for item in manifest['inputs']:
            path = cache / item['local_name']
            if args.cache is None:
                urllib.request.urlretrieve(item['url'], path)
            raw = path.read_bytes()
            if digest(raw) != item['sha256']:
                raise SystemExit('Source hash mismatch: ' + item['role'])
            inputs[item['role']] = raw.splitlines(keepends=True)

    pairs = inputs['pairs'][:20]
    pair_ids = {token for line in pairs for token in line.split()[:2]}
    group_indices = [i for i, line in enumerate(inputs['groups'], 1)
                     if 2 <= len(line.split()) - 1 <= 4][:3]
    groups = [inputs['groups'][i - 1] for i in group_indices]
    group_ids = {token for line in groups for token in line.split()[1:]}

    def positions(role, ids):
        indices = [i for i, line in enumerate(inputs[role], 1) if line.split()[1] in ids]
        selected = [inputs[role][i - 1] for i in indices]
        found = [line.split()[1] for line in selected]
        if set(found) != ids or len(found) != len(ids):
            raise SystemExit('Missing or ambiguous original feature IDs in ' + role)
        return b''.join(selected), indices

    products = {
        'plant-feature-pairs.homology': (b''.join(pairs), list(range(1, 21))),
        'plant-pair-features.gff': positions('pair_features', pair_ids),
        'plant-feature-groups.txt': (b''.join(groups), group_indices),
        'plant-group-features.gff': positions('group_features', group_ids),
        'FEATURE-RELATIONS-LICENSE.txt': (b''.join(inputs['license']), None),
    }
    # Verify every product and its original physical line index before writing.
    for item in manifest['outputs']:
        raw, indices = products[item['file']]
        if digest(raw) != item['sha256'] or indices != item['source_rows']:
            raise SystemExit('Subset rule/hash mismatch: ' + item['file'])
    for name, (raw, indices) in products.items():
        (folder / name).write_bytes(raw)
        print(name, len(indices or []), 'original lines; hash verified')


if __name__ == '__main__':
    main()
