#!/usr/bin/env python3
"""Check manual image provenance and optional copied build assets; no UI checks."""
import argparse
import hashlib
import json
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--built', action='store_true')
    args = parser.parse_args()
    rows = json.loads((ROOT/'website/manual-figures.json').read_text())
    render = json.loads((ROOT/'tools/site-audit/flagship-rendering.json').read_text())
    source_sha = sha(ROOT/'tools/site-flagships.R')
    retained = json.loads((ROOT/'tools/long-term-audit/v0.27.0-normal-layout-verification.json').read_text())
    profiles = json.loads((ROOT/'website/reference-catalog.json').read_text())
    complete = [r for r in rows if r.get('replaces')]
    assert len(complete) == len(profiles) == 13
    assert {r['replaces'] for r in complete} == {r['figure'] for r in profiles}
    old_records = {r['retained_file']: r for r in retained['cases'].values()}
    for row in rows:
        file = ROOT/row['file']
        digest = sha(file)
        if row['id'] in render['figures']:
            record = render['figures'][row['id']]
            assert record['rendered_source_sha256'] == source_sha, row['id']
            assert record['sha256'] == digest, row['id']
        else:
            assert row['sha256'] == digest, row['id']
        if row.get('replaces'):
            record = old_records[row['file']]
            assert record['complete_inputs_identical']
            assert record['source_labels'] == row['source_labels'], row['id']
            assert record['png_sha256'] == digest, row['id']
        else:
            # All seven homepage images use the retained flagship paper colour.
            with Image.open(file) as im:
                assert im.convert('RGB').getpixel((0, 0)) == (255, 252, 245), row['id']
        if args.built:
            assert sha(ROOT/'.site/assets'/f"{row['id']}.png") == digest, row['id']
    assert not (ROOT/'website/en/learn.md').exists()
    assert not (ROOT/'website/zh/learn.md').exists()
    result = dict(status='passed', manual_images=len(rows), new_flagship_pngs=len(render['figures']),
                  complete_reference_maps=len(complete), common_home_background='#FFFCF5',
                  built_assets_checked=args.built, helper_sha256=source_sha,
                  scope='Source and PNG provenance; reused reference maps retain earlier measured-layout evidence. No browser interaction or viewport acceptance.')
    (ROOT/'tools/site-audit/manual-assets-check.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f"PASS: {len(rows)} manual images, four new flagships, 13 complete reference maps" + (' and built copies' if args.built else ''))


if __name__ == '__main__':
    main()
