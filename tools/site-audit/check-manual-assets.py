#!/usr/bin/env python3
"""Verify current handbook images. Use --record-render LOG after the R renderer exits successfully."""
import argparse
import csv
import hashlib
import json
import re
import struct
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WEB = ROOT / 'website'
RECEIPT = ROOT / 'tools/site-audit/handbook-art-rendering.json'
RENDERER = 'tools/website/render-handbook-art.R'


def read(path):
    return json.loads(path.read_text())


def digest(path, algorithm='sha256'):
    return hashlib.new(algorithm, path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def image_record(row):
    path = ROOT / row['file']
    with path.open('rb') as stream:
        header = stream.read(24)
    require(header[:16] == b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR', f'Invalid PNG: {path}')
    pixels = struct.unpack('>II', header[16:])
    expected = (round(row['width'] * row['dpi']), round(row['height'] * row['dpi']))
    require(pixels == expected, f'{path}: expected {expected}, found {pixels}')
    return dict(file=row['file'], width=row['width'], height=row['height'],
                dpi=row['dpi'], sha256=digest(path))


def render_inputs():
    paths = [ROOT / name for name in (
        RENDERER, 'tools/home-previews.R', 'tools/site-course-examples.R',
        'tools/site-flagships.R', 'tools/documentation-examples.R',
        'website/zh/plasmid-catalog.md', 'website/home-previews.json',
        'website/reference-catalog.json', 'DESCRIPTION')]
    paths.extend(sorted((ROOT / 'R').glob('*.R')))
    paths.extend(sorted((ROOT / 'R').glob('*.rda')))
    paths.extend(sorted((ROOT / 'data').glob('*.rda')))
    return {str(path.relative_to(ROOT)): digest(path) for path in paths}


def verify(args):
    manual = read(WEB / 'manual-figures.json')
    home = read(WEB / 'home-previews.json')
    profiles = read(WEB / 'reference-catalog.json')
    manual_ids = {row['id'] for row in manual}
    home_ids = {row['id'] for row in home}
    require(len(manual_ids) == len(manual) and len(home_ids) == len(home) and
            manual_ids.isdisjoint(home_ids), 'Duplicate manual or homepage image ID')
    catalog = [row for row in manual if row.get('replaces')]
    flagship = [row for row in manual if not row.get('replaces')]
    require(len(catalog) == len(profiles) and
            {row['replaces'] for row in catalog} == {r['figure'] for r in profiles},
            'Reference catalog and full maps differ')
    require(all(row['full'] in manual_ids for row in home), 'Missing homepage full-size map')
    published = read(WEB / 'home.json')['modules']
    require({key for module in published for key in module['figures']} == home_ids,
            'Homepage modules and preview manifest differ')

    records = {r['id']: image_record(r) for r in manual + home}
    if args.record_render:
        # The R command draws every homepage preview before the catalog loop.
        # It prints each accession only after ggsave() finishes that map.
        written = re.findall(r'^Rendered (.+?)\s*$', args.record_render.read_text(), re.M)
        require(written == [row['accver'] for row in profiles],
                'Render log is incomplete or lists different catalog inputs')
    for row in manual:
        require(records[row['id']]['sha256'] == row['sha256'],
                f"{row['id']}: SHA-256 differs from website/manual-figures.json")

    flagship_log = read(ROOT / 'tools/site-audit/flagship-rendering.json')
    with (WEB / 'assets/plots.csv').open() as stream:
        plots = {row['file']: row for row in csv.DictReader(stream)}
    for row in flagship:
        identifier = row['id']
        if row['source'] == 'tools/render-site-flagships.R':
            logged = flagship_log['figures'][identifier]
            require(logged['rendered_source_sha256'] == digest(ROOT / 'tools/site-flagships.R')
                    and all(logged[k] == v for k, v in records[identifier].items()),
                    f'{identifier}: flagship render record differs')
        else:
            registered = plots[row['file']]
            require(registered['md5'] == digest(ROOT / row['file'], 'md5')
                    and registered['object'] == row['object']
                    and all(float(registered[k]) == row[k] for k in ('width', 'height', 'dpi')),
                    f'{identifier}: public plotting manifest differs')

    receipt = dict(renderer=RENDERER, sources=render_inputs(),
                   figures={r['id']: records[r['id']] for r in home + catalog})
    if args.record_render:
        RECEIPT.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + '\n')
    else:
        require(RECEIPT.exists() and read(RECEIPT) == receipt,
                'Handbook render receipt is missing or stale. After a successful render, '
                'run check-manual-assets.py --record-render <completed-render.log>.')
    if args.built:
        assets = args.site_dir / 'assets'
        for row in manual + home:
            name = row['id'] + '.png' if row['id'] in manual_ids else Path(row['file']).name
            require(digest(assets / name) == records[row['id']]['sha256'],
                    f'{assets / name}: copied PNG differs')
    print(f"Verified {len(flagship)} flagship images, {len(catalog)} catalog maps and "
          f"{len(home)} homepage previews" + (' and their built copies.' if args.built else '.'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--built', action='store_true', help='Also verify the final website PNGs.')
    parser.add_argument('--site-dir', type=Path, default=ROOT / '.site')
    parser.add_argument('--record-render', type=Path, metavar='LOG',
                        help='Write the image and source receipt after a successful full handbook render.')
    try:
        verify(parser.parse_args())
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(str(error)) from error
