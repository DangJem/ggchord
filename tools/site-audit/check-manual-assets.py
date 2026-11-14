#!/usr/bin/env python3
"""Check manual image provenance and optional copied build assets; no UI checks."""
import argparse
import hashlib
import json
import re
import subprocess
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
    revised_path=ROOT/'tools/site-audit/home-rendering.json'
    revised=json.loads(revised_path.read_text()) if revised_path.exists() else None
    if revised:
        assert revised['status']=='complete','Homepage image generation incomplete'
        assert revised['source_sha256']['tools/site-flagships.R']==source_sha
        before=subprocess.check_output(['git','show','3a6b412:tools/site-flagships.R'],cwd=ROOT)
        assert hashlib.sha256(before).hexdigest()==render['source_sha256']
        def chunks(text):
            return dict(re.findall(r'^# ---- ([a-z]+) ----\n(.*?)(?=^# ---- [a-z]+ ----)',text,re.M|re.S))
        old_chunks=chunks(before.decode());new_chunks=chunks((ROOT/'tools/site-flagships.R').read_text())
        assert {key for key in old_chunks.keys()|new_chunks.keys() if old_chunks.get(key)!=new_chunks.get(key)}=={'flagshipmulti'}
        assert all(row['original_records']>0 and row['all_source_fields_identical'] for row in revised['source_checks'].values())
        if revised.get('source_check_correction'):
            assert revised['source_check_correction']['verifier_sha256']==sha(ROOT/'tools/site-audit/render-home.R')

    retained = json.loads((ROOT/'tools/long-term-audit/v0.27.0-normal-layout-verification.json').read_text())
    profiles = json.loads((ROOT/'website/reference-catalog.json').read_text())
    complete = [r for r in rows if r.get('replaces')]
    assert len(complete) == len(profiles) == 13
    assert {r['replaces'] for r in complete} == {r['figure'] for r in profiles}
    old_records = {r['retained_file']: r for r in retained['cases'].values()}
    for row in rows:
        file = ROOT/row['file']
        digest = sha(file)
        if revised and row['id'] in ('flagradial','flagsankey'):
            record=revised['figures'][row['id']]
            assert record['status']=='actual PNG rendered'
            assert record['sha256']==digest
            assert all(float(record[k])==float(row[k]) for k in ('width','height','dpi'))
        elif row['id'] in render['figures']:
            record = render['figures'][row['id']]
            if revised:
                assert record['rendered_source_sha256']==render['source_sha256']
                assert all(old_chunks.get(key)==new_chunks.get(key) for key in row['snippets'] if key in old_chunks),row['id']
            else: assert record['rendered_source_sha256'] == source_sha, row['id']
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
    previews=json.loads((ROOT/'website/home-previews.json').read_text()) if revised else []
    for row in previews:
        record=revised['figures'][row['id']]
        assert record['sha256']==sha(ROOT/row['file']),row['id']
        assert record['status']=='actual PNG rendered'
        with Image.open(ROOT/row['file']) as im:
            assert im.size==(round(row['width']*row['dpi']),round(row['height']*row['dpi']))
            assert im.convert('RGB').getpixel((0,0))==(255,252,245)
        if args.built:assert sha(ROOT/'.site/assets'/Path(row['file']).name)==record['sha256']
    if revised:
        home=json.loads((ROOT/'website/home.json').read_text())
        assert [row['slug'] for row in home['modules']]==['multi-sequence','circular','dotplot','fusion']
        assert len(home['modules'][0]['figures'])==4
        assert len(home['capabilities'])==3
        assert home['capabilities'][1]['title']['zh']=='灵活的参数格式'
        assert home['modules'][2]['title']['zh']=='序列比对点图'
        css=(ROOT/'website/assets/site.css').read_text()
        assert 'svh' not in css and '.home-shell{height:' not in css
        assert '.home-logo{display:block;width:224px' in css
        assert 'body[data-page=index] > .container{max-width:1280px}' in css
        assert 'object-fit:cover' not in css
        if args.built:
            for lang in ['en','zh']:
                text=(ROOT/'.site'/lang/'index.html').read_text()
                assert len(re.findall(r'<h1\b',text))==1,lang
                assert 'home-introduction' in text and 'manual-brand' not in text
                assert text.index('home-introduction')<text.index('home-capabilities')
                fragments=re.findall(r'<section class="home-module"[^>]*>(.*?)</section>',text,re.S)
                assert len(fragments)==len(home['modules'])
                for module,fragment in zip(home['modules'],fragments):
                    assert fragment.index('<h2')<fragment.index('home-purpose')<fragment.index('home-module-art')<fragment.index('home-module-actions')
                    assert module['entry'] in fragment
                full_keys={row['full'] for row in previews}
                original=(ROOT/'.site'/lang/'index.html').read_text()
                for row in previews:
                    assert Path(row['file']).name in original
                    target=row.get('full_file',row['full']+'.png')
                    assert f'../assets/{target}' in original
                    assert (ROOT/'.site/assets'/target).exists()
    assert not (ROOT/'website/en/learn.md').exists()
    assert not (ROOT/'website/zh/learn.md').exists()
    result = dict(status='passed', manual_images=len(rows), new_flagship_pngs=len(render['figures']),
                  complete_reference_maps=len(complete), homepage_previews=len(previews),
                  current_flagship_rerenders=2 if revised else 0,
                  retained_original_source_sha256=render['source_sha256'],common_home_background='#FFFCF5',
                  built_assets_checked=args.built, helper_sha256=source_sha,
                  scope='Source and PNG provenance; reused reference maps retain earlier measured-layout evidence. No browser interaction or viewport acceptance.')
    (ROOT/'tools/site-audit/manual-assets-check.json').write_text(json.dumps(result, indent=2)+'\n')
    print(f"PASS: {len(rows)} manual images, {len(previews)} homepage previews, 13 complete reference maps" + (' and built copies' if args.built else ''))


if __name__ == '__main__':
    main()
