#!/usr/bin/env python3
"""Prepare shared offline reading assets; preserve all canonical plot PNGs."""
import hashlib
import json
import re
import shutil
from pathlib import Path
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / 'vignettes/tutorial-assets'
out.mkdir(exist_ok=True)
shutil.copy2(ROOT / 'man/figures/logo.svg', out / 'logo.svg')
records = []
refs = set()
for path in (ROOT / 'vignettes').glob('*.Rmd'):
    refs.update(re.findall(r'tutorial-assets/([a-z]+\.(?:png|jpg))', path.read_text()))
for name in sorted(refs):
    source = ROOT / 'man/figures' / (Path(name).stem + '.png')
    target = out / name
    with Image.open(source) as original:
        picture = original.convert('RGBA')
        picture.thumbnail((900, 900) if name != 'logo.png' else (320, 320))
        if target.suffix == '.jpg':
            background = Image.new('RGB', picture.size, '#fffcf5')
            background.paste(picture, mask=picture.getchannel('A'))
            background.save(target, quality=80, optimize=True)
        else:
            picture.save(target, optimize=True)
    records.append(dict(file=name, source=str(source.relative_to(ROOT)),
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        sha256=hashlib.sha256(target.read_bytes()).hexdigest(), bytes=target.stat().st_size))
(ROOT / 'website/assets/vignette-figures.json').write_text(json.dumps(records, indent=2)+'\n')
print(f'{len(records)} offline figures, {sum(r["bytes"] for r in records):,} bytes; canonical PNGs unchanged.')
