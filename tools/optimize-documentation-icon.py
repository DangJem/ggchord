#!/usr/bin/env python3
"""Losslessly encode supplied artwork; never redraw or quantize it.
Pillow and cwebp are optional asset-maintenance tools, not package/site build deps.
Usage: python3 tools/optimize-documentation-icon.py /absolute/input.png
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from PIL import Image, PngImagePlugin
ROOT = Path(__file__).resolve().parent.parent
source = Path(sys.argv[1])
source_data = source.read_bytes()
original = Image.open(source)
original.load()
assert original.mode == 'RGBA', 'Expected the supplied RGBA PNG'
pixels = original.tobytes()
output = ROOT/'man/figures/logo.png'
metadata = PngImagePlugin.PngInfo()
for key, value in original.info.items():
    if isinstance(value, str):
        metadata.add_text(key, value)
options = dict(optimize=True, compress_level=9, pnginfo=metadata)
for key in ['icc_profile','dpi','exif']:
    if key in original.info:
        options[key] = original.info[key]
original.save(output, **options)
# Keep original encoding if the PNG optimizer cannot make it smaller.
if output.stat().st_size >= len(source_data) and source.resolve() != output.resolve():
    output.write_bytes(source_data)
webp = ROOT/'website/assets/logo.webp'
subprocess.run(['cwebp','-lossless','-exact','-z','9','-metadata','all',str(output),'-o',str(webp)],check=True)
for path in [output,webp]:
    decoded = Image.open(path).convert('RGBA')
    assert decoded.size == original.size and decoded.tobytes() == pixels, f'Pixel mismatch: {path}'
# The favicon is an explicitly documented thumbnail; full artwork is never resized.
favicon = ROOT/'website/assets/favicon.png'
thumb = original.copy()
thumb.thumbnail((64,64), Image.Resampling.LANCZOS)
thumb.save(favicon,optimize=True)
record = dict(provided='User-supplied ggchord artwork, 2026-11-01',
    source_filename=source.name, source_sha256=hashlib.sha256(source_data).hexdigest(),
    source_bytes=len(source_data), width=original.width,height=original.height,
    mode=original.mode, rgba_sha256=hashlib.sha256(pixels).hexdigest(),
    encoding='PNG optimization and exact lossless WebP; no crop, palette reduction or resampling of full artwork',
    permission='User authorized use in ggchord website and documentation; no third-party license inferred',
    assets=[dict(file=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
        reduction_percent=round(100*(1-p.stat().st_size/len(source_data)),2)) for p in [output,webp]],
    favicon=dict(file=str(favicon.relative_to(ROOT)),size=[64,64],purpose='navigation thumbnail only',
        sha256=hashlib.sha256(favicon.read_bytes()).hexdigest()))
(ROOT/'website/assets/icon.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(record,ensure_ascii=False,indent=2))
