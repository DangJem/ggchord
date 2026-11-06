#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check the controlled design palette and record scoped contrast measurements.

Uses only Python's standard library. --write refreshes contrast.json after a
reviewed CSS change. This is a token/static check, not an accessibility audit or
WCAG conformance certification; browser computation and keyboard tests remain
separate evidence.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CSS = ROOT / 'website/assets/site.css'
REPORT = ROOT / 'tools/visual-audit/contrast.json'
NORMAL_TEXT = 4.5
LARGE_TEXT = 3.0
NECESSARY_UI = 3.0
SURFACES = ('--surface-paper', '--surface-canvas', '--surface-soft', '--surface-code')
TEXT = ('--text-main', '--text-strong', '--text-muted', '--color-brand',
        '--color-chord', '--color-collinear', '--color-circular')


def tokens(css):
    result = {}
    for name, selector in [('light', ':root'), ('dark', '[data-bs-theme=dark]')]:
        match = re.search(re.escape(selector) + r'\s*\{([^}]+)\}', css)
        if not match:
            raise ValueError(f'Missing {selector} token block')
        result[name] = dict(re.findall(r'(--[\w-]+)\s*:\s*([^;{}]+)\s*;', match[1]))
    return result


def color(values, token, trail=()):
    if token in trail or token not in values:
        raise ValueError(f'Unresolved or cyclic colour token: {token}')
    value = values[token].strip()
    variable = re.fullmatch(r'var\((--[\w-]+)\)', value)
    if variable:
        return color(values, variable[1], trail + (token,))
    if not re.fullmatch(r'#[0-9a-fA-F]{6}', value):
        raise ValueError(f'{token}: contrast colours must be opaque 6-digit sRGB values')
    return value.lower()


def luminance(hex_value):
    channels = [int(hex_value[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    linear = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in channels]
    return sum(value * weight for value, weight in zip(linear, (0.2126, 0.7152, 0.0722)))


def contrast(foreground, background):
    low, high = sorted((luminance(foreground), luminance(background)))
    return (high + 0.05) / (low + 0.05)


def measurements(css):
    themes = tokens(css)
    rows = []
    for theme, values in themes.items():
        def add(role, foreground, background, threshold):
            fg, bg = color(values, foreground), color(values, background)
            ratio = contrast(fg, bg)
            rows.append(dict(theme=theme, role=role, foreground_token=foreground,
                             background_token=background, foreground=fg, background=bg,
                             ratio=round(ratio, 4), minimum=threshold, passed=ratio >= threshold))
        for foreground in TEXT:
            for background in SURFACES:
                add('normal text / labelled module accents', foreground, background, NORMAL_TEXT)
        for background in SURFACES:
            add('custom focus indicator / necessary non-text', '--color-focus', background, NECESSARY_UI)
        add('primary button text', '--surface-paper', '--color-brand', NORMAL_TEXT)
        add('primary button hover text', '--surface-paper', '--text-strong', NORMAL_TEXT)
        # Decorative separators intentionally stay quiet. Necessary control
        # boundaries must be a separate token, or fall back to the actual rule.
        control = '--color-control' if '--color-control' in values else '--color-rule'
        for background in SURFACES:
            add('necessary control boundary', control, background, NECESSARY_UI)
    return themes, rows


def report():
    css = CSS.read_text()
    themes, rows = measurements(css)
    failures = [r for r in rows if not r['passed']]
    static = dict(
        no_scientific_image_inversion=not bool(re.search(r'filter\s*:[^;}]*invert\s*\(', css, re.I)),
        reduced_motion_rule='prefers-reduced-motion:reduce' in css.replace(' ', ''),
        visible_focus_rule=':focus-visible' in css,
        no_external_font_or_stylesheet_import=not bool(re.search(r'@import|@font-face|https?://', css, re.I)),
        reading_width_separate_from_wide_content='--content-reading' in css and '--content-wide' in css,
        legacy_palette_removed=not bool(re.search(r'--(?:accent|ink|body-ink|paper|canvas|gg-soft|gg-muted)\s*:', css)))
    return dict(schema=1, source='website/assets/site.css', source_sha256=hashlib.sha256(css.encode()).hexdigest(),
                method=dict(
                    reference='https://www.w3.org/TR/WCAG22/#contrast-minimum',
                    non_text_reference='https://www.w3.org/TR/WCAG22/#non-text-contrast',
                    srgb='c = channel / 255; c_linear = c / 12.92 for c <= 0.04045, otherwise ((c + 0.055) / 1.055)^2.4',
                    luminance='L = 0.2126 R_linear + 0.7152 G_linear + 0.0722 B_linear',
                    ratio='(max(L1,L2) + 0.05) / (min(L1,L2) + 0.05); compare the unrounded value',
                    thresholds=dict(normal_text=NORMAL_TEXT, large_text=LARGE_TEXT, necessary_ui=NECESSARY_UI),
                    large_text='At least 18pt (24 CSS px), or 14pt (18.667 CSS px) bold. This audit conservatively applies 4.5 even to heading/accent tokens.',
                    focus='Custom focus colour is tested at 3:1 against adjacent controlled surfaces; this does not claim AAA Focus Appearance conformance.',
                    decorative_rules='Quiet separators are not scored as necessary controls; input/control boundaries use --color-control when declared.'),
                tokens=themes, measurements=rows, static_checks=static,
                passed=not failures and all(static.values()),
                minimum_text_ratio=min(r['ratio'] for r in rows if r['minimum'] == NORMAL_TEXT),
                minimum_necessary_ui_ratio=min(r['ratio'] for r in rows if r['minimum'] == NECESSARY_UI),
                limits=['Controlled opaque token colours only; no claim of complete WCAG 2.2 AA certification',
                        'Browser overrides, native syntax-highlight colours, image text, states and focus occlusion require browser evidence',
                        'The tokens must actually be applied to their intended components; this script does not implement a CSS cascade engine',
                        '44 CSS px is a preferred interaction target, not a blanket WCAG AA minimum; 2.5.8 generally uses 24x24 with exceptions'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    data = report()
    serialized = json.dumps(data, ensure_ascii=False, indent=2) + '\n'
    if args.write:
        REPORT.parent.mkdir(exist_ok=True)
        REPORT.write_text(serialized)
    elif not REPORT.exists() or REPORT.read_text() != serialized:
        raise SystemExit('Contrast evidence is stale: python3 tools/check-site-design.py --write')
    failures = [r for r in data['measurements'] if not r['passed']]
    for row in failures:
        print(f"FAIL: {row['theme']} {row['role']}: {row['foreground_token']} / {row['background_token']} = {row['ratio']}:1 < {row['minimum']}:1")
    for key, passed in data['static_checks'].items():
        if not passed:
            print('FAIL static:', key)
    if not data['passed']:
        raise SystemExit(1)
    print(f"PASS: {len(data['measurements'])} controlled token contrast pairs and {len(data['static_checks'])} scoped static checks; not a WCAG certification")


if __name__ == '__main__':
    main()
