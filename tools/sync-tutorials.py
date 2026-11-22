#!/usr/bin/env python3
"""Generate the paired READMEs and check the separately maintained offline tutorial."""
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / 'website'


def check_offline_tutorial():
    directory = ROOT / 'vignettes'
    source = (directory / 'ggchord_vignette.Rmd').read_bytes()
    english = source.decode('utf-8')
    han = (directory / 'ggchord_guide_hans.Rmd').read_text()
    digest = hashlib.sha256(source).hexdigest()
    if f'<!-- Reviewed English source SHA-256: {digest} -->' not in han:
        raise ValueError('Review the Han offline tutorial against the current English Rmd')
    for name, pattern in [('R chunks', r'```\{r[^}]*\}\n.*?\n```'),
                          ('chapter anchors', r'\{#[\w-]+\}'),
                          ('local assets', r'tutorial-assets/[\w.-]+')]:
        expected = re.findall(pattern, english, re.S)
        if expected != re.findall(pattern, han, re.S):
            raise ValueError(f'Offline tutorial {name} differ between editions')
        if name == 'local assets':
            for asset in set(expected):
                if not (directory / asset).is_file():
                    raise ValueError(f'Missing offline tutorial asset: {asset}')


def readme(lang):
    template = (WEB / 'readme-template.md').read_text()
    translations = json.loads((WEB / 'readme-text.json').read_text())

    def replace(match):
        key = match[1]
        if key == 'design-references':
            fragment = (WEB / f'readme-inspirations-{lang}.md').read_text()
            return fragment.strip().replace('](../inst/NOTICE)', '](inst/NOTICE)')
        return lang if key == 'lang' else translations[key][lang]

    return re.sub(r'\{\{([a-z-]+)\}\}', replace, template)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    subprocess.run(['python3', 'tools/sync-website-languages.py', '--check'], cwd=ROOT, check=True)
    check_offline_tutorial()
    for lang, name in [('en', 'README.md'), ('zh', 'README-Hans.md')]:
        path = ROOT / name
        content = readme(lang)
        if args.check:
            if path.read_text() != content:
                raise SystemExit(f'Out-of-date README: {path}')
        elif content != path.read_text():
            path.write_text(content)
    print('README editions match the shared template; offline tutorial sources verified.')


if __name__ == '__main__':
    main()
