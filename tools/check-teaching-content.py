#!/usr/bin/env python3
"""Contextual teaching lint and paired-page inventory; not a browser certificate."""
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGES = json.loads((ROOT / 'website/pages.json').read_text())
PATTERNS = {
    'developer narration': r'本网站面向开发版|已实测的版本差异|尚未推送|本地 checkout|不要假定远端 main|源码 SHA|not yet pushed|tested version differences|verified differences|development version|local checkout|source SHA',
    'internal revision': r'\b[0-9a-f]{40}\b',
    'manual HTML break': r'<br\s*/?>',
    'classroom wording': r'课程|教程|本课|首课|食谱|学习|\b(?:lessons?|tutorials?|exercises?|teaching)\b',
}
series="v"+".".join(re.search(r"^Version: (.+)$",(ROOT/"DESCRIPTION").read_text(),re.M)[1].split(".")[:3])
rows, findings = [], []
for page in PAGES:
    for language in ('en', 'zh'):
        path = ROOT / 'website' / language / (page['slug'] + '.md')
        text = path.read_text()
        fenced = False
        hits = []
        table_columns = None
        for number, line in enumerate(text.splitlines(), 1):
            if re.match(r'^\s*(```|~~~)', line):
                fenced = not fenced
                continue
            if fenced:
                continue
            if line.startswith('|'):
                columns = len(re.split(r'(?<!\\)\|', line)) - 2
                if table_columns is not None and columns != table_columns:
                    raise SystemExit(f'Table column mismatch: {path}:{number}')
                table_columns = columns
            else:
                table_columns = None
            visible = re.sub(r'\{\{.*?\}\}|\{#[^}]+\}|<[^>]+>|`[^`]+`|\]\([^)]*\)', '', line)
            for category, pattern in PATTERNS.items():
                if re.search(pattern, visible if category == 'classroom wording' else line, re.I):
                    # Installation/migration explain versions intentionally.
                    exempt = page['slug'] in ('about', 'migration')
                    findings.append(dict(path=str(path.relative_to(ROOT)), line=number,
                        category=category, evidence=line.strip(),
                        disposition='version guidance; contextual exception' if exempt else 'review required'))
                    if not exempt:
                        hits.append(category)
        rows.append(dict(slug=page['slug'], language=language,
            title=page['title'][language], purpose=page['group'],
            source_sha256=hashlib.sha256(text.encode()).hexdigest(),
            source_review='review required' if hits else 'static editorial/structure scan passed',
            code='shared snippet/lesson execution recorded separately',
            rendering='HTML links/anchors checked by check-documentation.py',
            visual_scope='browser interaction acceptance cancelled by user; source/R/asset checks only',
            findings='; '.join(hits)))
out = ROOT / 'tools/teaching-audit'
for name, data, fields in [
    (f'content-pages-{series}.csv', rows, list(rows[0])),
    (f'content-findings-{series}.csv', findings,
     ['path', 'line', 'category', 'evidence', 'disposition'])]:
    with (out / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)
templates = ['tools/site_components.py', 'tools/finalize-website.R']
for name in templates:
    text = (ROOT / name).read_text()
    for category, pattern in PATTERNS.items():
        if category == 'classroom wording':
            continue  # Raw templates contain internal lesson/recipe identifiers.
        if re.search(pattern, text, re.I):
            raise SystemExit(f'Review template: {name}: {category}')
bad = [row for row in rows if row['findings']]
if bad:
    raise SystemExit('Editorial review required: ' + ', '.join(row['language']+'/'+row['slug'] for row in bad))
print(f'PASS: {len(rows)} paired source pages; contextual editorial lint (not browser verification).')
