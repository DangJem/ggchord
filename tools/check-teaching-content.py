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
}
rows, findings = [], []
for page in PAGES:
    for language in ('en', 'zh'):
        path = ROOT / 'website' / language / (page['slug'] + '.md')
        text = path.read_text()
        fenced = False
        hits = []
        for number, line in enumerate(text.splitlines(), 1):
            if re.match(r'^\s*(```|~~~)', line):
                fenced = not fenced
                continue
            if fenced:
                continue
            for category, pattern in PATTERNS.items():
                if re.search(pattern, line, re.I):
                    # Installation/migration explain versions intentionally.
                    exempt = page['slug'] in ('installation', 'migration')
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
            visual_scope='representative browser pages; not every page individually',
            findings='; '.join(hits)))
out = ROOT / 'tools/teaching-audit'
for name, data, fields in [
    ('content-pages-v0.17.0.csv', rows, list(rows[0])),
    ('content-findings-v0.17.0.csv', findings,
     ['path', 'line', 'category', 'evidence', 'disposition'])]:
    with (out / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)
templates = ['tools/site_components.py', 'tools/finalize-website.R']
for name in templates:
    text = (ROOT / name).read_text()
    for category, pattern in PATTERNS.items():
        if re.search(pattern, text, re.I):
            raise SystemExit(f'Review template: {name}: {category}')
bad = [row for row in rows if row['findings']]
if bad:
    raise SystemExit('Editorial review required: ' + ', '.join(row['language']+'/'+row['slug'] for row in bad))
print(f'PASS: {len(rows)} paired source pages; contextual editorial lint (not browser verification).')
