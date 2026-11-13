"""Regression check for the final native gallery, not package release tests."""
from html.parser import HTMLParser
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]


class Gallery(HTMLParser):
    def __init__(self):
        super().__init__()
        self.key = None
        self.links = {}

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'figure':
            self.key = attrs.get('data-figure')
        if tag == 'a' and self.key:
            if self.key in self.links:
                raise AssertionError(f'Duplicate figure link: {self.key}')
            self.links[self.key] = attrs

    def handle_endtag(self, tag):
        if tag == 'figure':
            self.key = None


expected = set(json.loads((ROOT/'website/figure-lessons.json').read_text()))
results = {}
for lang in ('en', 'zh'):
    path = ROOT/'.site'/lang/'gallery.html'
    page = Gallery()
    page.feed(path.read_text())
    assert set(page.links) == expected
    for key, attrs in page.links.items():
        assert attrs['href'] == f'../assets/{key}.png', (lang, key)
        assert attrs['data-code-url'] == '../assets/documentation-examples.R', (lang, key)
    results[lang] = dict(figures=len(page.links),
                        html_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
report = dict(status='passed', scope='Final native gallery figure URLs in two languages',
              pages=results, not_checked=['browser interaction and visual quality'])
(ROOT/'tools/long-term-audit/v0.27.0-gallery-source-link-results.json').write_text(
    json.dumps(report, indent=2)+'\n')
print('PASS: both native galleries retain shared R source links for all figures')
