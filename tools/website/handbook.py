"""Render independently authored R examples and their output."""
import html
import json
import re
import struct
from . import ROOT

WEB = ROOT / 'website'
EXAMPLE = re.compile(r'\{\{example:([a-z0-9-]+)\}\}')
PREVIEW = re.compile(r'\{\{preview:([a-z0-9-]+)\}\}')
PREVIEW_GROUP = re.compile(r'(?:\{\{preview:[a-z0-9-]+\}\}\s*)+')


def manifest():
    return json.loads((WEB / 'handbook-examples.json').read_text())


def image_dimensions(folder, name):
    with (folder / (name + '.png')).open('rb') as image:
        image.seek(16)
        return struct.unpack('>II', image.read(8))


def expand(text, lang, base='../assets', *, output=ROOT / '.site/assets/handbook'):
    from .source_links import site_document
    from .components import manual_home
    text = site_document(text)
    if '{{manualhome}}' in text:
        text = text.replace('{{manualhome}}', manual_home(lang, base))
    if '{{case-gallery}}' in text:
        text = text.replace('{{case-gallery}}', case_gallery(lang, base))
    specs = manifest()

    def previews(match):
        figures = []
        for name in PREVIEW.findall(match[0]):
            entry = specs[name]
            caption = html.escape(entry['caption'][lang], quote=True)
            image = f'{base}/handbook/{name}.png'
            target = f'{entry["preview_page"]}.html#example-{name}'
            width, height = image_dimensions(output, name)
            figures.append(f'<figure class="documentation-figure handbook-preview"><a href="{target}">'
                           f'<img src="{image}" alt="{caption}" width="{width}" height="{height}" loading="lazy"></a>'
                           f'<figcaption>{caption}</figcaption></figure>')
        return '<div class="handbook-previews">' + ''.join(figures) + '</div>\n\n'

    text = PREVIEW_GROUP.sub(previews, text)

    def example(match):
        name = match[1]
        entry = specs[name]
        code = (WEB / 'examples' / (name + '.R')).read_text().rstrip()
        kind = entry['kind']
        fence = '\n```{#example-' + name + ' .r .run-example data-example-kind="' + kind + '"}\n' + code + '\n```\n'
        if kind == 'plot':
            image = base + '/handbook/' + name + '.png'
            caption = entry.get('caption', {}).get(lang) or entry.get('label', {}).get(lang) or name
            width, height = image_dimensions(output, name)
            size_class = ' wide-output' if entry['width'] >= 10 else ''
            visual = (f'<figure class="documentation-figure handbook-output{size_class}"><a href="{image}" class="figure-link" data-code-url="#example-{name}">'
                      f'<img src="{image}" alt="{html.escape(caption, quote=True)}" width="{width}" height="{height}" loading="lazy"></a>'
                      f'<figcaption>{html.escape(caption)}</figcaption></figure>')
            return '\n' + visual + '\n\n' + fence
        result = (output / (name + '.txt')).read_text().rstrip()
        caption = html.escape(entry['caption'][lang])
        return fence + f'\n<figure class="data-output"><figcaption>{caption}</figcaption><pre>{html.escape(result)}</pre></figure>\n'

    return EXAMPLE.sub(example, text)


def case_figure(page):
    examples = EXAMPLE.findall((WEB / 'en' / (page + '.md')).read_text())
    if len(examples) != 1:
        raise ValueError(f'{page}: expected one complete figure')
    return examples[0]


def case_gallery(lang, base):
    pages = {p['slug']: p for p in json.loads((WEB / 'pages.json').read_text())}
    specs = manifest()
    entries = []
    for guide in ('chord', 'collinear', 'radial-synteny', 'sankey', 'circular', 'dotplot', 'fusion'):
        case = 'case-' + ('radial' if guide == 'radial-synteny' else guide)
        key = case_figure(case)
        title = html.escape(pages[guide]['title'][lang])
        caption = html.escape(specs[key]['caption'][lang], quote=True)
        label = html.escape(pages[case]['title'][lang])
        action = '试画简图' if lang == 'zh' else 'Draw a small example'
        entries.append(f'<section class="gallery-card"><h2 id="gallery-{guide}">{title}</h2>'
                       f'<a href="{case}.html"><img loading="lazy" src="{base}/handbook/{key}.png" alt="{caption}"></a>'
                       f'<p><a href="{case}.html">{label}</a></p><a href="{guide}.html">{action}</a></section>')
    return '<div class="handbook-gallery">' + ''.join(entries) + '</div>'
