"""Render independently authored R examples and the public parameter directory."""
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


def image_dimensions(name):
    with (ROOT / '.site/assets/handbook' / (name + '.png')).open('rb') as image:
        image.seek(16)
        return struct.unpack('>II', image.read(8))


def parameter_directory(lang):
    flexible = json.loads((WEB / 'flexible-parameters.json').read_text())
    inventory = json.loads((ROOT / '.site/parameter-directory.json').read_text())
    zh = lang == 'zh'
    records = []
    sections = {'sequence':'scalar', 'strand':'width', 'strand-offset':'strand',
                'fixed-label':'scope', 'band':'scope'}
    for row in inventory:
        api, param = row['api'], row['parameter']
        destination = 'reference/' + row['topic'] + '.html#' + row['anchor']
        kind = flexible.get(api, {}).get(param)
        badge = ''
        if kind:
            label = '【支持灵活参数】' if zh else '[Flexible parameter]'
            badge = f'<a class="flexible-parameter" href="parameters.html#{sections[kind]}">{label}</a>'
        default = f'<code>{html.escape(row["default"])}</code>' if row['default'] else '—'
        records.append(f'<tr id="param-{html.escape(api)}-{html.escape(param)}" data-flex="{int(bool(kind))}">'
                       f'<td><span class="argument-term"><code>{html.escape(param)}</code>{badge}</span></td>'
                       f'<td><a href="{destination}"><code>{html.escape(api)}()</code></a></td><td>{default}</td></tr>')
    labels = (('搜索参数或函数', '仅显示支持灵活参数的条目', '参数', '函数', '签名默认值') if zh else
              ('Find a parameter or function', 'Flexible parameters only', 'Parameter', 'Function', 'Signature default'))
    return (f'<div class="parameter-controls" hidden><label>{labels[0]} <input id="directory-query" type="search"></label> '
            f'<label><input id="directory-flex" type="checkbox"> {labels[1]}</label><output aria-live="polite"></output></div>'
            '<table class="parameter-directory"><thead><tr>' + ''.join('<th scope="col">'+s+'</th>' for s in labels[2:]) +
            '</tr></thead><tbody>' + ''.join(records) + '</tbody></table>')


def expand(text, lang, base='../assets'):
    from .source_links import site_document
    from .components import manual_home
    text = site_document(text)
    if '{{manualhome}}' in text:
        text = text.replace('{{manualhome}}', manual_home(lang, base))
    if '{{parameter-directory}}' in text:
        # Preserve literal R defaults such as "\\n" across Markdown conversion.
        text = text.replace('{{parameter-directory}}', '\n```{=html}\n' + parameter_directory(lang) + '\n```\n')
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
            width, height = image_dimensions(name)
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
            width, height = image_dimensions(name)
            size_class = ' wide-output' if entry['width'] >= 10 else ''
            visual = (f'<figure class="documentation-figure handbook-output{size_class}"><a href="{image}" class="figure-link" data-code-url="#example-{name}">'
                      f'<img src="{image}" alt="{html.escape(caption, quote=True)}" width="{width}" height="{height}" loading="lazy"></a>'
                      f'<figcaption>{html.escape(caption)}</figcaption></figure>')
            return '\n' + visual + '\n\n' + fence
        result = (ROOT / '.site/assets/handbook' / (name + '.txt')).read_text().rstrip()
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
