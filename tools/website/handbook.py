"""Render independently authored R examples and the public parameter directory."""
import html
import json
import re
from . import ROOT

WEB = ROOT / 'website'
EXAMPLE = re.compile(r'\{\{example:([a-z0-9-]+)\}\}')


def manifest():
    return json.loads((WEB / 'handbook-examples.json').read_text())


def parameter_directory(lang):
    flexible = json.loads((WEB / 'flexible-parameters.json').read_text())
    inventory = json.loads((ROOT / '.site/parameter-directory.json').read_text())
    zh = lang == 'zh'
    records = []
    for row in inventory:
        api, param = row['api'], row['parameter']
        destination = 'reference/' + row['topic'] + '.html#' + row['anchor']
        supported = param in flexible.get(api, {})
        flag = '<a class="flexible-label" href="parameters.html">' + ('灵活参数' if zh else 'Flexible parameters') + '</a>' if supported else ''
        records.append(f'<tr id="param-{html.escape(api)}-{html.escape(param)}" data-flex="{int(supported)}"><td><code>{html.escape(param)}</code></td><td><a href="{destination}"><code>{html.escape(api)}()</code></a></td><td>{flag}</td></tr>')
    labels = (('搜索参数或函数', '仅显示支持灵活参数的条目', '参数', '函数', '用法') if zh else
              ('Find a parameter or function', 'Flexible parameters only', 'Parameter', 'Function', 'Usage'))
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
        text = text.replace('{{parameter-directory}}', parameter_directory(lang))
    if '{{case-gallery}}' in text:
        text = text.replace('{{case-gallery}}', case_gallery(lang, base))
    specs = manifest()

    def example(match):
        name = match[1]
        entry = specs[name]
        code = (WEB / 'examples' / (name + '.R')).read_text().rstrip()
        kind = entry['kind']
        fence = '\n```{#example-' + name + ' .r .run-example data-example-kind="' + kind + '"}\n' + code + '\n```\n'
        if kind == 'plot':
            image = base + '/handbook/' + name + '.png'
            caption = entry.get('caption', {}).get(lang) or entry.get('label', {}).get(lang) or name
            width, height = (round(entry.get(key, default) * entry.get('dpi', 144)) for key, default in [('width',7),('height',5)])
            visual = (f'<figure class="documentation-figure handbook-output"><a href="{image}" class="figure-link" data-code-url="#example-{name}">'
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
