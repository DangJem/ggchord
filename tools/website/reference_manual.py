"""Render the paired v0.50.0 API manual into the documentation site."""
import argparse
import hashlib
import html
import json
import os
import re
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory

from . import ROOT

SOURCE = ROOT / 'website/reference-manual'
ASSETS = ROOT / 'website/assets'
LOCALES = ('en', 'zh')
NATIVE = {'facet_null', 'facet_wrap', 'facet_grid'}
LABELS = {
    'en': dict(reference='API reference',
        usage='Usage', arguments='Arguments', flexible='[Flexible parameter]',
        purpose='Description', value='Return value', details='Behavior', errors='Errors',
        roles='Data roles and aesthetics', examples='Examples', related='See also',
        source='Native ggplot2', download='Download example', filter='Filter by function name or purpose',
        console='Console output', figure='Plot from the example for',
        input='Input tables', setup='Shared example data',
        upstream='Native ggplot2 tools',
        archived='Historical API address', current='Open the current API reference'),
    'zh': dict(reference='API 参考',
        usage='调用形式', arguments='参数', flexible='【支持灵活参数】',
        purpose='说明', value='返回值', details='行为规则', errors='错误与限制',
        roles='数据角色与视觉映射', examples='示例', related='相关接口',
        source='原生 ggplot2', download='下载示例', filter='按函数名或用途筛选',
        console='控制台输出', figure='示例绘图',
        input='输入表', setup='共用示例数据',
        upstream='原生 ggplot2 工具',
        archived='历史 API 地址', current='打开当前 API 参考'),
}
FLEXIBLE_SECTIONS = {'sequence': 'scalar', 'strand': 'width', 'strand-offset': 'strand',
                     'fixed-label': 'scope', 'band': 'scope'}
UPSTREAM = (
    ('aes(), after_stat()', 'aes_eval.html'), ('vars()', 'vars.html'),
    ('scale_*()', 'index.html#scales'), ('guides()', 'guides.html'),
    ('position_*()', 'layer_positions.html'), ('theme()', 'theme.html'),
    ('labs()', 'labs.html'), ('ggplot()', 'ggplot.html'), ('ggsave()', 'ggsave.html'),
)
INPUT_ROLES = (
    ('Reference', '参考', ('accver', 'length', 'sequence_unit')),
    ('Interval feature', '区间特征', ('accver', 'start', 'end', 'strand')),
    ('Site', '点位', ('accver', 'pos', 'cut_top', 'cut_bottom')),
    ('Point link', '点连接', ('qaccver', 'qpos', 'saccver', 'spos')),
    ('Interval relationship', '区间关系', ('qaccver', 'qstart', 'qend', 'saccver', 'sstart', 'send')),
    ('Numeric track', '数值轨道', ('accver', 'pos', 'value')),
    ('Character cell', '字符单元', ('accver', 'alignment_column', 'symbol', 'source_position', 'is_gap')),
)


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding='utf-8')


def write_json(path, data):
    write(path, json.dumps(data, ensure_ascii=False, indent=2) + '\n')


def link(href, name):
    return f'<a href="{html.escape(href, quote=True)}">{html.escape(name)}</a>'


def code(s):
    return f'<code>{html.escape(s)}</code>'


def load():
    structure = read(SOURCE / 'structure.json')
    locales = {lang: read(SOURCE / f'{lang}.json') for lang in LOCALES}
    digest = hashlib.sha256((SOURCE / 'en.json').read_bytes()).hexdigest()
    if locales['zh']['source_sha256'] != digest:
        raise ValueError('The Han reference requires review against the current English source')
    topics = structure['topics']
    exports = set(re.findall(r'^export\(([^)]+)\)', (ROOT / 'NAMESPACE').read_text(), re.M))
    if set(topics) != exports | NATIVE:
        raise ValueError(f'Reference inventory mismatch: {sorted(set(topics) ^ (exports | NATIVE))}')
    grouped = [topic for group in structure['groups'].values() for topic in group['topics']]
    if len(grouped) != len(set(grouped)) or set(grouped) != set(topics):
        raise ValueError('Each reference topic must appear in exactly one category')
    for lang, locale in locales.items():
        if not locale['description'].strip():
            raise ValueError(f'{lang}: missing reference introduction')
        if set(locale['groups']) != set(structure['groups']) or set(locale['topics']) != set(topics):
            raise ValueError(f'{lang}: categories or topics missing')
        for key, part in topics.items():
            body = locale['topics'][key]
            if not re.fullmatch(r'[a-z][a-z0-9_]+', key) or not set(part['related']) <= set(topics):
                raise ValueError(f'{key}: invalid topic identity or related link')
            if not part['usage'].startswith(('ggplot2::' if part['native'] else '') + key + '('):
                raise ValueError(f'{key}: invocation does not match topic')
            if not part['examples'].strip() or not all(body[field].strip() for field in
                    ('title', 'intro', 'value', 'details')):
                raise ValueError(f'{lang}/{key}: example or core prose missing')
            if len(body['arguments']) != len(part['arguments']) or not all(
                    row['description'].strip() and row['format'].strip() for row in body['arguments']):
                raise ValueError(f'{lang}/{key}: argument descriptions missing')
            if len(body['roles']) != len(part['roles']) or not all(body['roles']):
                raise ValueError(f'{lang}/{key}: data roles missing')
            if not all(re.fullmatch(r'(?:[A-Za-z_.][A-Za-z_.0-9]*|\.\.\.)', name)
                       for row in part['arguments'] for name in row):
                raise ValueError(f'{key}: invalid parameter name')
            if set(part['defaults']) != {name for row in part['arguments'] for name in row}:
                raise ValueError(f'{key}: parameter defaults differ from the argument list')
    return structure, locales


@lru_cache(maxsize=4096)
def markdown(text, pandoc):
    # Translate handbook links to the corresponding page in the same locale.
    text = re.sub(r'\]\(#api-([a-z0-9_]+)\)', r'](\1.html)', text)
    return subprocess.check_output(
        [pandoc, '--from=markdown+pipe_tables', '--to=html5', '--wrap=none'],
        input=text, text=True)


def section(id, title, body):
    return f'<section aria-labelledby="{id}"><h2 id="{id}">{html.escape(title)}</h2>{body}</section>'


def snippet(markup, kind, caption):
    # Keep the prepared token spans through pkgdown's article rendering pass.
    return (f'<details class="code-disclosure" data-code-kind="{kind}" open>'
            f'<summary>{html.escape(caption)}</summary><pre class="downlit"><code>{markup}'
            '</code></pre></details>')


def topic_link(topic):
    return f'<a href="{topic}.html">{code(topic + "()")}</a>'


def invocation(topic):
    name = topic['usage'].split('(', 1)[0]
    lines = [f'  {key} = {value}' if value else f'  {key}'
             for key, value in topic['defaults'].items()]
    return name + '(\n' + ',\n'.join(lines) + '\n)' if lines else name + '()'


def prepare_examples(structure):
    setup = structure['setup'].strip().split('\n\n')
    # The fixture uses top-level assignments. Find the named objects an example
    # needs and include dependencies in the original order.
    names = re.findall(r'^([A-Za-z0-9_]+) <- ', structure['setup'], re.M)
    blocks = {}
    for match in re.finditer(r'^([A-Za-z0-9_]+) <- (.*?)(?=^[A-Za-z0-9_]+ <- |\Z)',
                             structure['setup'], re.M | re.S):
        blocks[match[1]] = match[0].strip()
    def words(s):
        return set(re.findall(r'\b[A-Za-z_][A-Za-z0-9_]*\b', s))
    scripts = {}
    for topic, page in structure['topics'].items():
        setup_source = setup[0]
        example = page['examples']
        requested = words(example)
        for name in reversed(names):
            if name in requested:
                requested |= words(blocks[name])
        selected = [blocks[name] for name in names if name in requested]
        scripts[topic] = '\n\n'.join([setup_source, *selected, example]) + '\n'
    return scripts


def render_examples(assets, structure):
    scripts = prepare_examples(structure)
    setup = structure['setup'].strip() + '\n'
    with TemporaryDirectory(prefix='ggchord-reference-examples-') as temporary:
        staging = Path(temporary)
        inputs = staging / 'inputs.json'
        generated = staging / 'output'
        write_json(inputs, dict(setup=setup, topics={
            name: dict(signature=invocation(topic), example=scripts[name])
            for name, topic in structure['topics'].items()}))
        subprocess.run(['Rscript', '--vanilla', 'tools/website/render-reference-examples.R',
                        str(inputs), str(generated)], cwd=ROOT, check=True)
        results = read(generated / 'results.json')
        folder = assets / 'reference-examples'
        for name, source in scripts.items():
            write(folder / f'{name}.R', source)
        write(folder / 'setup.R', setup)
        figures = {figure for entry in results['topics'].values() for figure in entry['figures']}
        for obsolete in folder.glob('*-[0-9][0-9].png'):
            if obsolete.name not in figures:
                obsolete.unlink()
        for figure in figures:
            shutil.copy2(generated / figure, folder / figure)
    return results


def example_output(rendered, name, titles, width, height):
    blocks = []
    if rendered['console']:
        blocks.append(f'<figure class="data-output"><figcaption>{titles["console"]}</figcaption>'
                      f'<pre class="downlit"><code>{html.escape(rendered["console"])}</code></pre></figure>')
    for index, figure in enumerate(rendered['figures'], 1):
        caption = f'{titles["figure"]} {name}()'
        if len(rendered['figures']) > 1:
            caption += f' ({index})'
        blocks.append(f'<figure class="documentation-figure"><img loading="lazy" '
                      f'src="../../assets/reference-examples/{figure}" alt="{html.escape(caption, quote=True)}" '
                      f'width="{width}" height="{height}"><figcaption>{html.escape(caption)}</figcaption></figure>')
    return ''.join(blocks)


def render_articles(root, pages):
    # Use the same pkgdown shell and navigation finalizer as the other chapters.
    with TemporaryDirectory(prefix='ggchord-reference-') as temporary:
        staging = Path(temporary)
        payload = staging / 'pages.json'
        write_json(payload, pages)
        # pkgdown protects files made by other generators. Render into a fresh
        # directory so old standalone pages cannot survive a successful build.
        options = [f'--site={staging}', f'--payload={payload}']
        for tool, flags in (('render-website.R', ['--part=articles']),
                            ('finalize-website.R', [f'--navigation-site={root}'])):
            result = subprocess.run(['Rscript', '--vanilla', f'tools/{tool}', *flags, *options],
                                    cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode:
                raise RuntimeError(result.stdout.rstrip())
        for lang, entries in pages.items():
            rendered = staging / lang
            target = root / lang
            for dependency in ('deps', 'lightswitch.js', 'pkgdown.js'):
                asset = rendered / dependency
                if asset.is_dir():
                    for source in (p for p in asset.rglob('*') if p.is_file()):
                        path = target / source.relative_to(rendered)
                        if not path.exists():
                            path.parent.mkdir(parents=True, exist_ok=True)
                            shutil.copy2(source, path)
                elif not (target / dependency).exists():
                    shutil.copy2(asset, target / dependency)
            for page in entries:
                relative = page['slug'] + '.html'
                path = target / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(rendered / relative, path)


def prepare(root, pandoc):
    """Prepare articles, search rows and evaluated examples for either build mode."""
    structure, locales = load()
    assets = root / 'assets'
    output = render_examples(assets, structure)
    search_by_language = {}
    flexible = read(ROOT / 'website/flexible-parameters.json')
    pages = {lang: [] for lang in LOCALES}
    for lang, locale in locales.items():
        topics, groups = locale['topics'], locale['groups']
        titles = LABELS[lang]
        search = []
        intro = locale['description']
        contents = ''.join(
            '<section class="reference-group">' +
            ''.join(f'<span id="{alias}"></span>' for alias in group['previous_anchors']) +
            f'<h2 id="family-{key}">{html.escape(groups[key]["title"])}</h2>'
            f'<p>{html.escape(groups[key]["intro"])}</p><dl>' +
            ''.join(f'<dt>{topic_link(name)}</dt><dd>{html.escape(topics[name]["title"])}</dd>'
                    for name in group['topics']) + '</dl></section>'
            for key, group in structure['groups'].items())
        field = ('<div class="reference-filter" hidden><label for="reference-filter">'
                 f'{titles["filter"]}</label><input id="reference-filter" type="search" autocomplete="off">'
                 '<output id="reference-filter-count" aria-live="polite"></output></div>')
        upstream = '<ul>' + ''.join(
            f'<li>{link("https://ggplot2.tidyverse.org/reference/" + href, name)}</li>'
            for name, href in UPSTREAM) + '</ul>'
        roles = '<dl>' + ''.join(
            f'<dt>{html.escape(names[lang == "zh"])}</dt>'
            f'<dd>{", ".join(code(key) for key in names[2])}</dd>' for names in INPUT_ROLES) + '</dl>'
        units = ('Closed source intervals are one-based. Documented cut boundaries use interbase coordinates from zero through the reference length.'
                 if lang == 'en' else '源区间使用一基闭区间。已知切口采用残基间坐标，范围为 0 到完整参考长度。')
        setup = (snippet(output['setup'], 'r', titles['setup']) +
                 f'<a href="../../assets/reference-examples/setup.R" download>{titles["download"]}</a>')
        index = (f'<p class="lead">{html.escape(intro)}</p>' +
                 field + contents + section('input-roles', titles['input'], roles + f'<p>{html.escape(units)}</p>') +
                 section('example-setup', titles['setup'], setup) + section('upstream', titles['upstream'], upstream))
        pages[lang].append(dict(slug='reference/index', title=titles['reference'], body=index))
        search.append(dict(url='reference/index.html',title=titles['reference'],kind='tool',
                           summary=intro, text=' '.join(groups[k]['intro'] for k in groups),keywords='API reference'))
        for family, group in structure['groups'].items():
            for name in group['topics']:
                shared = structure['topics'][name]
                text = topics[name]
                native = (f'<p>{link("https://ggplot2.tidyverse.org/reference/" + name + ".html", titles["source"])}</p>'
                          if shared['native'] else '')
                bodies = [native]
                def add(key, fragment):
                    if fragment:
                        bodies.append(section('ref-' + key, titles[key], fragment))
                add('purpose', '<span id="ref-description"></span>' + markdown(text['intro'], pandoc))
                rendered = output['topics'][name]
                add('usage', snippet(rendered['signature'], 'signature', titles['usage']))
                rows = []
                seen = set()
                for keys, param, original in zip(shared['arguments'], text['arguments'],
                                                 locales['en']['topics'][name]['arguments']):
                    description = markdown(param['description'], pandoc)
                    summary = ' '.join(html.unescape(re.sub(r'<[^>]*>', ' ', description)).split())
                    names = []
                    for key in keys:
                        anchor = 'para-ellipsis' if key == '...' else 'parameter-' + key
                        if key in seen:
                            names.append(code(key))
                            continue
                        seen.add(key)
                        alias = f'<span id="arg-{key}"></span>' if key != '...' else ''
                        badge = ''
                        kind = flexible.get(name, {}).get(key)
                        if kind or re.search(r'\bflexible\b', original['format']):
                            fragment = FLEXIBLE_SECTIONS.get(kind, 'scope')
                            badge = (f'<a class="flexible-parameter" href="../parameters.html#{fragment}">'
                                     f'{titles["flexible"]}</a>')
                        names.append(f'<span class="argument-term">{alias}'
                                     f'<code id="{anchor}" class="argument-anchor">{html.escape(key)}</code>'
                                     f'{badge}</span>')
                        search.append(dict(url=f'reference/{name}.html#{anchor}', kind='parameter',
                                           title=f'{name}() · {key}', keywords=f'{name} {key}',
                                           text=f'{param["format"]} {summary}', summary=summary))
                    rows.append(f'<dt>{", ".join(names)}</dt><dd>'
                                f'<p class="argument-format">{html.escape(param["format"])}</p>'
                                f'{description}</dd>')
                if rows:
                    add('arguments', f'<dl class="reference-arguments">{"".join(rows)}</dl>')
                if shared['roles']:
                    roles = ''.join(f'<dt>{html.escape(tag.replace(chr(96), ""))}</dt>'
                                    f'<dd>{markdown(desc,pandoc)}</dd>' for tag,desc in zip(shared['roles'], text['roles']))
                    add('roles', f'<dl>{roles}</dl>')
                add('value', markdown(text['value'], pandoc))
                add('details', markdown(text['details'], pandoc))
                if text['errors']:
                    add('errors', markdown(text['errors'], pandoc))
                example_note = markdown(text['example_setup'],pandoc) if text['example_setup'] else ''
                downloads = (f'<a href="../../assets/reference-examples/{name}.R" download>'
                             f'{titles["download"]}</a>')
                add('examples', example_note + snippet(rendered['example'], 'r', titles['examples']) +
                    downloads + example_output(rendered, name, titles, output['width'], output['height']))
                add('related', '<ul>' + ''.join(f'<li>{topic_link(target)}</li>' for target in shared['related']) + '</ul>')
                pages[lang].append(dict(slug=f'reference/{name}', title=name + '()',
                                        heading_id=f'api-{name}', body=''.join(bodies)))
                search.append(dict(url=f'reference/{name}.html',title=f'{name}() — {text["title"]}',
                                   kind='function',summary=text['intro'], text=text['details'],
                                   keywords=' '.join([name, family, *seen])))
        search_by_language[lang] = search
    version = re.search(r'^Version: (.+)$', (ROOT / 'DESCRIPTION').read_text(), re.M)[1]
    write_json(root / 'api-reference-manifest.json', {
        'version': version, 'source': structure['reference'], 'functions': len(structure['topics']) - len(NATIVE),
        'native_facets': len(NATIVE), 'categories': len(structure['groups']), 'locales': list(LOCALES),
        'topics': list(structure['topics']), 'deployment': 'local output',
    })
    write_json(root / 'reference-output.json', {
        'status': 'passed', 'topics': {name: {'console': entry['console'], 'figures': entry['figures']}
                                     for name, entry in output['topics'].items()},
    })
    return pages, search_by_language


def write_reference_redirect(out, path, destination, languages):
    labels = ({'en': 'Version selection', 'zh': '版本选择'} if destination == 'migration.html#api'
              else {'en': 'Current documentation', 'zh': '现行文档'})
    links = []
    for lang in languages:
        relative = os.path.relpath(out/lang/destination, path.parent)
        tag = 'zh-Hans' if lang == 'zh' else 'en'
        links.append('<a lang="' + tag + '" hreflang="' + tag + '" href="' +
                     html.escape(relative, quote=True) + '">' + labels[lang] + '</a>')
    target = os.path.relpath(out/languages[0]/destination, path.parent)
    script = '<script>location.replace(' + json.dumps(target) + ')</script>'
    tag = 'zh-Hans' if languages == ('zh',) else 'en'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('<!doctype html><html lang="' + tag + '"><head><meta charset="utf-8">'
                    '<title>' + labels[languages[0]] + ' · ggchord</title>' + script +
                    '</head><body><main>' + ' · '.join(links) + '</main></body></html>')


def build(root):
    root = root.resolve()
    assets = root / 'assets'
    assets.mkdir(parents=True, exist_ok=True)
    for name in ('site.css', 'site.js', 'favicon.png', 'logo.webp'):
        shutil.copy2(ASSETS / name, assets / name)
    pandoc = shutil.which('pandoc')
    if not pandoc:
        raise RuntimeError('Pandoc must be on PATH to build the API reference')
    pages, searches = prepare(root, pandoc)
    for lang, search in searches.items():
        home = root / lang / 'index.html'
        if not home.exists():
            write(home, redirect('reference/index.html'))
        search_file = root / lang / 'search.json'
        previous = read(search_file) if search_file.exists() else []
        previous = [row for row in previous if
                    row['url'].split('#')[0] != 'parameter-atlas.html' and
                    not row['url'].startswith('reference/')]
        write_json(search_file, previous + search)
    render_articles(root, pages)
    for folder in (root, *(root / lang for lang in LOCALES)):
        (folder / 'parameter-atlas.html').unlink(missing_ok=True)
    subprocess.run(['Rscript', '--vanilla', 'tools/finalize-website.R',
                    f'--site={root}', '--refresh-existing'], cwd=ROOT, check=True)
    for page in pages['en']:
        name = page['slug'] + '.html'
        write(root / name, redirect('../en/' + name))
    migrations = read(ROOT/'website/page-migrations.json')
    for path, mapping in migrations.items():
        if not path.startswith('reference/'):
            continue
        # A fresh standalone reference may not yet contain the handbook pages.
        for lang in LOCALES:
            destination = mapping['']
            if not (root / lang / destination.split('#')[0]).is_file():
                destination = 'reference/index.html'
            write_reference_redirect(root, root/lang/path, destination, (lang,))
        write(root/path, redirect('../en/'+path))
    if not (root / 'index.html').exists():
        write(root / 'index.html', redirect('en/reference/index.html'))
    print(f'Built {len(pages["en"]) - 1} API entries in each language at {root}')


def redirect(destination):
    return (f'<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" '
            f'content="0; url={destination}"><title>ggchord API reference</title>'
            f'<a href="{destination}">API reference</a>'
            '<script>if (location.hash) location.replace(document.querySelector("a").href + location.hash)</script>\n')


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument('--output', type=Path, default=ROOT / '.site')
    args.add_argument('--check', action='store_true')
    options = args.parse_args()
    if options.check:
        structure, _ = load()
        print(f'PASS: {len(structure["topics"])} paired API entries')
    else:
        build(options.output)


if __name__ == '__main__':
    main()
