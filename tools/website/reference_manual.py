"""Render the paired v0.50.0 API manual into the documentation site."""
import argparse
import hashlib
import html
import json
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
        parameter_index='Parameter index', usage='Usage', arguments='Arguments', format='Type / format', default='Default',
        purpose='Description', value='Return value', details='Behavior', errors='Errors',
        roles='Data roles and aesthetics', examples='Examples', related='See also',
        source='Native ggplot2', download='Download example', filter='Filter by function name or purpose',
        search_parameter='Filter by parameter or function', flex='Flexible parameters only',
        overview='How the API fits together', input='Input tables', setup='Shared example data',
        upstream='Native ggplot2 tools',
        archived='Historical API address', current='Open the current API reference'),
    'zh': dict(reference='API 参考',
        parameter_index='参数索引', usage='调用形式', arguments='参数', format='类型与格式', default='默认值',
        purpose='说明', value='返回值', details='行为规则', errors='错误与限制',
        roles='数据角色与视觉映射', examples='示例', related='相关接口',
        source='原生 ggplot2', download='下载示例', filter='按函数名或用途筛选',
        search_parameter='按参数或函数筛选', flex='只显示灵活参数',
        overview='接口如何配合', input='输入表', setup='共用示例数据',
        upstream='原生 ggplot2 工具',
        archived='历史 API 地址', current='打开当前 API 参考'),
}
GRAMMAR = {
    'en': [
        ('References and layers', 'Register a complete sequence dictionary with ggchord(data = seqs) or with an explicit geom_seq(data = seqs) layer on an empty plot. Feature, relationship and track observations are supplied to their own layers.'),
        ('Layout, coordinates and facets', 'Layout declares order, orientation and display windows. Coordinates project the reference paths. Native facets assign observations to panels. Circular uses one effective reference per panel; free native axes require coord_dotplot(ratio = NULL).'),
        ('Appearance and statistics', 'Aesthetic mappings and visual scales control appearance. Track height has its own numerical scale. Statistics use source coordinates and keep their contributing records, while display windows select what to show.'),
        ('Independent views and output', 'Fusion combines named panels. Layout export reports source, mapped, statistic and geometry stages; display queries record their size and measurement conditions.'),
    ],
    'zh': [
        ('参考与图层', '用 ggchord(data = seqs) 登记完整参考，或在空图上用 geom_seq(data = seqs) 显式登记。特征、关系和轨道观察表分别交给各自图层。'),
        ('布局、坐标与分面', 'Layout 声明顺序、方向和显示窗口，Coord 投影参考路径，原生 Facet 按数据分组。Circular 每个面板需要一条有效参考。原生自由轴使用 coord_dotplot(ratio = NULL)。'),
        ('视觉映射与统计', '视觉属性由映射和标度控制，轨道高度有独立数值标度。统计按源坐标计算并保留来源成员，显示窗口决定呈现哪些区间。'),
        ('独立视图与输出', 'Fusion 组合具名面板。布局导出区分来源、映射、统计和几何四阶段，显示位置查询记录尺寸与测量条件。'),
    ],
}
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
                raise ValueError(f'{key}: parameter defaults differ from the argument table')
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


def snippet(content, kind, caption):
    # pkgdown's article highlighter expands literal tabs, including tabs in
    # file-format fixture strings. Keep the code text intact for the copy action.
    return (f'<details class="code-disclosure" data-code-kind="{kind}" open>'
            f'<summary>{html.escape(caption)}</summary><pre class="downlit"><code>{html.escape(content)}'
            '</code></pre></details>')


def topic_link(topic):
    return f'<a href="{topic}.html">{code(topic + "()")}</a>'


def syntax_format(value):
    return ' '.join(value.split())


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


def build(root):
    structure, locales = load()
    pandoc = shutil.which('pandoc')
    if not pandoc:
        raise RuntimeError('Pandoc must be on PATH to build the API reference')
    root = root.resolve()
    version = re.search(r'^Version: (.+)$', (ROOT / 'DESCRIPTION').read_text(), re.M)[1]
    assets = root / 'assets'
    assets.mkdir(parents=True, exist_ok=True)
    shared_assets = ['site.css', 'site.js', 'favicon.png', 'logo.webp']
    for name in shared_assets:
        shutil.copy2(ASSETS / name, assets / name)
    scripts = prepare_examples(structure)
    for name, text in scripts.items():
        write(assets / 'reference-examples' / f'{name}.R', text)
    write(assets / 'reference-examples' / 'setup.R', structure['setup'].strip() + '\n')
    old_names = {p.stem for p in (ROOT / 'website/reference-zh').glob('*.json')}
    existing_help = {p.stem for p in (ROOT / 'man').glob('*.Rd')}
    retired = old_names - existing_help - set(structure['topics'])
    flexible = read(ROOT / 'website/flexible-parameters.json')
    pages = {lang: [] for lang in LOCALES}
    for lang, locale in locales.items():
        topics, groups = locale['topics'], locale['groups']
        titles = LABELS[lang]
        home = root / lang / 'index.html'
        if not home.exists():
            write(home, redirect('reference/index.html'))
        search = []
        overview = ''.join(f'<dt>{html.escape(k)}</dt><dd>{html.escape(v)}</dd>' for k,v in GRAMMAR[lang])
        total = len(topics)
        intro = (f'ggchord {version} · {total} API entries. The faceting entries use ggplot2 directly.'
                 if lang == 'en' else f'ggchord {version} · {total} 个接口条目。分面入口直接使用 ggplot2。')
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
        setup = (snippet(structure['setup'], 'r', titles['setup']) +
                 f'<a href="../../assets/reference-examples/setup.R" download>{titles["download"]}</a>')
        index = (f'<p class="lead">{html.escape(intro)}</p>' +
                 section('grammar', titles['overview'], f'<dl>{overview}</dl>') +
                 field + contents + section('input-roles', titles['input'], roles + f'<p>{html.escape(units)}</p>') +
                 section('example-setup', titles['setup'], setup) + section('upstream', titles['upstream'], upstream))
        pages[lang].append(dict(slug='reference/index', title=titles['reference'], body=index))
        search.append(dict(url='reference/index.html',title=titles['reference'],kind='api',
                           summary=intro, text=' '.join(groups[k]['intro'] for k in groups),keywords='API reference'))
        parameter_rows = []
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
                add('usage', snippet(shared['usage'], 'signature', titles['usage']))
                rows = []
                seen = set()
                for keys, param in zip(shared['arguments'], text['arguments']):
                    names = []
                    for key in keys:
                        anchor = 'para-ellipsis' if key == '...' else 'parameter-' + key
                        if key in seen:
                            names.append(code(key))
                            continue
                        seen.add(key)
                        names.append(f'<span class="parameter-name" id="{anchor}">{code(key)}</span>')
                        if key != '...':
                            names[-1] += f'<span id="arg-{key}"></span>'
                        is_flex = key in flexible.get(name, {}) or bool(re.search(
                            'flexible specification|灵活配置',param['format']))
                        parameter_rows.append((key, name, anchor, is_flex,
                                               syntax_format(param['format']), shared['defaults'][key]))
                    rows.append(f'<tr><th scope="row">{", ".join(names)}</th>'
                                f'<td>{html.escape(syntax_format(param["format"])) or "—"}</td>'
                                f'<td>{markdown(param["description"],pandoc)}</td></tr>')
                table = ('<table><thead><tr><th>'
                         f'{titles["arguments"]}</th><th>{titles["format"]}</th>'
                         f'<th>{titles["purpose"]}</th></tr></thead><tbody>{"".join(rows)}</tbody></table>')
                add('arguments', table)
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
                add('examples', example_note + snippet(scripts[name], 'r', titles['examples']) + downloads)
                add('related', '<ul>' + ''.join(f'<li>{topic_link(target)}</li>' for target in shared['related']) + '</ul>')
                pages[lang].append(dict(slug=f'reference/{name}', title=name + '()',
                                        heading_id=f'api-{name}', body=''.join(bodies)))
                search.append(dict(url=f'reference/{name}.html',title=f'{name}() — {text["title"]}',
                                   kind='api',summary=text['intro'], text=text['details'],
                                   keywords=' '.join([name, family, *seen])))
        ptitle = titles['parameter_index']
        controls = ('<div class="parameter-controls" hidden><label for="directory-query">'
                    f'{titles["search_parameter"]}</label><input id="directory-query" type="search">'
                    f'<label><input type="checkbox" id="directory-flex"> {titles["flex"]}</label>'
                    '<output aria-live="polite"></output></div>')
        table = ''.join(f'<tr data-flex="{int(flex)}"><th scope="row">{code(key)}</th>'
                        f'<td><a href="reference/{topic}.html#{anchor}">{code(topic + "()")}</a></td>'
                        f'<td>{html.escape(fmt)}</td><td>{code(default) if default else "—"}</td></tr>'
                        for key,topic,anchor,flex,fmt,default in sorted(parameter_rows))
        body = (f'{controls}<table class="parameter-directory">'
                f'<thead><tr><th>{titles["arguments"]}</th><th>API</th><th>{titles["format"]}</th>'
                f'<th>{titles["default"]}</th></tr></thead>'
                f'<tbody>{table}</tbody></table>')
        pages[lang].append(dict(slug='parameter-atlas', title=ptitle, body=body))
        for name in sorted(retired):
            title = titles['archived']
            body = f'<p>{title}</p><p>{link("index.html", titles["current"])}</p>'
            pages[lang].append(dict(slug=f'reference/{name}', title=name + '()', body=body))
        search_file = root / lang / 'search.json'
        previous = read(search_file) if search_file.exists() else []
        previous = [row for row in previous if not row['url'].startswith('reference/') or
                    Path(row['url'].split('#')[0]).stem in existing_help - set(topics)]
        write_json(search_file, previous + search)
    render_articles(root, pages)
    for name in ('index', *structure['topics'], *sorted(retired)):
        destination = f'../en/reference/{name}.html'
        write(root / 'reference' / f'{name}.html', redirect(destination))
    if not (root / 'index.html').exists():
        write(root / 'index.html', redirect('en/reference/index.html'))
    write_json(root / 'api-reference-manifest.json', {
        'version': version, 'source': structure['reference'], 'functions': len(structure['topics']) - len(NATIVE),
        'native_facets': len(NATIVE), 'categories': len(structure['groups']), 'locales': list(LOCALES),
        'topics': list(structure['topics']), 'deployment': 'local output',
    })
    print(f'Built {len(structure["topics"])} API entries in each language at {root}')


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
