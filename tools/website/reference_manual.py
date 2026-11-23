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
        upstream='Native ggplot2 tools',
        archived='Historical API address', current='Open the current API reference'),
    'zh': dict(reference='API 参考',
        usage='调用形式', arguments='参数', flexible='【支持灵活参数】',
        purpose='说明', value='返回值', details='行为规则', errors='错误与限制',
        roles='数据角色与视觉映射', examples='示例', related='相关接口',
        source='原生 ggplot2', download='下载示例', filter='按函数名或用途筛选',
        console='控制台输出', figure='示例绘图',
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
    families = read(SOURCE / 'families.json')
    members = [name for names in families.values() for name in names]
    if sorted(members) != sorted(topics):
        raise ValueError('Each API must belong to exactly one family')
    categories = {name: group for group, record in structure['groups'].items() for name in record['topics']}
    for group, record in structure['groups'].items():
        shortcuts = record.get('shortcuts', [])
        if len(shortcuts) != len(set(shortcuts)) or any(
                name not in topics or categories[name] == group for name in shortcuts):
            raise ValueError(f'{group}: invalid cross-category reference entries')
    for primary, names in families.items():
        if primary != names[0] or len({categories[name] for name in names}) != 1:
            raise ValueError(f'{primary}: invalid family landing page or mixed categories')
    structure['families'] = families
    previous = structure['previous_family_members']
    if any(head not in topics or any(member not in topics or member == head or member in families.get(head, [])
                                    for member in members)
           for head, members in previous.items()):
        raise ValueError('Former family members must have moved to another article')
    articles = structure['articles']
    if any(head not in families or len(families[head]) < 2 or not set(spec['related']) <= set(topics)
           for head, spec in articles.items()):
        raise ValueError('Shared articles must cover a family and link to known topics')
    for lang, locale in locales.items():
        if set(locale['families']) != set(families) or not all(locale['families'].values()):
            raise ValueError(f'{lang}: missing API family titles')
        if not locale['description'].strip():
            raise ValueError(f'{lang}: missing reference introduction')
        if set(locale['groups']) != set(structure['groups']) or set(locale['topics']) != set(topics):
            raise ValueError(f'{lang}: categories or topics missing')
        if set(locale['articles']) != set(articles):
            raise ValueError(f'{lang}: shared articles missing')
        for head, article in locale['articles'].items():
            members = families[head]
            arguments = {key for name in members for key in topics[name]['defaults']}
            if (set(article['examples']) != set(members) or
                    set(article['arguments']) != set(locales['en']['articles'][head]['arguments']) or
                    not set(article['arguments']) <= arguments or
                    not all(article[field].strip() for field in ('intro', 'roles', 'details', 'example_intro')) or
                    not all(note[field].strip() for note in article['examples'].values()
                            for field in ('title', 'description')) or
                    not all(entry.strip() for entry in article['arguments'].values())):
                raise ValueError(f'{lang}/{head}: incomplete or mismatched shared article')
            for key in arguments:
                entries = []
                for name in members:
                    part, body = topics[name], locale['topics'][name]
                    entries.extend((part['defaults'][key], param['format'], param['description'])
                                   for keys, param in zip(part['arguments'], body['arguments']) if key in keys)
                if len(set(entries)) > 1 and key not in article['arguments']:
                    raise ValueError(f'{lang}/{head}: reconcile shared argument {key} in the article')
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
    # Historical topic URLs resolve to their section on the family page.
    text = re.sub(r'\]\(#api-([a-z0-9_]+)\)', r'](\1.html)', text)
    return subprocess.check_output(
        [pandoc, '--from=markdown+pipe_tables', '--to=html5', '--wrap=none'],
        input=text, text=True)


def section(id, title, body, level=2):
    return f'<section aria-labelledby="{id}"><h{level} id="{id}">{html.escape(title)}</h{level}>{body}</section>'


def snippet(markup, kind, caption):
    # Keep the prepared token spans through pkgdown's article rendering pass.
    return (f'<details class="code-disclosure" data-code-kind="{kind}" open>'
            f'<summary>{html.escape(caption)}</summary><pre class="downlit"><code>{markup}'
            '</code></pre></details>')


def topic_link(topic, owners):
    return f'<a href="{owners[topic]}.html#api-{topic}">{code(topic + "()")}</a>'


def invocation(topic):
    name = topic['usage'].split('(', 1)[0]
    lines = [f'  {key} = {value}' if value else f'  {key}'
             for key, value in topic['defaults'].items()]
    return name + '(\n' + ',\n'.join(lines) + '\n)' if lines else name + '()'


def member_prefix(name, primary):
    return '' if name == primary else name + '-'


def parameter_anchor(key):
    return 'para-ellipsis' if key == '...' else 'parameter-' + key


def anchors(ids):
    return ''.join(f'<span id="{id}"></span>' for id in ids)


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


def prune_example_downloads(folder, keep):
    for pattern in ('*.R', '*-[0-9][0-9].png'):
        for path in folder.glob(pattern):
            if path.name not in keep:
                path.unlink()


def render_examples(assets, structure):
    scripts = prepare_examples(structure)
    with TemporaryDirectory(prefix='ggchord-reference-examples-') as temporary:
        staging = Path(temporary)
        inputs = staging / 'inputs.json'
        generated = staging / 'output'
        write_json(inputs, dict(topics={
            name: dict(signature=invocation(topic), example=scripts[name])
            for name, topic in structure['topics'].items()}))
        subprocess.run(['Rscript', '--vanilla', 'tools/website/render-reference-examples.R',
                        str(inputs), str(generated)], cwd=ROOT, check=True)
        results = read(generated / 'results.json')
        folder = assets / 'reference-examples'
        for name, source in scripts.items():
            write(folder / f'{name}.R', source)
        figures = {figure for entry in results['topics'].values() for figure in entry['figures']}
        prune_example_downloads(folder, figures | {f'{name}.R' for name in scripts})
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


def example_content(output, name, titles):
    rendered = output['topics'][name]
    return (snippet(rendered['example'], 'r', titles['examples']) +
            f'<a href="../../assets/reference-examples/{name}.R" download>{titles["download"]}</a>' +
            example_output(rendered, name, titles, output['width'], output['height']))


def example_search(url, name, title, anchor, lang, note):
    return dict(url=f'{url}#{anchor}', kind='tool', title=f'{name}() · {title}',
                keywords=f'{name} {LABELS[lang]["examples"]}', summary=note, text=note,
                queries=[f'{name} {word}' for word in ('example', 'examples', '示例')])


def bookmark_moves(structure, locale, primary, owners):
    # Match member bookmarks from the old landing page to their current
    # targets. The old landing page can now be a function alias itself.
    moves = {}
    for name in structure['previous_family_members'].get(primary, []):
        part, body = structure['topics'][name], locale['topics'][name]
        prefix = member_prefix(name, primary)
        target = owners[name]
        target_prefix = member_prefix(name, target)
        moves['api-' + name] = f'{target}.html#api-{name}'
        sections = ['purpose', 'description', 'usage', 'value', 'details', 'examples', 'related']
        if part['arguments']:
            sections.append('arguments')
        if part['roles'] or target in structure['articles']:
            sections.append('roles')
        if body['errors'] or target in structure['articles']:
            sections.append('errors')
        fragments = ['ref-' + key for key in sections]
        for key in part['defaults']:
            fragments.append(parameter_anchor(key))
            if key != '...':
                fragments.append('arg-' + key)
        moves.update({prefix + anchor: f'{target}.html#{target_prefix}{anchor}' for anchor in fragments})
    retired = structure.get('retired_family_members', {}).get(primary, [])
    if retired:
        migrations = read(ROOT / 'website/page-migrations.json')
        for name in retired:
            mapping = migrations[f'reference/{name}.html']
            moves['api-' + name] = mapping[''].removeprefix('reference/')
            moves.update({f'{name}-{anchor}': dest.removeprefix('reference/')
                          for anchor, dest in mapping.items() if anchor})
            moves['example-' + name] = mapping['ref-examples'].removeprefix('reference/')
    return moves


def render_bookmark_moves(moves):
    if not moves:
        return ''
    return ('<script data-reference-bookmarks>(function redirectReferenceBookmark(){'
            f'const targets={json.dumps(moves, separators=(",", ":"))};'
            'function redirect(){const hash=location.hash.slice(1);'
            'if(Object.hasOwn(targets,hash))location.replace(targets[hash]);}'
            'redirect();addEventListener("hashchange",redirect);})();</script>')


def render_shared_article(primary, structure, locale, lang, owners, output, flexible, pandoc):
    members = structure['families'][primary]
    article = locale['articles'][primary]
    topics = structure['topics']
    titles = LABELS[lang]
    url = f'reference/{primary}.html'
    prefixes = [member_prefix(name, primary) for name in members]
    search = []
    signatures, examples = [], []
    arguments = {}

    for name, prefix in zip(members, prefixes):
        shared = topics[name]
        text = locale['topics'][name]
        signatures.append(f'<div data-api-signature="{name}" id="api-{name}">' +
                          (anchors([prefix + 'ref-usage']) if prefix else '') +
                          snippet(output['topics'][name]['signature'], 'signature', name + '()') + '</div>')
        note = article['examples'][name]
        example_id = prefix + 'ref-examples' if prefix else 'example-' + name
        examples.append(section(example_id, note['title'], markdown(note['description'], pandoc) +
                                example_content(output, name, titles), level=3))
        search.append(example_search(url, name, note['title'], example_id, lang, note['description']))
        search.append(dict(url=f'{url}#api-{name}', title=f'{name}() — {text["title"]}',
                           kind='function', summary=text['intro'], text=article['intro'],
                           keywords=' '.join([name, locale['families'][primary], *shared['defaults']]),
                           queries=[name, name + '()']))
        for keys, param in zip(shared['arguments'], text['arguments']):
            for key in keys:
                arguments.setdefault(key, {})[name] = param

    rows = []
    for key, uses in arguments.items():
        param = next(iter(uses.values()))
        description = markdown(article['arguments'].get(key, param['description']), pandoc)
        summary = ' '.join(html.unescape(re.sub(r'<[^>]*>', ' ', description)).split())
        term_ids = []
        kinds = set()
        for name in uses:
            prefix = member_prefix(name, primary)
            anchor = prefix + parameter_anchor(key)
            term_ids.append(anchor)
            if key != '...':
                term_ids.append(prefix + 'arg-' + key)
            if kind := flexible.get(name, {}).get(key):
                kinds.add(kind)
            search.append(dict(url=f'{url}#{anchor}', kind='parameter', title=f'{name}() · {key}',
                               keywords=f'{name} {key}', text=summary, summary=summary,
                               queries=[f'{name} {key}']))
        if len(kinds) > 1:
            raise ValueError(f'{primary}: incompatible flexible parameter formats for {key}')
        badge = ''
        if kinds:
            fragment = FLEXIBLE_SECTIONS[next(iter(kinds))]
            badge = (f'<a class="flexible-parameter" href="../parameters.html#{fragment}">'
                     f'{titles["flexible"]}</a>')
        canonical, *aliases = term_ids
        rows.append(f'<dt><span class="argument-term">{anchors(aliases)}'
                    f'<code id="{canonical}" class="argument-anchor">{html.escape(key)}</code>'
                    f'{badge}</span></dt><dd>{description}</dd>')

    def common(key, content):
        aliases = anchors(prefix + 'ref-' + key for prefix in prefixes[1:])
        return section('ref-' + key, titles[key], aliases + content)

    intro_aliases = anchors(prefix + 'ref-' + key for prefix in prefixes for key in ('purpose', 'description'))
    details_aliases = anchors(prefix + 'ref-' + key for prefix in prefixes for key in ('value', 'errors'))
    related = structure['articles'][primary]['related']
    body = (
        intro_aliases + markdown(article['intro'], pandoc) +
        section('ref-usage', titles['usage'], ''.join(signatures)) +
        common('arguments', f'<dl class="reference-arguments">{"".join(rows)}</dl>') +
        common('roles', markdown(article['roles'], pandoc)) +
        common('details', details_aliases + markdown(article['details'], pandoc)) +
        section('ref-examples', titles['examples'], markdown(article['example_intro'], pandoc) + ''.join(examples)) +
        common('related', '<ul>' + ''.join(f'<li>{topic_link(name, owners)}</li>' for name in related) + '</ul>')
    )
    return f'<div data-api-article="{primary}">{body}</div>', search


def render_articles(root, pages, navigation_site):
    # Use the same pkgdown shell and navigation finalizer as the other chapters.
    with TemporaryDirectory(prefix='ggchord-reference-') as temporary:
        staging = Path(temporary)
        payload = staging / 'pages.json'
        write_json(payload, pages)
        # pkgdown protects files made by other generators. Render into a fresh
        # directory so old standalone pages cannot survive a successful build.
        options = [f'--site={staging}', f'--payload={payload}']
        for tool, flags in (('render-website.R', ['--part=articles']),
                            ('finalize-website.R', [f'--navigation-site={navigation_site}'])):
            result = subprocess.run(['Rscript', '--vanilla', f'tools/{tool}', *flags, *options],
                                    cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if result.returncode:
                raise RuntimeError(result.stdout.rstrip())
        for lang, entries in pages.items():
            rendered = staging / lang
            target = root / lang
            target.mkdir(parents=True, exist_ok=True)
            for dependency in ('deps', 'lightswitch.js', 'pkgdown.js'):
                asset = rendered / dependency
                if asset.is_dir():
                    for source in (p for p in asset.rglob('*') if p.is_file()):
                        path = target / source.relative_to(rendered)
                        if not (navigation_site / lang / source.relative_to(rendered)).exists():
                            path.parent.mkdir(parents=True, exist_ok=True)
                            shutil.copy2(source, path)
                elif not (navigation_site / lang / dependency).exists():
                    shutil.copy2(asset, target / dependency)
            for page in entries:
                relative = page['slug'] + '.html'
                path = target / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(rendered / relative, path)


def reference_index(structure, locale, lang, owners):
    groups = locale['groups']
    titles = LABELS[lang]
    families = structure['families']

    def family_index(primary):
        terms = ' · '.join(topic_link(name, owners) for name in families[primary])
        return (f'<dt class="col-sm-5 mb-sm-2">{link(primary + ".html", locale["families"][primary])}</dt>'
                f'<dd class="col-sm-7 mb-2">{terms}</dd>')

    contents = ''.join(
        '<section class="reference-group">' +
        ''.join(f'<span id="{alias}"></span>' for alias in group['previous_anchors']) +
        f'<h2 id="family-{key}">{html.escape(groups[key]["title"])}</h2>'
        f'<p>{html.escape(groups[key]["intro"])}</p><dl class="row">' +
        ''.join(family_index(primary) for primary in dict.fromkeys(owners[name] for name in group['topics'])) +
        '</dl>' +
        ('<p class="reference-shortcuts">' + ' · '.join(topic_link(name, owners) for name in group['shortcuts']) +
         '</p>' if group.get('shortcuts') else '') + '</section>'
        for key, group in structure['groups'].items())
    field = ('<div class="reference-filter" hidden><label for="reference-filter">'
             f'{titles["filter"]}</label><input id="reference-filter" type="search" autocomplete="off">'
             '<output id="reference-filter-count" aria-live="polite"></output></div>')
    upstream = '<ul>' + ''.join(
        f'<li>{link("https://ggplot2.tidyverse.org/reference/" + href, name)}</li>'
        for name, href in UPSTREAM) + '</ul>'
    body = (f'<p class="lead">{html.escape(locale["description"])}</p>' +
            field + contents + section('upstream', titles['upstream'], upstream))
    return dict(slug='reference/index', title=titles['reference'], body=body)


def prepare(root, pandoc):
    """Prepare articles, search rows and evaluated examples for either build mode."""
    structure, locales = load()
    assets = root / 'assets'
    output = render_examples(assets, structure)
    search_by_language = {}
    flexible = read(ROOT / 'website/flexible-parameters.json')
    former = structure['previous_family_members']
    revised = set(former) | {name for names in former.values() for name in names}
    revised.update(name for head in structure['articles'] for name in structure['families'][head])
    pages = {lang: [] for lang in LOCALES}
    for lang, locale in locales.items():
        topics, groups = locale['topics'], locale['groups']
        titles = LABELS[lang]
        search = []
        intro = locale['description']
        families = structure['families']
        owners = {name: primary for primary, names in families.items() for name in names}
        pages[lang].append(reference_index(structure, locale, lang, owners))
        search.append(dict(url='reference/index.html',title=titles['reference'],kind='tool',
                           summary=intro, text=' '.join(groups[k]['intro'] for k in groups),keywords='API reference'))
        for primary, members in families.items():
            moves = render_bookmark_moves(bookmark_moves(structure, locale, primary, owners))
            if primary in structure['articles']:
                body, rows = render_shared_article(primary, structure, locale, lang, owners, output, flexible, pandoc)
                pages[lang].append(dict(slug=f'reference/{primary}', title=locale['families'][primary],
                                        body=body + moves))
                search.extend(rows)
                continue
            sections = []
            for name in members:
                prefix = member_prefix(name, primary)
                url = f'reference/{primary}.html'
                shared = structure['topics'][name]
                text = topics[name]
                native = (f'<p>{link("https://ggplot2.tidyverse.org/reference/" + name + ".html", titles["source"])}</p>'
                          if shared['native'] else '')
                bodies = [native]
                def add(key, fragment):
                    if fragment:
                        bodies.append(section(prefix + 'ref-' + key, titles[key], fragment, level=3))
                add('purpose', f'<span id="{prefix}ref-description"></span>' + markdown(text['intro'], pandoc))
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
                        anchor = prefix + parameter_anchor(key)
                        if key in seen:
                            names.append(code(key))
                            continue
                        seen.add(key)
                        alias = f'<span id="{prefix}arg-{key}"></span>' if key != '...' else ''
                        badge = ''
                        kind = flexible.get(name, {}).get(key)
                        if kind or re.search(r'\bflexible\b', original['format']):
                            fragment = FLEXIBLE_SECTIONS.get(kind, 'scope')
                            badge = (f'<a class="flexible-parameter" href="../parameters.html#{fragment}">'
                                     f'{titles["flexible"]}</a>')
                        names.append(f'<span class="argument-term">{alias}'
                                     f'<code id="{anchor}" class="argument-anchor">{html.escape(key)}</code>'
                                     f'{badge}</span>')
                        hit = dict(url=f'{url}#{anchor}', kind='parameter',
                                   title=f'{name}() · {key}', keywords=f'{name} {key}',
                                   text=f'{param["format"]} {summary}', summary=summary)
                        if name in revised:
                            hit['queries'] = [f'{name} {key}']
                        search.append(hit)
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
                add('examples', example_note + example_content(output, name, titles))
                if name in revised:
                    search.append(example_search(url, name, titles['examples'], prefix + 'ref-examples',
                                                 lang, text['example_setup'] or text['intro']))
                add('related', '<ul>' + ''.join(f'<li>{topic_link(target, owners)}</li>' for target in shared['related']) + '</ul>')
                sections.append(f'<section data-api="{name}"><h2 id="api-{name}">{code(name + "()")}</h2>'
                                + ''.join(bodies) + '</section>')
                row = dict(url=f'{url}#api-{name}', title=f'{name}() — {text["title"]}',
                           kind='function', summary=text['intro'], text=text['details'],
                           keywords=' '.join([name, locale['families'][primary], *sorted(seen)]))
                if name in revised:
                    row['queries'] = [name, name + '()']
                search.append(row)
            pages[lang].append(dict(slug=f'reference/{primary}', title=locale['families'][primary],
                                    body=''.join(sections) + moves))
        search_by_language[lang] = search
    version = re.search(r'^Version: (.+)$', (ROOT / 'DESCRIPTION').read_text(), re.M)[1]
    write_json(root / 'api-reference-manifest.json', {
        'version': version, 'source': structure['reference'], 'functions': len(structure['topics']) - len(NATIVE),
        'native_facets': len(NATIVE), 'categories': len(structure['groups']), 'locales': list(LOCALES),
        'families': structure['families'], 'shared_articles': list(structure['articles']),
        'topics': list(structure['topics']), 'deployment': 'local output',
    })
    write_json(root / 'reference-output.json', {
        'status': 'passed', 'topics': {name: {'console': entry['console'], 'figures': entry['figures']}
                                     for name, entry in output['topics'].items()},
    })
    return pages, search_by_language


def write_reference_redirect(out, path, destination, languages):
    mapping = destination if isinstance(destination, dict) else {'': destination}
    destination = mapping['']
    labels = ({'en': 'Version selection', 'zh': '版本选择'} if destination == 'migration.html#api'
              else {'en': 'Current documentation', 'zh': '现行文档'})
    links = []
    for lang in languages:
        relative = os.path.relpath(out/lang/destination, path.parent)
        tag = 'zh-Hans' if lang == 'zh' else 'en'
        links.append('<a lang="' + tag + '" hreflang="' + tag + '" href="' +
                     html.escape(relative, quote=True) + '">' + labels[lang] + '</a>')
    targets = {fragment: os.path.relpath(out/languages[0]/target, path.parent)
               for fragment, target in mapping.items()}
    script = ('<script>const targets=' + json.dumps(targets, separators=(',', ':')) + ';'
              'const key=location.hash.slice(1);'
              'location.replace(Object.hasOwn(targets,key)?targets[key]:targets[""]);</script>')
    tag = 'zh-Hans' if languages == ('zh',) else 'en'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('<!doctype html><html lang="' + tag + '"><head><meta charset="utf-8">'
                    '<title>' + labels[languages[0]] + ' · ggchord</title>' + script +
                    '</head><body><main>' + ' · '.join(links) + '</main></body></html>')


def write_family_redirects(root, families):
    # Old per-function links occur in published articles and user bookmarks.
    # A hash on those URLs belonged to the individual function's article.
    structure = read(SOURCE / 'structure.json')
    owners = {name: primary for primary, names in families.items() for name in names}
    for lang in LOCALES:
        locale = read(SOURCE / f'{lang}.json')
        for primary, members in families.items():
            for name in members[1:]:
                fallback = primary + '.html#api-' + name
                moves = bookmark_moves(structure, locale, name, owners)
                script = ('const old = location.hash.slice(1);'
                          f'const former = {json.dumps(moves, separators=(",", ":"))};'
                          f'const local = old === "example-{name}" ? "ref-examples" : old;'
                          f'const anchor = /^(ref-|arg-|parameter-|para-)/.test(local) ? "{name}-" + local : "api-{name}";'
                          f'location.replace(Object.hasOwn(former, old) ? former[old] : "{primary}.html#" + anchor);')
                title = html.escape(name + '() · ggchord')
                tag = 'zh-Hans' if lang == 'zh' else 'en'
                marker = ' data-reference-bookmarks' if moves else ''
                write(root/lang/'reference'/f'{name}.html',
                      f'<!doctype html><html lang="{tag}"><head><meta charset="utf-8">'
                      f'<title>{title}</title><script{marker}>{script}</script><noscript>'
                      f'<meta http-equiv="refresh" content="0;url={fallback}"></noscript></head>'
                      f'<body><a href="{fallback}">{name}()</a></body></html>')


def build(root, write_package_download):
    root = root.resolve()
    pandoc = shutil.which('pandoc')
    if not pandoc:
        raise RuntimeError('Pandoc must be on PATH to build the API reference')
    # Build and render before replacing any files in the user's local preview.
    with TemporaryDirectory(prefix='ggchord-api-update-') as temporary:
        stage = Path(temporary)
        subprocess.run(['Rscript', '--vanilla', 'tools/website/parameter-directory.R', str(ROOT),
                        str(stage/'parameter-directory.json')], cwd=ROOT, check=True)
        pages, searches = prepare(stage, pandoc)
        write_package_download(stage/'assets')
        render_articles(stage, pages, navigation_site=root)
        manifest = read(stage/'api-reference-manifest.json')
        write_family_redirects(stage, manifest['families'])
        for lang, search in searches.items():
            search_file = root/lang/'search.json'
            previous = read(search_file) if search_file.exists() else []
            retained = [row for row in previous if not row['url'].startswith('reference/')]
            write_json(stage/lang/'search.json', retained + search)
            if not (root/lang/'index.html').exists():
                write(stage/lang/'index.html', redirect('reference/index.html'))
        # Keep unlocalized entry points, including their original fragments.
        routes_file = root/'legacy-routes.json'
        routes = read(routes_file) if routes_file.is_file() else {}
        for topic in ['index', *manifest['topics']]:
            name = f'reference/{topic}.html'
            routes[name] = 'en/' + name
            write(stage/name, redirect('../en/' + name))
        for name, mapping in read(ROOT/'website/page-migrations.json').items():
            if not name.startswith('reference/'):
                continue
            routes[name] = 'en/' + name
            for languages in [('en',), ('zh',), LOCALES]:
                prefix = languages[0] if len(languages) == 1 else ''
                destination = mapping['']
                target = Path(languages[0]) / destination.split('#', 1)[0]
                if not (stage/target).is_file() and not (root/target).is_file():
                    write_reference_redirect(stage, stage/prefix/name, 'reference/index.html', languages)
                    continue
                write_reference_redirect(stage, stage/prefix/name, mapping, languages)
        write_json(stage/'legacy-routes.json', routes)
        if not (root/'index.html').exists():
            write(stage/'index.html', redirect('en/reference/index.html'))
        for name in ('site.css', 'site.js', 'favicon.png', 'logo.webp'):
            if not (root/'assets'/name).exists():
                shutil.copy2(ASSETS/name, stage/'assets'/name)
        for source in stage.rglob('*'):
            if source.is_file():
                dest = root/source.relative_to(stage)
                if dest.is_file() and dest.read_bytes() == source.read_bytes():
                    continue
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, dest)
        generated = stage/'assets/reference-examples'
        prune_example_downloads(root/'assets/reference-examples', {path.name for path in generated.iterdir()})
    print(f'Built {len(manifest["topics"])} API entries in {len(manifest["families"])} families per language at {root}')


def redirect(destination):
    href = html.escape(destination, quote=True)
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>ggchord API reference</title>'
            f'<script>location.replace({json.dumps(destination)} + location.hash)</script>'
            f'<noscript><meta http-equiv="refresh" content="0;url={href}"></noscript></head>'
            f'<body><a href="{href}">API reference</a></body></html>\n')


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
