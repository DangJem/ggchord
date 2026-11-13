"""Shared, dependency-free source expansion for the site and paired tutorials."""
import csv
import html
import json
import re
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / 'website'
PAGES = json.loads((WEB / 'pages.json').read_text())
CAPTIONS = json.loads((WEB / 'assets/captions.json').read_text())
RECIPES = json.loads((WEB / 'showcase-recipes.json').read_text())

def snippets():
    text = (ROOT / 'tools/documentation-examples.R').read_text()
    text += '\n' + (ROOT / 'tools/site-course-examples.R').read_text()
    return dict((name, code.strip()) for name, code in re.findall(
        r'^# ---- (\w+) ----\n(.*?)(?=^# ---- \w+ ----)', text, re.M | re.S))


def localize_r_comments(code, lang):
    if lang != 'zh':
        return code
    translations = json.loads((WEB / 'r-comments-zh.json').read_text())
    def translate(match):
        indent, comment = match.groups()
        if comment in translations:
            return indent + '# ' + translations[comment]
        if comment.startswith('----') or re.match(r'^(?:ggplot2::)?ggsave\(', comment):
            return match[0]
        if re.search('[A-Za-z]{3}', comment) and not re.search('[\u3400-\u9fff]', comment):
            raise ValueError('Missing Han code comment: ' + comment)
        return match[0]
    return re.sub(r'(?m)^([ \t]*)# ?(.*)$', translate, code)

def table(headers, rows):
    def cell(x):
        return str(x).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(map(cell, headers)) + ' |',
        '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
        ['| ' + ' | '.join(map(cell, row)) + ' |' for row in rows])

def expand(text, lang, *, figures='assets', reference=''):
    from site_components import home_hero, task_cards, responsive_image, lesson_urls
    text=text.replace('{{sitehero}}',home_hero(lang,figures)).replace('{{taskcards}}',task_cards(lang))
    chunks = {name: localize_r_comments(code, lang) for name, code in snippets().items()}
    with (WEB / 'assets/plots.csv').open() as f:
        assets = list(csv.DictReader(f))
    dimensions={a['id']:(round(float(a['width'])*float(a['dpi'])),
        round(float(a['height'])*float(a['dpi']))) for a in assets}
    def preparation(match):
        keys=match[1].split(',')
        return '\n'.join(
            '\n```r\n' + re.sub(r'\np_\w+$', '', chunks[key]) + '\n```\n'
            for key in keys)
    text=re.sub(r'\{\{prepare:([\w,]+)\}\}',preparation,text)
    def recipe(match):
        key = match[1]
        if key not in RECIPES:
            raise ValueError(f'Unknown showcase recipe {key}')
        name = f'prepare-{key}' + ('-zh' if lang == 'zh' else '') + '.R'
        label = 'Download the preparation script' if lang == 'en' else '下载准备脚本'
        return f'[{label}]({figures}/{name})\n\n```r\nsource("{name}")\n```\n'
    text = re.sub(r'\{\{recipe:(\w+)\}\}', recipe, text)
    def figure(key, source_url=None):
        caption=html.escape(CAPTIONS[key][lang])
        code_url=source_url or lesson_urls().get(key,'gallery.html')
        if figures!='../assets':
            hint='Click to enlarge' if lang=='en' else '点击放大查看'
            return (f'<figure class="documentation-figure"><a class="figure-link" href="{figures}/{key}.png" '
                f'title="{hint}"><img src="{figures}/{key}.png" alt="{caption}" width="{dimensions[key][0]}" height="{dimensions[key][1]}" loading="lazy"></a>'
                f'<figcaption>{caption} · <span class="figure-hint">{hint}</span></figcaption></figure>')
        hint='Open full-resolution figure' if lang=='en' else '打开高清原图'
        preview=responsive_image(key,lang,figures,alt=CAPTIONS[key][lang]) if figures=='../assets' else f'<img src="{figures}/{key}.png" alt="{caption}" width="{dimensions[key][0]}" height="{dimensions[key][1]}" loading="lazy">'
        return (f'<figure class="documentation-figure" data-figure="{key}"><a class="figure-link" href="{figures}/{key}.png" '
            f'title="{hint}" data-code-url="{code_url}">{preview}</a>'
            f'<figcaption>{caption} · <span class="figure-hint">{hint}</span></figcaption></figure>')
    previews = {a['object']: (float(a['width']), float(a['height'])) for a in assets}
    def preview_code(code):
        lines = code.splitlines()
        if lines[-1] in previews:
            obj = lines[-1]
            width, height = previews[obj]
            lines[-1] = f'view_ggchord({obj}, width = {width:g}, height = {height:g}, units = "in")'
        elif lines[-1] == 'label_plots$auto':
            lines[-1] = 'view_ggchord(label_plots$auto, width = 12, height = 10, units = "in")'
        elif re.fullmatch(r'p_\w+', lines[-1]):
            lines[-1] = f'view_ggchord({lines[-1]})'
        return '\n'.join(lines)

    def render_code(match):
        keys = match[1].split(',')
        parts = []
        for index, key in enumerate(keys):
            code = chunks[key]
            if index < len(keys) - 1:
                code = re.sub(r'\np_\w+$', '', code)
            else:
                code = preview_code(code)
            parts.append(code)
        return '\n```r\n' + '\n\n'.join(parts) + '\n```\n'
    text = re.sub(r'\{\{code:([\w,]+)\}\}', render_code, text)
    text = re.sub(r'\{\{figure:(\w+)\}\}', lambda m: figure(m[1]), text)
    modes = '<div class="mode-grid">' + ''.join(
        f'<section class="mode-card"><h3 id="preview-{k}">{html.escape(CAPTIONS[k][lang])}</h3>'
        f'<img src="{figures}/{k}.png" alt="{html.escape(CAPTIONS[k][lang])}" loading="lazy">'
        f'<a href="modes.html#{k}">{"Code" if lang == "en" else "代码"}</a></section>'
        for k in ['chord','collinear','circular'])+'</div>'
    text = text.replace('{{modes}}', modes)
    showcase_specs = [
        ('showcasechord','chord-start','Chord comparisons','多序列弦图'),
        ('showcasecollinear','collinear-start','Collinear comparisons','共线性图'),
        ('showcaseplasmid','plasmid-start','Plasmid maps','质粒图')]
    overview = '<div class="showcase-overview">' + ''.join(
        f'<section class="showcase-card"><a href="{page}.html">'
        f'<img src="{figures}/{key}.png" alt="{html.escape(CAPTIONS[key][lang])}" '
        f'width="{dimensions[key][0]}" height="{dimensions[key][1]}">'
        f'<strong>{title_en if lang=="en" else title_zh} →</strong></a></section>'
        for key,page,title_en,title_zh in showcase_specs)+'</div>'
    text = text.replace('{{showcases}}', overview)
    if '{{plasmidcatalog}}' in text:
        with (WEB/'assets/plasmid-catalog.csv').open() as source:
            catalog = list(csv.DictReader(source))
        text = text.replace('{{plasmidcatalog}}', table(
            ['Reference', 'Length (bp)', 'Features'] if lang == 'en' else ['参考序列', '长度（bp）', '识别特征数'],
            [[f'[{row["accver"]}](#{row["figure"]})', row['length'], row['features']] for row in catalog]))
    migration = json.loads((WEB / 'migrations.json').read_text())
    rows = [[f'`{r[0]}`', f'`{r[1]}`', r[2], r[3 if lang == 'en' else 4]] for r in migration]
    text = text.replace('{{migration}}', table(
        ['Old (do not run)', 'Current', 'Removal / guard', 'Meaning'] if lang == 'en' else
        ['旧写法（勿执行）', '当前写法', '删除 / 守卫', '含义'], rows))
    text = text.replace('{{terms}}', table(['English', '现代汉语', 'API / units'],
        json.loads((WEB / 'terms.json').read_text())))
    gallery = []
    destinations=json.loads((WEB/'figure-lessons.json').read_text())
    for asset in assets:
        key = asset['id']
        page,anchor=destinations[key].split('#')
        label='Related lesson →' if lang=='en' else '相关教程 →'
        source_label='Shared R source' if lang=='en' else '共享 R 源码'
        source_url=f'{figures}/documentation-examples.R' if figures=='../assets' else None
        source_link=f' · <a href="{source_url}">{source_label}</a>' if source_url else ''
        gallery.append(f'<section class="gallery-card"><h3 id="{key}">{html.escape(CAPTIONS[key][lang])}</h3>'
            f'{figure(key,source_url)}<p><a href="{page}.html#{anchor}">{label}</a>{source_link}</p>'
            f'<p class="asset-size"><code>{asset["object"]}</code> · {asset["width"]} × {asset["height"]} in · '
            f'{asset["dpi"]} dpi</p></section>')
    results_path = WEB/'assets/example-results.json'
    if results_path.is_file():
        results = json.loads(results_path.read_text())
        result_labels=json.loads((WEB/'assets/result-labels.json').read_text())
        text = re.sub(r'\{\{results:(\w+)\}\}', lambda m: table(
            ['Metric','Observed output'] if lang=='en' else ['项目','结果'],
            [[result_labels.get(k,{}).get(lang,k), ('依赖未安装，已跳过' if lang=='zh' and isinstance(v,str) and v.startswith('SKIP:') else v)] for k,v in results[m[1]].items()]),text)
    if '{{output:' in text or '{{console:' in text:
        outputs = json.loads((WEB / 'assets/lesson-outputs.json').read_text())
        text = re.sub(r'\{\{output:(\w+)\}\}', lambda m:
            '<div class="data-output" tabindex="0" role="region" aria-label="R output">\n\n' +
            table(outputs[m[1]]['columns'], outputs[m[1]]['rows']) + '\n\n</div>', text)
        text = re.sub(r'\{\{console:(\w+)\}\}', lambda m:
            '\n```text\n' + outputs[m[1]]['console'] + '\n```\n', text)
    coverage=json.loads((WEB/'coverage.json').read_text())
    titles={p['slug']:p['title'][lang] for p in PAGES}
    navigation=json.loads((WEB/'navigation.json').read_text())
    article_groups=[]
    for group in navigation:
        entries=[slug for slug in group['pages'] if slug not in ['index','about','reference']]
        if entries:
            article_groups.append('### '+group['title'][lang]+' {#path-'+group['id']+'}\n\n'+
                '\n'.join(f'- [{titles[k]}]({k}.html)' for k in entries))
    text=text.replace('{{articles}}','\n\n'.join(article_groups))
    text=text.replace('{{coverage}}',table(
        ['Capability','Read','Executable source / output'] if lang=='en' else ['功能','使用教程','代码与结果'],
        [[c['title'][lang],', '.join(f'[{titles[k]}]({k}.html)' for k in c['articles']),
          ', '.join('`'+k+'`' for k in c['snippets'])+(' · '+', '.join(c['figures']) if c['figures'] else '')] for c in coverage]))
    text = text.replace('{{gallery}}', '<div class="gallery-grid">'+''.join(gallery)+'</div>').replace('{{reference}}', reference)
    return text

def category(name):
    if name.startswith(('coord_', 'ggchord_api_', 'ggchord_layout_')) or name in ['ggchord']:
        return 'construction'
    if name.startswith(('geom_track_', 'scale_track_', 'position_track')):
        return 'tracks'
    if name.startswith('geom_'):
        return 'layers'
    if name.startswith(('scale_', 'guide_', 'theme_', 'position_')):
        return 'style'
    if name.startswith(('read_', 'as_ggchord', 'validate_', 'clean_', 'prepare_', 'ggchord_data_')):
        return 'data'
    if name.startswith(('find_', 'filter_restriction', 'primer_catalog')):
        return 'recognition'
    if name.startswith(('focus_', 'segment_', 'map_', 'build_', 'project_', 'link_', 'filter_', 'merge_', 'deduplicate_', 'bundle_', 'stat_')):
        return 'workflow'
    if name.startswith(('get_', 'export_', 'view_', 'diagnose_', 'optimize_')):
        return 'output'
    return 'other'

CATEGORIES = {
 'construction': ('Construction and contracts', '构建与契约'),
 'data': ('Input and validation', '输入与校验'),
 'layers': ('Geometry and labels', '几何与标签'),
 'style': ('Positions, scales, themes and guides', '位置、标度、主题与图例'),
 'tracks': ('Quantitative tracks', '定量轨道'),
 'recognition': ('Sequence recognition', '序列识别'),
 'workflow': ('Projection and comparison workflows', '投影与比较工作流'),
 'output': ('Output and diagnostics', '输出与诊断'),
 'other': ('Methods, teaching data and upstream reexports', '方法、教学数据与上游重导出')}
