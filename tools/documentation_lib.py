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

def snippets():
    text = (ROOT / 'tools/documentation-examples.R').read_text()
    return dict((name, code.strip()) for name, code in re.findall(
        r'^# ---- (\w+) ----\n(.*?)(?=^# ---- \w+ ----)', text, re.M | re.S))

def table(headers, rows):
    def cell(x):
        return str(x).replace('|', '\\|').replace('\n', ' ')
    return '\n'.join(['| ' + ' | '.join(map(cell, headers)) + ' |',
        '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
        ['| ' + ' | '.join(map(cell, row)) + ' |' for row in rows])

def expand(text, lang, *, figures='assets', reference=''):
    chunks = snippets()
    with (WEB / 'assets/plots.csv').open() as f:
        assets = list(csv.DictReader(f))
    dimensions={a['id']:(round(float(a['width'])*float(a['dpi'])),
        round(float(a['height'])*float(a['dpi']))) for a in assets}
    def preparation(match):
        keys=match[1].split(',')
        label='Preparation code — run once on this page' if lang=='en' else '准备代码 · 在本页运行一次'
        code='\n\n'.join(chunks[key] for key in keys)
        return '<details class="preparation"><summary>'+label+'</summary>\n\n```r\n'+code+'\n```\n\n</details>'
    text=re.sub(r'\{\{prepare:([\w,]+)\}\}',preparation,text)
    def figure(key):
        caption=html.escape(CAPTIONS[key][lang])
        hint='Click to enlarge' if lang=='en' else '点击放大查看'
        width,height=dimensions[key]
        return (f'<figure class="documentation-figure"><a class="figure-link" href="{figures}/{key}.png" '
            f'title="{hint}"><img src="{figures}/{key}.png" alt="{caption}" width="{width}" height="{height}" loading="lazy"></a>'
            f'<figcaption>{caption} · <span class="figure-hint">{hint}</span></figcaption></figure>')
    text = re.sub(r'\{\{code:(\w+)\}\}', lambda m: '\n```r\n' + chunks[m[1]] + '\n```\n', text)
    text = re.sub(r'\{\{figure:(\w+)\}\}', lambda m: figure(m[1]), text)
    modes = '<div class="mode-grid">' + ''.join(
        f'<section class="mode-card"><h3 id="preview-{k}">{html.escape(CAPTIONS[k][lang])}</h3>'
        f'<img src="{figures}/{k}.png" alt="{html.escape(CAPTIONS[k][lang])}" loading="lazy">'
        f'<a href="modes.html#{k}">{"Code" if lang == "en" else "代码"}</a></section>'
        for k in ['chord','collinear','circular'])+'</div>'
    text = text.replace('{{modes}}', modes)
    showcase_specs = [
        ('showcasechord','chord-data','Chord comparisons','多序列弦图'),
        ('showcasecollinear','collinear-data','Collinear comparisons','共线性图'),
        ('showcaseplasmid','plasmid-data','Plasmid maps','质粒图')]
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
    destinations={'showcasechord':'chord#compose','showcasecollinear':'collinear#compose','showcaseplasmid':'plasmids#complex','chord':'modes#chord','circular':'modes#circular','collinear':'modes#collinear',
        'first':'installation#first-plot','realchord':'real-data#chord','reallinear':'real-data#linear',
        'plasmid':'plasmids#input','homologs':'homology#membership','tracks':'tracks#domain',
        'multitrack':'tracks#multiple','density':'ribbons#density','pointlinks':'ribbons#points',
        'featurestack':'annotations#stack','decorations':'annotations#decoration','focus':'views#focus',
        'labels':'labels#layout','labelauto':'labels#comparison','labelradial':'labels#comparison',
        'labelarc':'labels#comparison','styling':'styling#roles','paf':'formats#zero-based'}
    destinations.update(dict(zip(
        ['learnseq','learnlinks','learngenes','learnlabels','linearseq','linearlinks','linearlabels','plasmidbackbone','plasmidfeatures','plasmidtext'],
        ['chord-build#backbones','chord-build#links','chord-build#genes','chord-build#labels','collinear-build#backbones','collinear-build#links','collinear-build#labels','plasmid-build#backbone','plasmid-build#features','plasmid-build#labels'])))
    for key in CAPTIONS:
        if key.startswith('catalog'): destinations[key] = 'plasmid-catalog#'+key
    for asset in assets:
        key = asset['id']
        page,anchor=destinations[key].split('#')
        label='Data, code and interpretation →' if lang=='en' else '数据、代码与读图说明 →'
        gallery.append(f'<section class="gallery-card"><h3 id="{key}">{html.escape(CAPTIONS[key][lang])}</h3>'
            f'{figure(key)}<p><a href="{page}.html#{anchor}">{label}</a></p>'
            f'<p class="asset-size"><code>{asset["object"]}</code> · {asset["width"]} × {asset["height"]} in · '
            f'{asset["dpi"]} dpi</p></section>')
    results_path = WEB/'assets/example-results.json'
    if results_path.is_file():
        results = json.loads(results_path.read_text())
        result_labels=json.loads((WEB/'assets/result-labels.json').read_text())
        text = re.sub(r'\{\{results:(\w+)\}\}', lambda m: table(
            ['Metric','Observed output'] if lang=='en' else ['项目','结果'],
            [[result_labels.get(k,{}).get(lang,k), ('依赖未安装，已跳过' if lang=='zh' and isinstance(v,str) and v.startswith('SKIP:') else v)] for k,v in results[m[1]].items()]),text)
    coverage=json.loads((WEB/'coverage.json').read_text())
    titles={p['slug']:p['title'][lang] for p in PAGES}
    navigation=json.loads((WEB/'navigation.json').read_text())
    article_groups=[]
    for group in navigation:
        entries=[slug for slug in group['pages'] if slug not in ['index','articles','reference','news']]
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
