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
    def figure(key):
        caption=html.escape(CAPTIONS[key][lang])
        hint='Click to enlarge' if lang=='en' else '点击放大查看'
        return (f'<figure class="documentation-figure"><a class="figure-link" href="{figures}/{key}.png" '
            f'title="{hint}"><img src="{figures}/{key}.png" alt="{caption}" loading="lazy"></a>'
            f'<figcaption>{caption} · <span class="figure-hint">{hint}</span></figcaption></figure>')
    text = re.sub(r'\{\{code:(\w+)\}\}', lambda m: '\n```r\n' + chunks[m[1]] + '\n```\n', text)
    text = re.sub(r'\{\{figure:(\w+)\}\}', lambda m: figure(m[1]), text)
    modes = '<div class="mode-grid">' + ''.join(
        f'<section class="mode-card"><h3 id="preview-{k}">{html.escape(CAPTIONS[k][lang])}</h3>'
        f'<img src="{figures}/{k}.png" alt="{html.escape(CAPTIONS[k][lang])}" loading="lazy">'
        f'<a href="modes.html#{k}">{"Code" if lang == "en" else "代码"}</a></section>'
        for k in ['chord','circular','collinear'])+'</div>'
    text = text.replace('{{modes}}', modes)
    migration = json.loads((WEB / 'migrations.json').read_text())
    rows = [[f'`{r[0]}`', f'`{r[1]}`', r[2], r[3 if lang == 'en' else 4]] for r in migration]
    text = text.replace('{{migration}}', table(
        ['Old (do not run)', 'Current', 'Removal / guard', 'Meaning'] if lang == 'en' else
        ['旧写法（勿执行）', '当前写法', '删除 / 守卫', '含义'], rows))
    text = text.replace('{{terms}}', table(['English', '现代汉语', 'API / units'],
        json.loads((WEB / 'terms.json').read_text())))
    with (WEB / 'assets/plots.csv').open() as f:
        assets = list(csv.DictReader(f))
    gallery = []
    for asset in assets:
        key = asset['id']
        gallery.append(f'### {CAPTIONS[key][lang]} {{#{key}}}\n\n{figure(key)}\n\n'
            f'`{asset["object"]}` · {asset["width"]} × {asset["height"]} in · '
            f'{asset["dpi"]} dpi · PNG · `{key}`')
    results_path = WEB/'assets/example-results.json'
    if results_path.is_file():
        results = json.loads(results_path.read_text())
        result_labels=json.loads((WEB/'assets/result-labels.json').read_text())
        text = re.sub(r'\{\{results:(\w+)\}\}', lambda m: table(
            ['Metric','Observed output'] if lang=='en' else ['项目','结果'],
            [[result_labels.get(k,{}).get(lang,k), ('依赖未安装，已跳过' if lang=='zh' and isinstance(v,str) and v.startswith('SKIP:') else v)] for k,v in results[m[1]].items()]),text)
    coverage=json.loads((WEB/'coverage.json').read_text())
    titles={p['slug']:p['title'][lang] for p in PAGES}
    article_slugs=['real-data','plasmids','modes','data','formats','alignments','parameters','ribbons','annotations','labels','tracks','homology','projection','views','styling','adjustments','output','migration','faq']
    text=text.replace('{{articles}}','\n'.join(f'- [{titles[k]}]({k}.html)' for k in article_slugs))
    text=text.replace('{{coverage}}',table(
        ['Capability','Read','Executable source / output'] if lang=='en' else ['功能','使用教程','代码与结果'],
        [[c['title'][lang],', '.join(f'[{titles[k]}]({k}.html)' for k in c['articles']),
          ', '.join('`'+k+'`' for k in c['snippets'])+(' · '+', '.join(c['figures']) if c['figures'] else '')] for c in coverage]))
    text = text.replace('{{gallery}}', '\n\n'.join(gallery)).replace('{{reference}}', reference)
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
