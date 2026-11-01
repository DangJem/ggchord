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
        return f'![{CAPTIONS[key][lang]}]({figures}/{key}.png)'
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
