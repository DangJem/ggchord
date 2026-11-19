"""Browse the plotting inventory by its recorded coordinate system."""
import html
import json
import re
from . import ROOT
from .registry import canonical

# Order and representatives are editorial. Membership comes from the plot manifest.
GROUPS = {
    'chord': ('弦图', 'Chord diagrams', ['showcasechord']),
    'collinear': ('共线性图', 'Collinear maps', ['showcasecollinear']),
    'circular': ('单序列圆图', 'Single-sequence maps', ['firstcircle', 'showcaseplasmid']),
    'dotplot': ('序列比对点图', 'Alignment dotplots', ['viewdotplot']),
    'radial_synteny': ('径向共线性图', 'Radial synteny', ['viewradial']),
    'sankey': ('双列区间关系图', 'Sankey diagrams', ['foundationsankey']),
    'fusion': ('多视图组合', 'Multi-view compositions', ['fusionview']),
}


REPRESENTATIVES = {
    'showcasechord': ('四条基因组间的连接带', 'Ribbons across four genomes'),
    'showcasecollinear': ('在直线骨架上比较基因方向', 'Gene orientation on linear backbones'),
    'firstcircle': ('pUC19 的环形骨架', 'Circular backbone of pUC19'),
    'showcaseplasmid': ('pETDuet-1 的密集注释', 'Dense annotation on pETDuet-1'),
    'viewdotplot': ('比对区间的正向与反向斜率', 'Forward and reverse alignments'),
    'viewradial': ('用辐条比较多条序列', 'Compare sequences on radial spokes'),
    'foundationsankey': ('在左右两列追踪区间关系', 'Trace intervals between two columns'),
    'fusionview': ('连接圆图与线性面板', 'Connect circular and linear panels'),
}


def gallery(assets, lang, render_figure):
    def word(cn, en):
        return cn if lang == 'zh' else en

    web = ROOT / 'website'
    captions = json.loads((web / 'assets/captions.json').read_text())
    destinations = json.loads((web / 'figure-lessons.json').read_text())
    comparisons = {row['id'] for row in json.loads((web / 'parameter-comparisons.json').read_text())}
    grouped = {key: [] for key in GROUPS}
    for row in assets:
        grouped[row['coordinate']].append(row)

    def card(row):
        key = row['id']
        slug, _, anchor = destinations[key].partition('#')
        url = canonical(slug + '.html' + ('#' + anchor if anchor else ''))
        caption = captions[key][lang]
        title = html.escape(word(*REPRESENTATIVES[key]) if key in REPRESENTATIVES else
                            re.split(r'[:：。；]|; |[.] ', caption, maxsplit=1)[0])
        return (f'<article class="gallery-card" data-comparison="{str(key in comparisons).lower()}">'
                f'<h3 id="{key}">{title}</h3>{render_figure(key, url)}'
                f'<a href="{url}">{word("查看对应条目", "Open related entry")}</a></article>')

    options = ''.join(f'<option value="{key}">{word(cn, en)}</option>'
                      for key, (cn, en, _) in GROUPS.items())
    controls = ('<div class="gallery-controls" hidden>'
                f'<label>{word("图形", "Plot type")}<select name="family">'
                f'<option value="all">{word("所有图形", "All plot types")}</option>{options}</select></label>'
                f'<label>{word("显示", "Show")}<select name="scope">'
                f'<option value="featured">{word("代表图", "Selected figures")}</option>'
                f'<option value="all">{word("全部案例", "All figures")}</option>'
                f'<option value="comparison">{word("参数对照", "Parameter comparisons")}</option>'
                '</select></label><p class="gallery-count" role="status" aria-live="polite"></p></div>')
    groups = []
    for key, (cn, en, selected) in GROUPS.items():
        rows = {row['id']: row for row in grouped[key]}
        featured = ''.join(card(rows.pop(asset_id)) for asset_id in selected)
        remainder = ''
        if rows:
            label = word(f'展开其余 {len(rows)} 张图', f'Browse {len(rows)} more figures')
            remainder = (f'<details class="gallery-more"><summary>{label}</summary>'
                         '<div class="gallery-grid">' + ''.join(map(card, rows.values())) + '</div></details>')
        groups.append(f'<section class="gallery-family" data-family="{key}">'
                      f'<h2 id="gallery-{key}">{word(cn, en)}</h2>'
                      f'<div class="gallery-grid">{featured}</div>{remainder}</section>')
    return controls + '<div class="gallery-families">' + ''.join(groups) + '</div>'
