"""Shared bilingual presentation components; no scientific data transformations."""
import html,json,re,csv,struct
from functools import lru_cache
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent

@lru_cache(maxsize=1)
def assets():
    icon=json.loads((ROOT/'website/assets/icon.json').read_text())
    result={'logo':dict(file='logo.webp',width=icon['width'],height=icon['height'])}
    with (ROOT/'website/assets/plots.csv').open(newline='') as handle:
        for row in csv.DictReader(handle):
            path=ROOT/row['file']
            width,height=struct.unpack('>II',path.read_bytes()[16:24])
            result[row['id']]=dict(file=path.name,width=width,height=height)
    return result

@lru_cache(maxsize=1)
def lesson_urls():
    mapping=json.loads((ROOT/'website/figure-lessons.json').read_text())
    return {key:target.replace('#','.html#',1) for key,target in mapping.items()}


def motif(mode, css='mode-motif'):
    paths=json.loads((ROOT/'website/module-motifs.json').read_text())
    return '<svg class="'+css+'" viewBox="0 0 96 48" aria-hidden="true" focusable="false">'+paths[mode]+'</svg>'

def responsive_image(key,lang,base='../assets',*,eager=False,css='',alt=None):
    entry=assets()[key]
    return '<img class="'+css+'" src="'+base+'/'+entry['file']+'" alt="'+html.escape(alt or ('ggchord' if key=='logo' else key),quote=True)+'" width="'+str(entry['width'])+'" height="'+str(entry['height'])+'" decoding="async" '+('loading="eager" fetchpriority="high"' if eager else 'loading="lazy"')+'>'

def home_hero(lang,base):
    zh=lang=='zh'
    return ('<section class="home-hero"><div><p class="eyebrow">ggchord · R / ggplot2</p><h1 class="home-lead">'+
        ('<span class="hero-phrase">真实坐标，</span><span class="hero-phrase">四个绘图模块。</span>' if zh else 'Real coordinates. Four plotting modules.')+'</h1><p>'+
        ('比较多序列关系、组织单序列注释、阅读 HSP 双轴点图，并用真实源区间连接独立面板。' if zh else 'Compare sequence relationships, annotate one sequence, inspect HSP endpoints, and connect independent views through real source intervals.')+
        '</p><div class="home-actions"><a class="primary-action" href="installation.html">'+('从安装开始' if zh else 'Start with installation')+' →</a><a href="#choose">'+('选择绘图任务' if zh else 'Choose your task')+'</a></div></div>'+responsive_image('logo',lang,base,eager=True,css='brand-art')+'</section>')

def task_cards(lang):
    zh=lang=='zh'
    specs=[
      ('multi-sequence','多序列关系','Multi-sequence relationships','Chord · Collinear · Radial · Sankey','viewradial','序列 ID、完整长度；连接和特征按需加入。','Sequence IDs, complete lengths; optional links and features.'),
      ('circular','单序列圆图','Single-sequence circular maps','Circular','showcaseplasmid','一条序列长度与已有注释表；计算或搜索另需 DNA。','One sequence length and prepared annotations; computations need DNA.'),
      ('dotplot','HSP 点图','HSP Dotplot','Dotplot','viewdotplot','已有 HSP 端点与完整来源长度；不运行比对。','Existing HSP endpoints and complete lengths; no alignment run.'),
      ('fusion','多视图组合','Multi-view composition','Fusion','fusionview','独立有效的已注册图与明确来源区间；不是新的 Coord。','Independently valid registered plots and explicit source intervals.')]
    result='<div class="task-cards">'
    for slug,cn,en,views,fig,input_cn,input_en in specs:
        result+=f'<section class="task-card" data-module="{slug}"><a href="{slug}.html">'+responsive_image(fig,lang,css='module-preview')+f'<h3>{cn if zh else en}</h3></a><p>{views}</p><p class="minimum-input">'+(input_cn if zh else input_en)+f'</p><a href="{slug}.html">'+('进入模块' if zh else 'Explore module')+' →</a></section>'
    return result+'</div>'
