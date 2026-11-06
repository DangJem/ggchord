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
        ('真实坐标，三种视角。' if zh else 'Real coordinates. Three perspectives.')+'</h1><p>'+
        ('用弦图看序列联系，用共线性图比较结构，用圆图组织单序列注释。独立图层，同一套可复现的绘图语法。' if zh else 'Trace relationships with chord plots, compare structure with collinear views, and annotate a single sequence with circular maps. Independent layers. Reproducible figures.')+
        '</p><div class="home-actions"><a class="primary-action" href="installation.html">'+('从安装开始' if zh else 'Start with installation')+' →</a><a href="#choose">'+('选择绘图任务' if zh else 'Choose your task')+'</a></div></div>'+responsive_image('logo',lang,base,eager=True,css='brand-art')+'</section>')

def task_cards(lang):
    zh=lang=='zh';items=[('chord','弦图','Chord','比较多条序列中对应区间的分布与联系。','Compare corresponding intervals across multiple sequences.','序列 ID 与完整长度；添加连接需区间比对表。','Sequence IDs and full lengths; interval alignments for links.','chord'),('collinear','共线性图','Collinear','按共同尺度比较基因位置、方向与局部结构。','Compare gene positions, directions and local structure at a common scale.','序列 ID 与完整长度；可逐步加入基因和比对。','Sequence IDs and full lengths; add genes and alignments as layers.','collinear'),('plasmid','单序列圆图','Circular','组织一条序列的特征、标签、引物和定量轨道。','Organize features, labels, primers and quantitative tracks on one sequence.','一条序列的 ID 与完整长度；识别特征需要 DNA。','One sequence ID and full length; recognition needs DNA.','plasmids')]
    result='<div class="task-cards">'
    for mode,cn,en,desc_cn,desc_en,input_cn,input_en,case in items:
        result+=f'<section class="task-card" data-module="{mode}">'+motif(mode)+f'<h3>{cn if zh else en}</h3><p>{desc_cn if zh else desc_en}</p><p class="minimum-input">'+('最少输入：'+input_cn if zh else 'Minimum input: '+input_en)+f'</p><div class="task-links"><a href="{mode}-start.html">'+('第一张图' if zh else 'First plot')+f' →</a><a href="{case}.html">'+('完整案例' if zh else 'Complete case')+'</a></div></section>'
    return result+'</div>'
