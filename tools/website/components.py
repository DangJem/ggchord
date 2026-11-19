"""Shared bilingual presentation components; no scientific data transformations."""
import html,json,re,csv,struct
from functools import lru_cache
from . import ROOT

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
    return {key:target.replace('#','.html#',1) if '#' in target else target+'.html' for key,target in mapping.items()}


def motif(mode, css='mode-motif'):
    paths=json.loads((ROOT/'website/module-motifs.json').read_text())
    return '<svg class="'+css+'" viewBox="0 0 96 48" aria-hidden="true" focusable="false">'+paths[mode]+'</svg>'

def responsive_image(key,lang,base='../assets',*,eager=False,css='',alt=None):
    entry=assets()[key]
    captions=json.loads((ROOT/'website/assets/captions.json').read_text())
    description=alt or ('ggchord' if key=='logo' else captions[key][lang])
    return '<img class="'+css+'" src="'+base+'/'+entry['file']+'" alt="'+html.escape(description,quote=True)+'" width="'+str(entry['width'])+'" height="'+str(entry['height'])+'" decoding="async" '+('loading="eager" fetchpriority="high"' if eager else 'loading="lazy"')+'>'

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
      ('dotplot','序列比对点图','Sequence alignment Dotplot','Dotplot','viewdotplot','已有 HSP 端点与完整来源长度；不运行比对。','Existing HSP endpoints and complete lengths; no alignment run.'),
      ('fusion','多视图组合','Multi-view composition','Fusion','fusionview','独立有效的已注册图与明确来源区间；不是新的 Coord。','Independently valid registered plots and explicit source intervals.')]
    result='<div class="task-cards">'
    for slug,cn,en,views,fig,input_cn,input_en in specs:
        result+=f'<section class="task-card" data-module="{slug}"><a href="{slug}.html">'+responsive_image(fig,lang,css='module-preview')+f'<h3>{cn if zh else en}</h3></a><p>{views}</p><p class="minimum-input">'+(input_cn if zh else input_en)+f'</p><a href="{slug}.html">'+('进入模块' if zh else 'Explore module')+' →</a></section>'
    return result+'</div>'

@lru_cache(maxsize=1)
def manual_figures():
    return {r['id']:r for r in json.loads((ROOT/'website/manual-figures.json').read_text())}

def manual_figure(key,lang,base='../assets'):
    r=manual_figures()[key];caption=html.escape(r['caption'][lang])
    target=base+'/'+key+'.png'
    code=html.escape(r['page'],quote=True)
    return (f'<figure class="documentation-figure" data-figure="{key}"><a class="figure-link" href="{target}" data-code-url="{code}">'
      f'<img src="{target}" alt="{caption}" width="{round(r["width"]*r["dpi"])}" height="{round(r["height"]*r["dpi"])}" loading="lazy" decoding="async"></a>'
      f'<figcaption>{caption} · <a href="{code}">'+('完整 R 源码' if lang=='zh' else 'Complete R source')+'</a></figcaption></figure>')

def manual_home(lang,base):
    zh=lang=='zh'
    content=json.loads((ROOT/'website/home.json').read_text())
    previews={row['id']:row for row in json.loads((ROOT/'website/home-previews.json').read_text())}
    complete_figures=manual_figures()
    escape=lambda value:html.escape(value,quote=True)
    def picture(key, title):
        row=previews[key]
        from .handbook import case_figure
        code=escape(complete_figures[row['full']]['page'])
        complete='handbook/'+case_figure(code.split('.')[0])+'.png'
        name=escape(row['name'])
        text=escape(title+' · '+row['name'])
        label=escape(('查看完整图 ' if zh else 'View complete figure for ')+row['name'])
        return (f'<figure class="home-preview"><a class="figure-link" href="{base}/{complete}" '
          f'data-code-url="{code}" aria-label="{label}"><img src="{base}/{key}.png" alt="{text}" '
          f'width="{round(row["width"]*row["dpi"])}" height="{round(row["height"]*row["dpi"])}" loading="lazy" decoding="async"></a>'
          f'<figcaption><a href="{code}">{name}</a></figcaption></figure>')
    cards=[]
    for module in content['modules']:
        title=module['title'][lang]
        purpose=escape(module['purpose'][lang])
        entry=escape(module['entry'])
        case=escape(module['case'])
        slug=escape(module['slug'])
        grid=' home-four' if len(module['figures'])==4 else ''
        art=''.join(picture(key, title) for key in module['figures'])
        cards.append(f'<section class="home-module" aria-labelledby="home-{slug}"><h2 id="home-{slug}">{escape(title)}</h2>'
          f'<p class="home-purpose">{purpose}</p><div class="home-module-art{grid}">{art}</div>'
          f'<div class="home-module-actions"><a class="module-entry" href="{entry}">'+('试画' if zh else 'Try it')+'</a>'
          f'<a class="home-source" href="{case}">'+('完整示例和代码' if zh else 'Full example and code')+'</a></div></section>')
    features=''.join('<a href="'+item['url']+'"><strong>'+escape(item['title'][lang])+'</strong><span>'+escape(item['description'][lang])+'</span></a>' for item in content['capabilities'])
    hero=content['hero']
    return ('<div class="home-shell"><section class="home-introduction"><div class="home-intro-copy"><h1>'+
        escape(hero['title'][lang])+'</h1><p class="home-description">'+escape(hero['description'][lang])+
        '</p><div class="home-actions"><a class="primary-action" href="installation.html">'+
        ('开始使用' if zh else 'Get started')+'</a><a class="secondary-action" href="gallery.html">'+
        ('浏览案例' if zh else 'Browse examples')+'</a></div><nav class="home-features" aria-label="'+
        ('文档入口' if zh else 'Documentation')+'">'+features+'</nav></div>'+
        responsive_image('logo',lang,base,eager=True,css='home-logo')+
        '</section><div class="home-capabilities" id="modules">'+''.join(cards)+'</div><p class="home-image-note">'+
        ('点击预览图可查看完整图。' if zh else 'Select a preview to open the complete figure.')+'</p></div>')
