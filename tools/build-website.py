#!/usr/bin/env python3
"""Build a local-only static site from paired Markdown and canonical Rd."""
import argparse
import hashlib
import html
import json
import re
import shutil
import subprocess
from pathlib import Path
from documentation_lib import ROOT, WEB, PAGES, CATEGORIES, CAPTIONS, category, expand

def run(*args):
    return subprocess.check_output(args, cwd=ROOT, text=True).strip()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='.site')
    parser.add_argument('--publication', action='store_true', help='Build for an explicitly authorized website deployment; does not deploy or publish an R release')
    args = parser.parse_args()
    out = (ROOT / args.output).resolve()
    # Only the conventional generated output can be replaced automatically.
    if out != ROOT / '.site':
        raise SystemExit('Use the repository .site output; it is separate from frozen docs/.')
    if out.is_dir():
        shutil.rmtree(out)
    out.mkdir()
    assets = out / 'assets'
    assets.mkdir()
    for path in (WEB / 'assets').iterdir():
        if path.is_file():
            shutil.copy2(path, assets / path.name)
    for path in (ROOT / 'man/figures').glob('*.png'):
        if path.stem in CAPTIONS:
            shutil.copy2(path, assets / path.name)
    shutil.copy2(ROOT / 'man/figures/logo.svg', assets / 'logo.svg')
    shutil.copy2(ROOT / 'man/figures/logo.png', assets / 'logo.png')
    shutil.copy2(ROOT / 'tools/documentation-examples.R', assets / 'documentation-examples.R')
    shutil.copy2(ROOT / 'NEWS.md', assets / 'NEWS.md')
    refs = out / '_reference'
    subprocess.run(['Rscript', 'tools/render-reference.R', str(refs)], cwd=ROOT, check=True)
    topics = json.loads((refs / 'index.json').read_text())
    for topic in topics:
        if isinstance(topic['aliases'], str):
            topic['aliases'] = [topic['aliases']]
    pandoc = Path(run('Rscript', '-e', 'cat(rmarkdown::find_pandoc()$dir)')) / 'pandoc'
    version = re.search(r'^Version: (.+)$', (ROOT / 'DESCRIPTION').read_text(), re.M)[1]
    commit = run('git', 'rev-parse', 'HEAD')
    dirty = bool(run('git', 'status', '--porcelain'))
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for base in [WEB, ROOT / 'R', ROOT / 'man', ROOT / 'tools']
              for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    for name in ['DESCRIPTION','NAMESPACE','NEWS.md','README.md','README-Hans.md','MIGRATION.md','MIGRATION-Hans.md','_pkgdown.yml','.github/workflows/documentation.yaml']:
        path = ROOT/name
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    provenance = dict(version=version, api_candidate='0.16.4', commit=commit,
        dirty=dirty, deployment=('prepared for authorized website deployment; consult Pages run for actual status' if args.publication else 'local preview; not deployed'), source_sha256=hashes)
    (out / 'build-manifest.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2)+'\n')
    labels = {
      'en': dict(skip='Skip to content', search='Search documentation', theme='Theme', toc='On this page',
          menu='Browse', home='Home', status='Local development documentation · not deployed', api='Canonical English API text',
          source='Source', result='Search results', empty='No matching pages', open='Code', build='Build source'),
      'zh': dict(skip='跳至正文', search='搜索文档', theme='主题', toc='本页目录', menu='浏览目录', home='首页',
          status='本地开发文档 · 尚未部署', api='API 正文为完整英语原文', source='源码', result='搜索结果', empty='没有匹配页面',
          open='代码', build='构建来源')}
    if args.publication:
        labels['en']['status'] = 'Development documentation · unreleased API candidate'
        labels['zh']['status'] = '开发文档 · 尚未正式发布的 API 候选'
    search = {'en': [], 'zh': []}
    def write(lang, slug, title, body, *, topic=None):
        L = labels[lang]
        depth = '../' if topic else ''
        base = '../' + depth
        counterpart = f'{base}{"zh" if lang == "en" else "en"}/{slug}.html'
        nav = ''.join(f'<a {"aria-current=page" if p["slug"] == slug else ""} href="{depth}{p["slug"]}.html">{html.escape(p["title"][lang])}</a>' for p in PAGES)
        headings = []
        def heading(m):
            level, attrs, text = m.groups()
            name = re.search(r'\bid=["\']([^"\']+)', attrs)
            anchor = name[1] if name else 'section-' + str(len(headings)+1)
            headings.append((anchor, re.sub('<[^>]*>', '', text)))
            return f'<h{level} id="{html.escape(anchor)}">{text}</h{level}>'
        body = re.sub(r'<h([23])([^>]*)>(.*?)</h\1>', heading, body, flags=re.S)
        body = re.sub(r'(<pre[^>]*>.*?</pre>)', lambda m: '<details open class="code"><summary>'+L['open']+'</summary>'+m[1]+'</details>', body, flags=re.S)
        toc = ''.join(f'<a href="#{html.escape(a)}">{html.escape(t)}</a>' for a,t in headings)
        source = topic['source'] if topic else f'website/{lang}/{slug}.md'
        clean = re.sub('<[^>]+>', ' ', body)
        # A single local index contains complete text, including han substrings.
        search[lang].append(dict(url=slug+'.html', title=title, text=html.unescape(re.sub(r'\s+', ' ', clean))))
        brand_art = ('<picture class="brand-art"><source srcset="'+base+'assets/logo.webp" type="image/webp"><img src="'+base+'assets/logo.png" alt="ggchord" width="1254" height="1254"></picture>') if slug == 'index' else ''
        canonical_notice = f'<p class="notice">{L["api"]} · <a href="{depth}reference.html">{L["menu"]}</a></p>' if topic else ''
        page = f'''<!doctype html>
<html lang="{'zh-Hans' if lang == 'zh' else 'en'}"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)} · ggchord</title>
<meta name="description" content="{html.escape(title)} — ggchord {version}">
<link rel="icon" href="{base}assets/favicon.png" type="image/png">
<link rel="alternate" hreflang="{'en' if lang=='zh' else 'zh-Hans'}" href="{counterpart}">
<link rel="stylesheet" href="{base}assets/site.css"><script src="{base}assets/site.js" defer></script></head>
<body data-lang="{lang}" data-search="{depth}search.json" data-empty="{L['empty']}">
<a class="skip" href="#content">{L['skip']}</a>
<header><a class="brand" href="{depth}index.html"><img src="{base}assets/favicon.png" alt="" width="44" height="44">ggchord <span>{version}</span></a>
<div class="header-controls"><a class="language" href="{counterpart}" lang="{'en' if lang=='zh' else 'zh-Hans'}">{'English' if lang=='zh' else '现代汉语'}</a>
<button id="site-theme" type="button" aria-label="{L['theme']}">{L['theme']}</button></div></header>
<div class="status">{L['status']} · API 0.16.4</div>
<div class="layout"><aside class="sidebar"><details class="navigation" open><summary>{L['menu']}</summary><nav aria-label="{L['menu']}">{nav}</nav></details>
<label for="site-search">{L['search']}</label><input id="site-search" type="search" autocomplete="off" placeholder="accver, track, 标签…">
<div id="site-results" role="region" aria-label="{L['result']}" aria-live="polite"></div></aside>
<main id="content" tabindex="-1"><h1>{html.escape(title)}</h1>{canonical_notice}{brand_art}{body}
<footer>{L['source']}: <code>{source}</code> ·
<a href="{base}build-manifest.json">{L['build']}</a> · {commit[:7]}{' + local changes' if dirty else ''} · MIT</footer></main>
<aside class="toc"><nav aria-label="{L['toc']}"><strong>{L['toc']}</strong>{toc}</nav></aside></div></body></html>'''
        dest = out / lang / (slug + '.html')
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(page)
    base_packages = {'base','compiler','datasets','graphics','grDevices','grid','methods','parallel','splines','stats','stats4','tcltk','tools','utils'}
    def external_help(match):
        branch = 'R' if match[1] in base_packages else 'CRAN'
        return f'href="https://search.r-project.org/{branch}/refmans/{match[1]}/html/{match[2]}{match[3] or ""}"'
    for lang in ['en', 'zh']:
        refsections = []
        for cat, names in CATEGORIES.items():
            matches = [t for t in topics if category(t['slug']) == cat]
            if matches:
                refsections.append(f'## {names[0 if lang=="en" else 1]} {{#{cat}}}\n')
                refsections += [f'- ['+ ' · '.join('`'+a+'`' for a in t['aliases']) + f'](reference/{t["slug"]}.html) — {t["title"]}' for t in matches]
        for p in PAGES:
            text = expand((WEB/lang/(p['slug']+'.md')).read_text(), lang, figures='../assets', reference='\n'.join(refsections))
            body = subprocess.check_output([str(pandoc), '--from=markdown+pipe_tables+header_attributes', '--to=html5', '--wrap=none'], input=text, text=True)
            write(lang, p['slug'], p['title'][lang], body)
        for t in topics:
            full = (refs/(t['slug']+'.html')).read_text()
            body = full.split('<main>', 1)[1].split('</main>', 1)[0]
            body = re.sub(r'<table style="width: 100%;">.*?</table>', '', body, count=1, flags=re.S)
            body = re.sub(r'<hr>.*$', '', body, flags=re.S)
            body = body.replace('href="00Index.html"', 'href="../reference.html"')
            # R's external package links refer to an R installation tree; use public manuals.
            body = re.sub(r'href="(?:\.\./)+([^/]+)/html/([^"#]+)(#[^"]*)?"',
                external_help, body)
            write(lang, 'reference/'+t['slug'], t['title'], body, topic=t)
        (out/lang/'search.json').write_text(json.dumps(search[lang],ensure_ascii=False)+'\n')
    # Frozen historical URLs get paired landing routes only after a future cutover.
    routes = {'index.html': 'en/index.html', 'articles/ggchord_vignette.html': 'en/modes.html',
              'articles/ggchord_guide_hans.html': 'zh/modes.html', 'reference/index.html': 'en/reference.html',
              'news/index.html': 'en/news.html'}
    for t in topics:
        routes['reference/'+t['slug']+'.html'] = 'en/reference/'+t['slug']+'.html'
    for old, target in routes.items():
        path = out/old; path.parent.mkdir(parents=True,exist_ok=True)
        dest = Path(__import__('os').path.relpath(out/target,path.parent)).as_posix()
        path.write_text(f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ggchord documentation</title><meta http-equiv="refresh" content="0;url={dest}"></head><body><a href="{dest}">Continue / 继续阅读</a></body></html>')
    retired = []
    for old in sorted((ROOT/'docs').rglob('*.html')):
        relative = str(old.relative_to(ROOT/'docs'))
        if relative in routes:
            continue
        path = out/relative
        path.parent.mkdir(parents=True, exist_ok=True)
        depth = '../' * (len(path.relative_to(out).parts)-1)
        target = depth+'en/reference.html'
        migration = depth+'en/migration.html'
        path.write_text(f'<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Historical documentation · ggchord</title></head><body><h1>Historical documentation / 历史文档</h1><p>This URL belonged to the frozen site and has no one-to-one current API page. 此地址属于冻结旧站，没有一一对应的当前 API 页面。</p><p><a href="{target}">Current reference / 当前函数参考</a> · <a href="{migration}">Migration / 迁移</a></p></body></html>')
        retired.append(relative)
    (out/'retired-routes.json').write_text(json.dumps(retired,ensure_ascii=False,indent=2)+'\n')
    (out/'legacy-routes.json').write_text(json.dumps(routes,ensure_ascii=False,indent=2)+'\n')
    shutil.rmtree(refs)
    print(f'Built {len(PAGES)*2} guide pages and {len(topics)*2} reference pages in {out}. No deployment performed.')

if __name__ == '__main__':
    main()
