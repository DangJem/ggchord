#!/usr/bin/env python3
"""Build paired package-native documentation; never execute R examples/vignettes."""
import argparse
import hashlib
import html
import json
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path
from site_search import complete_examples, search_records
from documentation_lib import ROOT, WEB, PAGES, CAPTIONS, RECIPES, expand, snippets, localize_r_comments

class PageText(__import__('html.parser',fromlist=['HTMLParser']).HTMLParser):
    def __init__(self,text):
        super().__init__();self.parts=[];self.stack=[];self.ignored=0;self.feed(text)
    def handle_starttag(self,tag,attrs):
        values=dict(attrs)
        excluded=bool(self.ignored) or tag in ['nav','script','style','pre','summary','h1'] or bool(
            {'chapter-drawer','page-outline','candidate-note','lesson-stage','breadcrumb','example-actions','section-number'} & set(values.get('class','').split()))
        if tag not in ['area','base','br','col','embed','hr','img','input','link','meta','param','source','track','wbr']:
            self.stack.append((tag,excluded));self.ignored+=int(excluded)
    def handle_endtag(self,tag):
        for i in range(len(self.stack)-1,-1,-1):
            if self.stack[i][0]==tag:
                self.ignored-=sum(int(skip) for _,skip in self.stack[i:])
                del self.stack[i:];break
    def handle_data(self,data):
        if not self.ignored:self.parts.append(data)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--publication',action='store_true')
    args=parser.parse_args()
    out=ROOT/'.site'
    if out.is_symlink():raise SystemExit('Refusing a symlink output directory')
    if out.exists():shutil.rmtree(out)
    out.mkdir();assets=out/'assets';assets.mkdir()
    for p in (WEB/'assets').iterdir():
        if p.is_file():shutil.copy2(p,assets/p.name)
    complete_examples(assets)
    for p in (ROOT/'man/figures').iterdir():
        if p.is_file() and (p.stem in CAPTIONS or p.stem=='logo'):shutil.copy2(p,assets/p.name)
    for name in ['documentation-examples.R','documentation-plots.R']:
        shutil.copy2(ROOT/'tools'/name,assets/name)
    shutil.copy2(ROOT/'NEWS.md',assets/'NEWS.md')
    for name in ['MIGRATION.md', 'MIGRATION-Hans.md']:
        shutil.copy2(ROOT/name, assets/name)
    shutil.copy2(WEB/'ASSET-SOURCES.md',assets/'ASSET-SOURCES.md')
    shutil.copytree(ROOT/'R',assets/'source')
    shutil.copytree(ROOT/'tools/teaching-audit',assets/'source/tools/teaching-audit')
    inputs=sorted([*ROOT.glob('examples/fasta/*.fna'),*ROOT.glob('examples/gff3/*.gff3'),
        *ROOT.glob('examples/blastn/*.o7'),*ROOT.glob('examples/blastn-standard/*'),
        *ROOT.glob('examples/format-inputs/*'),*ROOT.glob('examples/protein/*'),*ROOT.glob('examples/alignment/*'),*ROOT.glob('examples/plasmid/*.fna')])
    # New lessons need their actual auxiliary inputs/scripts in the download,
    # not merely code links into a maintainer checkout. Preserve historical
    # provenance in sources.json while excluding its inactive data file.
    historical = {item['file'] for item in json.loads(
        (ROOT/'examples/genome-alignments/sources.json').read_text())
        if item.get('status') == 'historical-only'}
    inputs += [p for p in ROOT.glob('examples/genome-alignments/*')
        if p.is_file() and p.name not in historical]
    inputs += [p for p in ROOT.glob('examples/cds/*') if p.is_file()]
    for name in ['m2-detail.R', 'm3-nested.R', 'origin-detail.R', 'fusion-auto.R', 'fusion-routing.R']:
        inputs.append(ROOT/'examples'/name)
    inputs = sorted(set(inputs))
    input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    (assets/'input-hashes.json').write_text(json.dumps(input_hashes,indent=2)+'\n')
    shared_code = snippets()
    recipes = {}
    for name, keys in RECIPES.items():
        for lang in ['en', 'zh']:
            filename = f'prepare-{name}' + ('-zh' if lang == 'zh' else '') + '.R'
            recipes[filename] = '\n\n'.join(localize_r_comments(shared_code[key], lang) for key in keys) + '\n'
            (assets / filename).write_text(recipes[filename])
    with zipfile.ZipFile(assets/'tutorial-inputs.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
        for p in inputs:archive.write(p,p.relative_to(ROOT))
        archive.write(ROOT/'tools/documentation-examples.R','tools/documentation-examples.R')
        for name, code in recipes.items(): archive.writestr(name, code)
        archive.writestr('input-hashes.json',json.dumps(input_hashes,indent=2)+'\n')
    pandoc_dir=subprocess.check_output(['Rscript','-e','cat(rmarkdown::find_pandoc()$dir)'],cwd=ROOT,text=True).strip()
    payload={lang:[] for lang in ['en','zh']}
    for lang in payload:
        for page in PAGES:
            if page['slug']=='reference':continue
            md=expand((WEB/lang/(page['slug']+'.md')).read_text(),lang,figures='../assets')
            body=subprocess.check_output([str(Path(pandoc_dir)/'pandoc'),'--from=markdown+pipe_tables+header_attributes',
                '--to=html5','--wrap=none'],input=md,text=True)
            headings=re.findall(r'<h[23] id="([^"]+)"[^>]*>(.*?)</h[23]>',body,re.S)
            toc=''.join(f'<li><a href="#{a}">{re.sub("<[^>]+>","",t)}</a></li>' for a,t in headings)
            sidebar=('<h2>本页目录</h2>' if lang=='zh' else '<h2>On this page</h2>')+'<ul class="toc-list">'+toc+'</ul>'
            sidebar+=('<h2>项目资料</h2>' if lang=='zh' else '<h2>Project</h2>')+'<ul><li><a href="https://github.com/DangJem/ggchord">GitHub</a></li><li><a href="about.html#modules">'+('章节索引' if lang=='zh' else 'Chapter index')+'</a></li><li>MIT</li></ul>'
            payload[lang].append(dict(slug=page['slug'],title=page['title'][lang],body=body,sidebar=sidebar))
    (out/'page-content.json').write_text(json.dumps(payload,ensure_ascii=False)+'\n')
    # Locales write separate destinations and share only immutable source files.
    jobs=[]
    for lang in ['en','zh']:
        with (out/(lang+'-build.log')).open('w') as log:
            jobs.append(subprocess.Popen(['Rscript','tools/render-website.R',lang],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT))
    statuses=[job.wait() for job in jobs]
    combined='\n'.join((out/(lang+'-build.log')).read_text() for lang in ['en','zh'])
    (out/'build.log').write_text(combined)
    if any(statuses):raise SystemExit(combined[-7000:])
    # Post-process using the same XML library used by pkgdown; no hand-built page shell.
    subprocess.run(['Rscript','tools/finalize-website.R', 'publication' if args.publication else 'preview'],cwd=ROOT,check=True)
    # Keep project presentation files identical to the source digested by finalization.
    for name in ['site.css','site.js']: shutil.copy2(WEB/'assets'/name,assets/name)
    for lang in ['en','zh']:
        def redirect(old,target):
            path=out/lang/old;path.parent.mkdir(parents=True,exist_ok=True)
            dest=os.path.relpath(out/lang/target,path.parent)
            path.write_text(f'<!doctype html><html lang="{"zh-Hans" if lang=="zh" else "en"}"><head><meta charset="utf-8"><meta http-equiv="refresh" content="0;url={dest}"><title>ggchord</title></head><body><a href="{dest}">Continue / 继续</a></body></html>')
        redirect('reference.html','reference/index.html')
        index=[]
        for path in sorted((out/lang).rglob('*.html')):
            if path.name=='reference.html':continue
            content=path.read_text()
            title=html.unescape(re.search(r'<title>(.*?)</title>',content,re.S)[1]).split(' • ')[0]
            main=re.search(r'<main\b[^>]*>(.*?)</main>',content,re.S)
            text=' '.join(PageText(main[1] if main else content).parts)
            index.append(dict(url=str(path.relative_to(out/lang)),title=title,text=re.sub(r'\s+',' ',text)))
        (out/lang/'search.json').write_text(json.dumps(search_records(out,lang,index),ensure_ascii=False)+'\n')
    # Stable historical paths remain usable after a future, separately authorized cutover.
    topics=[p.stem for p in (ROOT/'man').glob('*.Rd')]
    routes={'index.html':'en/index.html','articles/ggchord_vignette.html':'en/real-data.html',
        'articles/ggchord_guide_hans.html':'zh/real-data.html','reference/index.html':'en/reference/index.html','news/index.html':'assets/NEWS.md'}
    routes.update({'reference/'+s+'.html':'en/reference/'+s+'.html' for s in topics})
    for old,target in routes.items():
        path=out/old;path.parent.mkdir(parents=True,exist_ok=True);dest=os.path.relpath(out/target,path.parent)
        extra='<script>try{if(localStorage.getItem("ggchord-language")==="zh")location.replace("zh/index.html");else location.replace("en/index.html")}catch(e){location.replace("en/index.html")}</script>' if old=='index.html' else f'<meta http-equiv="refresh" content="0;url={dest}">'
        path.write_text(f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>ggchord</title>{extra}</head><body><a href="{dest}">English</a> · <a href="{os.path.relpath(out/"zh/index.html",path.parent)}">现代汉语</a></body></html>')
    retired=[]
    for rel in json.loads((WEB/'historical-pages.json').read_text()):
        if rel in routes:continue
        path=out/rel;path.parent.mkdir(parents=True,exist_ok=True)
        en=os.path.relpath(out/'en/reference/index.html',path.parent);zh=os.path.relpath(out/'zh/reference/index.html',path.parent)
        path.write_text(f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Historical documentation</title></head><body><h1>Historical documentation / 历史文档</h1><p>This historical URL has no exact current equivalent. 此旧地址没有一一对应的当前页面。</p><a href="{en}">Current reference</a> · <a href="{zh}">现代汉语参考</a></body></html>')
        retired.append(rel)
    for name,data in [('legacy-routes',routes),('retired-routes',retired)]:
        (out/(name+'.json')).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    files=[p for base in [WEB,ROOT/'R',ROOT/'man',ROOT/'tools',ROOT/'vignettes'] for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    files += inputs + sorted((ROOT/'data').glob('*.rda'))
    files += [ROOT/p for p in ['DESCRIPTION','NAMESPACE','NEWS.md','README.md','README-Hans.md','MIGRATION.md','MIGRATION-Hans.md','_pkgdown.yml','.github/workflows/documentation.yaml']]
    manifest=dict(version=re.search(r'^Version: (.+)$',(ROOT/'DESCRIPTION').read_text(),re.M)[1],api_inventory='current NAMESPACE/formals; tools/teaching-audit/teaching-summary-v'+'.'.join(re.search(r'^Version: (.+)$',(ROOT/'DESCRIPTION').read_text(),re.M)[1].split('.')[:3])+'.json',generator='pkgdown + paired Markdown + shared static assets',
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        deployment='prepared for authorized deployment; builder does not deploy' if args.publication else 'local preview; not deployed',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (out/'build-manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    (out/'page-content.json').unlink()
    shutil.rmtree(out/'reference-zh')
    print(f'Built {len(PAGES)} paired entries with pkgdown and {len(topics)} reference topics; no deployment, examples or vignette execution.')

if __name__=='__main__':main()
