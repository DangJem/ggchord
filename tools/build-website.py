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
from documentation_lib import ROOT, WEB, PAGES, CAPTIONS, expand

class PageText(__import__('html.parser',fromlist=['HTMLParser']).HTMLParser):
    def __init__(self,text):
        super().__init__();self.parts=[];self.feed(text)
    def handle_data(self,data):self.parts.append(data)

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
    for p in (ROOT/'man/figures').iterdir():
        if p.is_file() and (p.stem in CAPTIONS or p.stem=='logo'):shutil.copy2(p,assets/p.name)
    for name in ['documentation-examples.R','documentation-plots.R']:
        shutil.copy2(ROOT/'tools'/name,assets/name)
    shutil.copy2(ROOT/'NEWS.md',assets/'NEWS.md')
    shutil.copy2(WEB/'ASSET-SOURCES.md',assets/'ASSET-SOURCES.md')
    shutil.copytree(ROOT/'R',assets/'source')
    inputs=sorted([*ROOT.glob('examples/fasta/*.fna'),*ROOT.glob('examples/gff3/*.gff3'),
        *ROOT.glob('examples/blastn/*.o7'),*ROOT.glob('examples/blastn-standard/*'),
        *ROOT.glob('examples/format-inputs/*'),ROOT/'examples/plasmid/pUC19.fna'])
    input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    (assets/'input-hashes.json').write_text(json.dumps(input_hashes,indent=2)+'\n')
    with zipfile.ZipFile(assets/'tutorial-inputs.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
        for p in inputs:archive.write(p,p.relative_to(ROOT))
        archive.write(ROOT/'tools/documentation-examples.R','tools/documentation-examples.R')
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
            sidebar+=('<h2>项目资料</h2>' if lang=='zh' else '<h2>Project</h2>')+'<ul><li><a href="https://github.com/DangJem/ggchord">GitHub</a></li><li><a href="articles.html">'+('全部专题' if lang=='zh' else 'All articles')+'</a></li><li>MIT</li></ul>'
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
        (out/lang/'search.json').write_text(json.dumps(index,ensure_ascii=False)+'\n')
    # Stable historical paths remain usable after a future, separately authorized cutover.
    topics=[p.stem for p in (ROOT/'man').glob('*.Rd')]
    routes={'index.html':'en/index.html','articles/ggchord_vignette.html':'en/real-data.html',
        'articles/ggchord_guide_hans.html':'zh/real-data.html','reference/index.html':'en/reference/index.html','news/index.html':'en/news.html'}
    routes.update({'reference/'+s+'.html':'en/reference/'+s+'.html' for s in topics})
    for old,target in routes.items():
        path=out/old;path.parent.mkdir(parents=True,exist_ok=True);dest=os.path.relpath(out/target,path.parent)
        extra='<script>try{if(localStorage.getItem("ggchord-language")==="zh")location.replace("zh/index.html");else location.replace("en/index.html")}catch(e){location.replace("en/index.html")}</script>' if old=='index.html' else f'<meta http-equiv="refresh" content="0;url={dest}">'
        path.write_text(f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>ggchord</title>{extra}</head><body><a href="{dest}">English</a> · <a href="{os.path.relpath(out/"zh/index.html",path.parent)}">现代汉语</a></body></html>')
    retired=[]
    for old in sorted((ROOT/'docs').rglob('*.html')):
        rel=str(old.relative_to(ROOT/'docs'))
        if rel in routes:continue
        path=out/rel;path.parent.mkdir(parents=True,exist_ok=True)
        en=os.path.relpath(out/'en/reference/index.html',path.parent);zh=os.path.relpath(out/'zh/reference/index.html',path.parent)
        path.write_text(f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Historical documentation</title></head><body><h1>Historical documentation / 历史文档</h1><p>This historical URL has no exact current equivalent. 此旧地址没有一一对应的当前页面。</p><a href="{en}">Current reference</a> · <a href="{zh}">现代汉语参考</a></body></html>')
        retired.append(rel)
    for name,data in [('legacy-routes',routes),('retired-routes',retired)]:
        (out/(name+'.json')).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    files=[p for base in [WEB,ROOT/'R',ROOT/'man',ROOT/'tools',ROOT/'vignettes'] for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    files += inputs + [ROOT/'data'/n for n in ['seq_data_example.rda','gene_data_example.rda','ribbon_data_example.rda','plasmid_example_pUC19.rda','plasmid_example_pETDuet_1.rda']]
    files += [ROOT/p for p in ['DESCRIPTION','NAMESPACE','NEWS.md','README.md','README-Hans.md','MIGRATION.md','MIGRATION-Hans.md','_pkgdown.yml','.github/workflows/documentation.yaml']]
    manifest=dict(version='0.16.5',api_candidate='0.16.4',generator='pkgdown + paired Markdown + shared static assets',
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        deployment='prepared for authorized deployment; builder does not deploy' if args.publication else 'local preview; not deployed',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    (out/'build-manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    (out/'page-content.json').unlink()
    shutil.rmtree(out/'reference-zh')
    print(f'Built {len(PAGES)} paired entries with pkgdown and {len(topics)} reference topics; no deployment, examples or vignette execution.')

if __name__=='__main__':main()
