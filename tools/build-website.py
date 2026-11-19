#!/usr/bin/env python3
"""Build paired package-native documentation; evaluate reference examples, never knit vignettes."""
import argparse
import hashlib
import html
import json
import os
import re
import shutil
import subprocess
import time
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from site_cache import Cache, runtime, key as cache_key, shell_inputs, reference_inputs, reference_cache_allowed
from website.search import complete_examples, search_records
from website.source_links import tracked_paths
from website import ROOT
from website.handbook import WEB, expand, manifest
from website.registry import pages

PAGES = pages()

class PageText(__import__('html.parser',fromlist=['HTMLParser']).HTMLParser):
    def __init__(self,text):
        super().__init__();self.parts=[];self.stack=[];self.ignored=0;self.feed(text)
    def handle_starttag(self,tag,attrs):
        values=dict(attrs)
        excluded=bool(self.ignored) or tag in ['nav','script','style','pre','summary','h1'] or bool(
            {'chapter-drawer','page-outline','candidate-note','breadcrumb','example-actions'} & set(values.get('class','').split()))
        if tag not in ['area','base','br','col','embed','hr','img','input','link','meta','param','source','track','wbr']:
            self.stack.append((tag,excluded));self.ignored+=int(excluded)
    def handle_endtag(self,tag):
        for i in range(len(self.stack)-1,-1,-1):
            if self.stack[i][0]==tag:
                self.ignored-=sum(int(skip) for _,skip in self.stack[i:])
                del self.stack[i:];break
    def handle_data(self,data):
        if not self.ignored:self.parts.append(data)

def public_files(paths):
    """Select regular, repository-owned downloads without following local symlinks."""
    tracked = tracked_paths()
    for path in sorted(set(paths)):
        relative = path.relative_to(ROOT)
        if relative.as_posix() not in tracked or any(part.startswith('.') for part in relative.parts):
            continue
        if path.resolve() != path:
            raise ValueError(f'Expected a regular repository file at {relative}')
        if path.is_file():
            yield path


def copy_source_downloads(assets):
    source = ROOT/'R'
    for path in public_files(source.rglob('*')):
        target = assets/'source'/path.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)


def render_handbook_examples(assets):
    output = assets / 'handbook'
    subprocess.run(['Rscript', '--vanilla', 'tools/website/draw-handbook.R',
                    str(ROOT), str(output)], cwd=ROOT, check=True)
    hashes = {}
    for name, spec in manifest().items():
        if spec['kind'] == 'plot':
            script = WEB / 'examples' / (name + '.R')
            png = output / (name + '.png')
            hashes[name] = dict(source_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
                                png_sha256=hashlib.sha256(png.read_bytes()).hexdigest())
    (output / 'rendered.json').write_text(json.dumps(hashes, indent=2) + '\n')


def render_source_snapshot():
    """Inputs consumed by rendering; audit outputs can still be written."""
    roots = [WEB, ROOT/'R', ROOT/'man', ROOT/'data', ROOT/'vignettes']
    files = {p for root in roots for p in root.rglob('*') if p.is_file()
             and '__pycache__' not in p.parts and
             not any(part.startswith('.') for part in p.relative_to(ROOT).parts)}
    files.update(public_files(ROOT.glob('*.md')))
    files.update(ROOT/row['file'] for row in json.loads((WEB/'home-previews.json').read_text()))
    files.update(p for p in (ROOT/'tools/website').iterdir() if p.suffix in {'.py', '.R'})
    files.update(p for p in (ROOT/'inst/extdata').rglob('*') if p.is_file())
    files.update(public_files((ROOT/'examples').rglob('*')))
    files.update(ROOT/name for name in [
        'DESCRIPTION', 'NAMESPACE', 'LICENSE', '_pkgdown.yml',
        'tools/build-website.py', 'tools/site_cache.py', 'tools/documentation_lib.py',
        'tools/render-website.R', 'tools/finalize-website.R', 'tools/reference-translations.py',
        'tools/website-dependencies.R', 'tools/check-documentation.py',
        'tools/sync-website-languages.py', 'man/figures/logo.png', 'man/figures/logo.svg',
        'data-raw/website-data-sources.tsv', 'data-raw/flagship-data-sources.tsv'])
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(files)}


def assert_render_source_snapshot(snapshot, input_hashes):
    current_inputs = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                      for name in input_hashes}
    if render_source_snapshot() != snapshot or current_inputs != input_hashes:
        raise RuntimeError('Website source inputs changed during the build; '
                           'retain partial files but do not mark the build complete.')


def audit_reference_output(out):
    """Require native output beside every topic's code, rather than code alone."""
    records = {}
    errors = []
    for lang in ('en', 'zh'):
        rows = []
        for rd in sorted((ROOT/'man').glob('*.Rd')):
            path = out/lang/'reference'/(rd.stem+'.html')
            text = path.read_text() if path.exists() else ''
            at = text.find('id="ref-examples"')
            examples = text[at:] if at >= 0 else ''
            figures = re.findall(r'<img\b[^>]*src="([^"]+)"', examples)
            text_output = 'class="r-out' in examples
            if not figures and not text_output:
                errors.append(f'{lang}/{rd.stem}: missing evaluated output')
            if 'class="r-err' in examples:
                errors.append(f'{lang}/{rd.stem}: example execution error')
            if not re.search(r'<h1\b[^>]*>\s*[^<\s]', text):
                errors.append(f'{lang}/{rd.stem}: missing topic title')
            for figure in figures:
                if not (path.parent/figure).is_file():
                    errors.append(f'{lang}/{rd.stem}: missing figure {figure}')
            rows.append(dict(topic=rd.stem, figures=figures, text_output=text_output))
        records[lang] = rows
    (out/'reference-output.json').write_text(json.dumps(dict(
        status='failed' if errors else 'passed', topics=records, errors=errors), indent=2)+'\n')
    if errors:
        raise RuntimeError('Reference output audit failed:\n'+'\n'.join(errors))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--publication',action='store_true')
    parser.add_argument('--render-jobs',type=int,choices=(1,2),default=1,
                        help='Maximum concurrent native locale renderers (default: 1)')
    parser.add_argument('--render-batch-size',type=int,choices=range(1,13),default=4,
                        help='Topics/articles per independent R worker (default: 4)')
    parser.add_argument('--no-cache',action='store_true',help='Force fresh native rendering; do not read or stage cached results')
    args=parser.parse_args()
    os.environ['TZ'] = 'UTC'
    build_started = time.monotonic()
    subprocess.run(["python3", "tools/check-documentation.py", "--source-only"], cwd=ROOT, check=True)
    render_snapshot = render_source_snapshot()
    out=ROOT/'.site'
    if out.is_symlink():raise SystemExit('Refusing a symlink output directory')
    if out.exists():shutil.rmtree(out)
    out.mkdir();assets=out/'assets';assets.mkdir()
    cache = Cache(out,enabled=not args.no_cache)
    environment = runtime() if not args.no_cache else {}
    shell_identity = dict(schema=1,runtime=cache_key(environment),inputs=shell_inputs())
    reference_identity = dict(shell=shell_identity,inputs=reference_inputs())
    diagnostics = ROOT/'.cache/website-diagnostics'
    diagnostics.mkdir(parents=True, exist_ok=True)
    (diagnostics/'cache-runtime.json').write_text(json.dumps(environment,sort_keys=True)+'\n')
    (diagnostics/'cache-inputs.json').write_text(json.dumps(dict(shell=shell_identity,reference=reference_identity),sort_keys=True)+'\n')
    cache.write()
    timings = dict(status='in_progress', phases={}, workers={})
    def save_timings():
        timings['elapsed_seconds'] = round(time.monotonic() - build_started, 3)
        (out/'build-timings.json').write_text(json.dumps(timings, indent=2)+'\n')
    def phase(name, started):
        elapsed = round(time.monotonic() - started, 3)
        timings['phases'][name] = elapsed
        save_timings()
        print(f'Completed {name}: {elapsed:.1f}s', flush=True)
    phase('source_preflight', build_started)
    preparation_started = time.monotonic()
    metadata=json.loads((WEB/'version.json').read_text())
    metadata.update(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        build_time=datetime.now(timezone.utc).isoformat(),
        example_runtime=subprocess.check_output(['Rscript','--vanilla','-e',
          'cat(R.version.string, "; ggplot2 ", as.character(packageVersion("ggplot2")), sep="")'],cwd=ROOT,text=True).strip(),
        deployment='prepared for authorized deployment; builder does not deploy' if args.publication else 'local preview; not deployed')
    (out/'build-metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    for name in ('logo.webp', 'favicon.png', 'icon.json', 'site.css', 'site.js'):
        shutil.copy2(WEB/'assets'/name, assets/name)
    for name in ('logo.png', 'logo.svg'):
        shutil.copy2(ROOT/'man/figures'/name, assets/name)
    for row in json.loads((WEB/'home-previews.json').read_text()):
        source = ROOT/row['file']
        shutil.copy2(source, assets/source.name)
    for name in ('website-data-sources.tsv', 'flagship-data-sources.tsv'):
        shutil.copy2(ROOT/'data-raw'/name, assets/name)
    shutil.copy2(ROOT/'NEWS.md', assets/'NEWS.md')
    shutil.copy2(WEB/'ASSET-SOURCES.md', assets/'ASSET-SOURCES.md')
    copy_source_downloads(assets)
    subprocess.run(['Rscript', '--vanilla', 'tools/website/parameter-directory.R',
                    str(ROOT), str(out/'parameter-directory.json')], cwd=ROOT, check=True)
    render_handbook_examples(assets)
    complete_examples(assets)
    input_hashes = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in sorted((WEB/'examples').glob('*.R'))}
    assert_render_source_snapshot(render_snapshot, input_hashes)
    (out/'source-snapshot.json').write_text(json.dumps(dict(
        status='started; completion unverified', render_inputs=render_snapshot,
        archive_inputs=input_hashes),indent=2)+'\n')
    pandoc_dir=subprocess.check_output(['Rscript','-e','cat(rmarkdown::find_pandoc()$dir)'],cwd=ROOT,text=True).strip()
    phase('asset_preparation', preparation_started)
    markdown_started = time.monotonic()
    payload={lang:[] for lang in ['en','zh']}
    for lang in payload:
        for page in PAGES:
            if page['slug']=='reference' or page.get('kind')=='alias':continue
            md=expand((WEB/lang/(page['slug']+'.md')).read_text(),lang,base='../assets')
            md_identity = dict(schema=1,language=lang,slug=page['slug'],markdown=md,pandoc=environment.get('pandoc'),generator='markdown+pipe_tables+header_attributes/html5/wrap-none')
            body_dir=out/'.native-markdown'/lang/page['slug']
            if cache.restore('markdown',md_identity,body_dir):
                body=(body_dir/'body.txt').read_text()
            else:
                body=subprocess.check_output([str(Path(pandoc_dir)/'pandoc'),'--from=markdown+pipe_tables+header_attributes',
                    '--to=html5','--wrap=none'],input=md,text=True)
                body_dir.mkdir(parents=True,exist_ok=True);(body_dir/'body.txt').write_text(body)
                cache.stage('markdown',md_identity,body_dir)
            headings=re.findall(r'<h[23] id="([^"]+)"[^>]*>(.*?)</h[23]>',body,re.S)
            toc=''.join(f'<li><a href="#{a}">{re.sub("<[^>]+>","",t)}</a></li>' for a,t in headings)
            sidebar=('<h2>本页目录</h2>' if lang=='zh' else '<h2>On this page</h2>')+'<ul class="toc-list">'+toc+'</ul>'
            payload[lang].append(dict(slug=page['slug'],title=page['title'][lang],body=body,sidebar=sidebar))
    article_identities={lang:{page['slug']:dict(shell=shell_identity,language=lang,page=page) for page in payload[lang]} for lang in payload}
    all_payload={lang:list(rows) for lang,rows in payload.items()}
    # Initialization refreshes legal/home/template assets before restoring raw pages.
    for lang in payload:
        subprocess.run(['Rscript','--vanilla','tools/render-website.R',lang,'--part=init'],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
        missing=[]
        for page in payload[lang]:
            target=out/'.native-article'/lang/page['slug']
            if cache.restore('article',article_identities[lang][page['slug']],target):
                shutil.copy2(target/'page.txt',out/lang/(page['slug']+'.html'))
            else:missing.append(page)
        payload[lang]=missing
    allow_reference_cache=reference_cache_allowed()
    reference_hits={lang:cache.restore('reference',dict(reference_identity,language=lang),out/lang/'reference') if allow_reference_cache else False for lang in payload}
    (out/'page-content.json').write_text(json.dumps(payload,ensure_ascii=False)+'\n')
    phase('markdown_conversion', markdown_started)
    rendering_started = time.monotonic()
    # Locales write separate destinations and share only immutable source files.
    languages = ['en', 'zh']
    def render_locale(lang):
        # Independent bounded R workers avoid accumulating a whole site's
        # template state in one process. Every actual return code is retained.
        jobs=[['--part=reference','--offset='+str(i),'--count='+str(args.render_batch_size)]
              for i in range(0,len(list((ROOT/'man').glob('*.Rd'))),args.render_batch_size)] if not reference_hits[lang] else []
        jobs += [['--part=articles','--offset='+str(i),'--count='+str(args.render_batch_size)]
              for i in range(0,len(payload[lang]),args.render_batch_size)]
        results=[];logs=[]
        (out/(lang+'-workers.json')).write_text('[]\n')
        (out/(lang+'-build.log')).write_text('')
        for index, options in enumerate(jobs):
            snapshot_started = time.monotonic()
            assert_render_source_snapshot(render_snapshot,input_hashes)
            snapshot_seconds = time.monotonic() - snapshot_started
            worker_started = time.monotonic()
            print(f'Starting {lang} worker {index + 1}/{len(jobs)}: '+ ' '.join(options), flush=True)
            attempts = []
            for attempt in (1, 2):
                suffix = '' if attempt == 1 else '-retry'
                logfile = out/f'{lang}-worker-{index}{suffix}.log'
                with logfile.open('w') as log:
                    result = subprocess.run(['Rscript','--vanilla','tools/render-website.R',lang,*options],
                                            cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
                logs.append(logfile.read_text())
                attempts.append(dict(return_code=result.returncode, log=logfile.name))
                if not result.returncode:
                    break
                if attempt == 1:
                    assert_render_source_snapshot(render_snapshot, input_hashes)
                    print(f'Retrying {lang} worker {index+1} after exit {result.returncode}', flush=True)
            results.append(dict(options=options,return_code=result.returncode,log=logfile.name,attempts=attempts,
                                elapsed_seconds=round(time.monotonic()-worker_started,3),
                                snapshot_check_seconds=round(snapshot_seconds,3)))
            print(f'Completed {lang} worker {index + 1}/{len(jobs)}: '
                  f'{results[-1]["elapsed_seconds"]:.1f}s, exit {result.returncode}', flush=True)
            (out/(lang+'-workers.json')).write_text(json.dumps(results,indent=2)+'\n')
            (out/(lang+'-build.log')).write_text('\n'.join(logs))
            if result.returncode:break
        return lang,0 if len(results)==len(jobs) and all(r['return_code']==0 for r in results) else 1
    with ThreadPoolExecutor(max_workers=args.render_jobs) as pool:
        statuses=dict(pool.map(render_locale,languages))
    (out/'render-processes.json').write_text(json.dumps(dict(
        render_jobs=args.render_jobs, return_codes=statuses),indent=2)+'\n')
    phase('native_rendering', rendering_started)
    timings['workers'] = {lang: json.loads((out/(lang+'-workers.json')).read_text()) for lang in languages}
    save_timings()
    combined='\n'.join((out/(lang+'-build.log')).read_text() for lang in languages)
    (out/'build.log').write_text(combined)
    if any(statuses.values()):
        timings['status'] = 'renderers_failed'
        save_timings()
        raise SystemExit('Native locale renderers failed: '+json.dumps(statuses)+'\n'+combined[-7000:])
    # Stage raw output before finalization mutates menus, anchors and page shells.
    for lang in languages:
        if allow_reference_cache and not reference_hits[lang]:cache.stage('reference',dict(reference_identity,language=lang),out/lang/'reference')
        for page in payload[lang]:
            folder=out/'.native-article'/lang/page['slug'];folder.mkdir(parents=True,exist_ok=True)
            shutil.copy2(out/lang/(page['slug']+'.html'),folder/'page.txt')
            cache.stage('article',article_identities[lang][page['slug']],folder)
    shutil.rmtree(out/'.native-article',ignore_errors=True)
    shutil.rmtree(out/'.native-markdown',ignore_errors=True)
    finishing_started = time.monotonic()
    try:
        audit_reference_output(out)
        assert_render_source_snapshot(render_snapshot, input_hashes)
        finalize(out,args,render_snapshot,input_hashes)
    except Exception:
        timings['status'] = 'audit_or_finalization_failed'
        save_timings()
        raise
    phase('output_audit_and_finalization', finishing_started)
    timings['status'] = 'complete'
    save_timings()

def page_redirect_script(mapping, default_for_unknown):
    fallback = 'm[""]' if default_for_unknown else 'null'
    return ('<script>(function(){const m='+json.dumps(mapping, ensure_ascii=False)+';'
            'function go(){let f=location.hash.slice(1);try{f=decodeURIComponent(f)}catch(e){}'
            'const d=Object.prototype.hasOwnProperty.call(m,f)?m[f]:'+fallback+';'
            'if(d && d!==location.pathname.split("/").pop()+location.hash)location.replace(d)}'
            'go();addEventListener("hashchange",go)})()</script>')


def write_reference_redirect(out, path, destination, languages):
    labels = {'en': 'Version selection', 'zh': '版本选择'}
    links = []
    for lang in languages:
        relative = os.path.relpath(out/lang/destination, path.parent)
        tag = 'zh-Hans' if lang == 'zh' else 'en'
        links.append('<a lang="' + tag + '" hreflang="' + tag + '" href="' +
                     html.escape(relative, quote=True) + '">' + labels[lang] + '</a>')
    target = os.path.relpath(out/languages[0]/destination, path.parent)
    script = '<script>location.replace(' + json.dumps(target) + ')</script>'
    tag = 'zh-Hans' if languages == ('zh',) else 'en'
    path.write_text('<!doctype html><html lang="' + tag + '"><head><meta charset="utf-8">'
                    '<title>' + labels[languages[0]] + ' · ggchord</title>' + script +
                    '</head><body><main>' + ' · '.join(links) + '</main></body></html>')


def finalize(out,args,render_snapshot,input_hashes):
    inputs=[ROOT/name for name in input_hashes]
    assets=out/'assets'
    metadata=json.loads((out/'build-metadata.json').read_text())
    # Post-process using the same XML library used by pkgdown; no hand-built page shell.
    subprocess.run(['Rscript','tools/finalize-website.R', 'publication' if args.publication else 'preview'],cwd=ROOT,check=True)
    migrations=json.loads((WEB/'page-migrations.json').read_text())
    shutil.copy2(WEB/'page-migrations.json',assets/'page-migrations.json')
    for lang in ['en','zh']:
        for page in PAGES:
            name=page['slug']+'.html'
            mapping=migrations.get(name,{'':page['canonical_path']})
            alias = page.get('kind')=='alias' or page['slug']=='reference'
            script=page_redirect_script(mapping, default_for_unknown=alias)
            path=out/lang/name
            if alias:
                links=''.join('<li><a href="'+html.escape(dest,quote=True)+'">'+html.escape(anchor or page['title'][lang])+'</a></li>' for anchor,dest in mapping.items())
                path.write_text('<!doctype html><html lang="'+('zh-Hans' if lang=='zh' else 'en')+'"><meta charset="utf-8"><title>Moved · ggchord</title><h1>'+('章节已迁移' if lang=='zh' else 'Section moved')+'</h1><p>'+('旧章节的具体新位置如下。' if lang=='zh' else 'Find the exact destination of each old section below.')+'</p><ul>'+links+'</ul>'+script+'</html>')
            elif name in migrations:
                path.write_text(path.read_text().replace('</head>',script+'</head>',1))
    # Keep project presentation files identical to the source digested by finalization.
    for name in ['site.css','site.js']: shutil.copy2(WEB/'assets'/name,assets/name)
    extra_migrations=set(migrations)-{p['slug']+'.html' for p in PAGES}
    retired_references = {name for name in extra_migrations
                          if name.startswith('reference/')}
    for lang in ['en','zh']:
        for name in extra_migrations:
            mapping=migrations[name]
            if not mapping.get(''):continue
            path = out/lang/name
            path.parent.mkdir(parents=True, exist_ok=True)
            targets = {anchor: os.path.relpath(out/lang/dest, path.parent)
                       for anchor, dest in mapping.items()}
            if name in retired_references:
                write_reference_redirect(out, path, mapping[''], (lang,))
                continue
            script=page_redirect_script(targets, default_for_unknown=True)
            path.write_text('<!doctype html><html lang="'+('zh-Hans' if lang=='zh' else 'en')+'"><meta charset="utf-8"><title>ggchord</title>'+script+'<a href="'+html.escape(targets[''], quote=True)+'">Continue / 继续</a></html>')
    for lang in ['en','zh']:
        index=[]
        for path in sorted((out/lang).rglob('*.html')):
            if str(path.relative_to(out/lang)) in extra_migrations or path.name=='reference.html' or any(p['slug']+'.html'==str(path.relative_to(out/lang)) and p.get('kind')=='alias' for p in PAGES):continue
            content=path.read_text()
            title=html.unescape(re.search(r'<title>(.*?)</title>',content,re.S)[1]).split(' • ')[0]
            main=re.search(r'<main\b[^>]*>(.*?)</main>',content,re.S)
            text=' '.join(PageText(main[1] if main else content).parts)
            index.append(dict(url=str(path.relative_to(out/lang)),title=title,text=re.sub(r'\s+',' ',text)))
        (out/lang/'search.json').write_text(json.dumps(search_records(out,lang,index),ensure_ascii=False)+'\n')
    # Stable historical paths remain usable after a future, separately authorized cutover.
    topics=[p.stem for p in (ROOT/'man').glob('*.Rd')]
    routes={'index.html':'en/index.html','articles/ggchord_vignette.html':'en/formats.html',
        'articles/ggchord_guide_hans.html':'zh/formats.html','reference/index.html':'en/reference/index.html','news/index.html':'assets/NEWS.md'}
    routes.update({'reference/'+s+'.html':'en/reference/'+s+'.html' for s in topics})
    # Keep the retired public names out of canonical search and reference lists.
    # Their localized redirects lead to version selection.
    routes.update({name: 'en/'+name for name in extra_migrations if name.startswith('reference/')})
    aliases={a:p.stem for p in (ROOT/'man').glob('*.Rd')
             for a in re.findall(r'\\alias\{([^}]+)\}',p.read_text())}
    for rel in json.loads((WEB/'historical-pages.json').read_text()):
        if rel.startswith('reference/') and Path(rel).stem in aliases:
            routes[rel]='en/reference/'+aliases[Path(rel).stem]+'.html#ref-usage'
    routes.update({name:'en/'+name for name in ['authors.html','LICENSE.html','LICENSE-text.html']})
    routes.update({'README-Hans.html':'zh/index.html','DESIGN-ROADMAP.html':'en/about.html',
                   'articles/index.html':'en/modes.html'})
    for old,target in routes.items():
        path=out/old;path.parent.mkdir(parents=True,exist_ok=True);dest=os.path.relpath(out/target,path.parent)
        if old in retired_references:
            write_reference_redirect(out, path, migrations[old][''], ('en', 'zh'))
            continue
        extra='<script>try{if(localStorage.getItem("ggchord-language")==="zh")location.replace("zh/index.html");else location.replace("en/index.html")}catch(e){location.replace("en/index.html")}</script>' if old=='index.html' else '<script>location.replace('+json.dumps(dest)+'+('+json.dumps('' if '#' in dest else '#')+'===""?"":location.hash))</script>'+f'<noscript><meta http-equiv="refresh" content="0;url={dest}"></noscript>'
        path.write_text(f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>ggchord</title>{extra}</head><body><a href="{dest}">English</a> · <a href="{os.path.relpath(out/"zh/index.html",path.parent)}">现代汉语</a></body></html>')
    (out/'404.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>Page not found / 未找到页面</title><h1>Page not found / 未找到页面</h1><p><a href="en/index.html">English home</a> · <a href="zh/index.html">现代汉语首页</a></p></html>')
    routes['404.html']='404.html'
    retired=[]
    for rel in json.loads((WEB/'historical-pages.json').read_text()):
        if rel in routes:continue
        path=out/rel;path.parent.mkdir(parents=True,exist_ok=True)
        if rel.startswith('reference/'):
            write_reference_redirect(out, path, 'migration.html#api', ('en', 'zh'))
            retired.append(rel)
            continue
        en=os.path.relpath(out/'en/reference/index.html',path.parent);zh=os.path.relpath(out/'zh/reference/index.html',path.parent)
        path.write_text(f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Historical documentation</title></head><body><h1>Historical documentation / 历史文档</h1><p>This historical URL has no exact current equivalent. 此旧地址没有一一对应的当前页面。</p><a href="{en}">Current reference</a> · <a href="{zh}">现代汉语参考</a></body></html>')
        retired.append(rel)
    for name,data in [('legacy-routes',routes),('retired-routes',retired)]:
        (out/(name+'.json')).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    files = list(public_files(p for base in [WEB, ROOT/'R', ROOT/'man', ROOT/'vignettes']
                              for p in base.rglob('*')))
    files += list(public_files((ROOT/'tools').glob('*')))
    files += list(public_files((ROOT/'tools/site-audit').glob('*.py')))
    files += list(public_files((ROOT/'tools/site-audit').glob('*.R')))
    files += [ROOT/name for name in render_snapshot]
    files += inputs + sorted((ROOT/'data').glob('*.rda'))
    files += [ROOT/p for p in ['DESCRIPTION','NAMESPACE','NEWS.md','README.md','README-Hans.md','MIGRATION.md','MIGRATION-Hans.md','_pkgdown.yml','.github/workflows/documentation.yaml']]
    manifest=dict(version=re.search(r'^Version: (.+)$',(ROOT/'DESCRIPTION').read_text(),re.M)[1],api_inventory='current NAMESPACE and man topics; reference-output.json records evaluated examples',generator='pkgdown + paired Markdown + shared static assets',
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)),
        deployment='prepared for authorized deployment; builder does not deploy' if args.publication else 'local preview; not deployed',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    manifest.update(metadata)
    assert_render_source_snapshot(render_snapshot, input_hashes)
    (out/'source-snapshot.json').write_text(json.dumps(dict(
        status='verified unchanged across rendering and finalization',
        render_inputs=render_snapshot,archive_inputs=input_hashes),indent=2)+'\n')
    (out/'build-manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    (out/'page-content.json').unlink()
    shutil.rmtree(out/'reference-zh')
    print(f'Built {len(PAGES)} paired entries with pkgdown and {len(topics)} reference topics; reference examples evaluated; no deployment or vignette execution.')

if __name__=='__main__':main()
