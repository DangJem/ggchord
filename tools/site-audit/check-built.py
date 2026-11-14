#!/usr/bin/env python3
"""Check audited URLs, declared section destinations and actual download assets."""
import hashlib,html,json,re,subprocess,sys,zipfile
from pathlib import Path
from urllib.parse import unquote,urlsplit
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
from site_registry import canonical
from documentation_lib import snippets

def main():
    site=ROOT/'.site';baseline=json.loads((ROOT/'tools/site-audit/baseline.json').read_text())
    failures=[]
    def require(ok,message):
        if not ok:failures.append(message)
    def ids(path):return set(re.findall(r'\bid="([^"]+)"',path.read_text())) if path.exists() else set()
    for row in baseline['urls']:require((site/row['path']).is_file(),'Missing audited URL: '+row['path'])
    sections=[]
    for row in baseline['pages']:
        slug=row['old']
        for lang in ['en','zh']:
            text=subprocess.check_output(['git','show',baseline['source_commit']+':website/'+lang+'/'+slug+'.md'],cwd=ROOT,text=True)
            anchors={a or b for a,b in re.findall(r'\{#([^}]+)\}|\bid="([^"]+)"',text)}
            for anchor in sorted(anchors):
                target=canonical(slug+'.html#'+anchor);url=urlsplit(target);path=site/lang/url.path
                valid=path.is_file() and (not url.fragment or unquote(url.fragment) in ids(path))
                require(valid,lang+'/'+slug+'#'+anchor+' -> missing '+target)
                sections.append(dict(language=lang,old=slug+'.html#'+anchor,target=target,valid=valid))
    old_nav=json.loads(subprocess.check_output(['git','show',baseline['source_commit']+':website/navigation.json'],cwd=ROOT,text=True))
    for lang in ['en','zh']:
        for group in old_nav:
            target=canonical('about.html#path-'+group['id']);url=urlsplit(target)
            require((site/lang/url.path).is_file() and (not url.fragment or url.fragment in ids(site/lang/url.path)),lang+': generated old chapter anchor '+group['id'])
    aliases={a:p.stem for p in (ROOT/'man').glob('*.Rd') for a in re.findall(r'\\alias\{([^}]+)\}',p.read_text())}
    legacy=json.loads((site/'legacy-routes.json').read_text());restored=[]
    for row in baseline['urls']:
        if 'public alias' not in (row.get('action') or ''):continue
        name=Path(row['path']).stem;expected='en/reference/'+aliases[name]+'.html#ref-usage'
        require(legacy.get(row['path'])==expected,'Incorrect public alias destination: '+name);restored.append(name)
    require(len(restored)==14,'Expected 14 restored public aliases')
    captions=json.loads((ROOT/'website/assets/captions.json').read_text())
    for lang in ['en','zh']:
        gallery_ids=ids(site/lang/'gallery.html')
        require(set(captions)<=gallery_ids,lang+': missing original figure ID')
    manifest=json.loads((site/'assets/input-hashes.json').read_text())
    with zipfile.ZipFile(site/'assets/tutorial-inputs.zip') as archive:
        for name,digest in manifest.items():
            require(hashlib.sha256(archive.read(name)).hexdigest()==digest,'Archive bytes: '+name)
            require(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'Archive source: '+name)
    for record in json.loads((site/'assets/examples.json').read_text()):
        path=site/record['file'];require(hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256'],'Example download: '+record['file'])
    courses=json.loads((ROOT/'website/course-figures.json').read_text())
    rendered=json.loads((ROOT/'tools/site-audit/course-rendering.json').read_text())
    chunks=snippets()
    for row in courses:
        code='\n\n'.join(chunks[key] for key in row['snippets'])
        require(hashlib.sha256(code.encode()).hexdigest()==rendered['figures'][row['id']]['code_sha256'],'Course code changed: '+row['id'])
        digest=rendered['figures'][row['id']]['sha256']
        require(hashlib.sha256((site/'assets'/Path(row['file']).name).read_bytes()).hexdigest()==digest,'Course figure: '+row['id'])
    # Verify visible section links against one canonical order, including exits.
    pages={p['slug']:p for p in json.loads((ROOT/'website/pages.json').read_text())}
    navigation=json.loads((ROOT/'website/navigation.json').read_text())
    navigation_checks=0
    for lang in ['en','zh']:
        for group in navigation:
            keys=group['pages']
            for i,key in enumerate(keys):
                page=pages[key]
                if len(keys)<2 or page['canonical_path'].startswith('reference/'):continue
                raw=(site/lang/page['canonical_path']).read_text()
                fragment=re.search(r'<nav class="lesson-pagination"[^>]*>(.*?)</nav>',raw,re.S)
                require(bool(fragment),lang+'/'+key+': missing section navigation')
                if not fragment:continue
                links=[html.unescape(link) for link in re.findall(r'href="([^"]+)"',fragment[1])]
                expected=([pages[keys[i-1]]['canonical_path']] if i else [])+[pages[keys[i+1]]['canonical_path'] if i+1<len(keys) else 'modes.html']
                require(links==expected and page['canonical_path'] not in links,lang+'/'+key+': wrong/self section link')
                require(page['next_steps']==keys[i+1:i+2],key+': navigation metadata differs')
                navigation_checks+=1
        for slug in ['multi-sequence','circular','dotplot','fusion','parameters']:
            raw=(site/lang/(slug+'.html')).read_text()
            crumb=re.search(r'<nav class="breadcrumb".*?</nav>',raw,re.S)
            require(bool(crumb) and not any(bad in crumb[0] for bad in ['installation.html','styling.html']),lang+'/'+slug+': misleading category breadcrumb')
        gallery=(site/lang/'gallery.html').read_text()
        require(not re.search(r'<h3[^>]*>\s*<span class="section-number"',gallery),lang+': gallery titles carry section numbers')
    # Execute the exact current ranking expressions over final indexes. This is
    # an algorithm check, not browser interaction or a simulated DOM.
    search_code=r"""
const fs=require('fs');
const source=fs.readFileSync('website/assets/site.js','utf8');
const tokens=source.split('\n').find(line=>line.includes('const tokens=q=>'));
const start=source.indexOf('        const ranked=records.filter');
const end=source.indexOf('        status.textContent=grouped',start);
if(!tokens||start<0||end<0)throw Error('Ranking source unavailable');
const rank=new Function('records','query','filter',tokens+'\nconst terms=tokens(query);\n'+source.slice(start,end)+'\nreturn grouped.map(({r,score})=>({url:r.url,title:r.title,score}));');
const queries=JSON.parse(fs.readFileSync('website/search-intents.json')).find(row=>row.url==='parameters.html').queries;
const checks=[];
for(const language of ['en','zh']){
 const records=JSON.parse(fs.readFileSync('.site/'+language+'/search.json'));
 for(const query of queries){
  const ranked=rank(records,query.toLocaleLowerCase(),'all');
  if(!ranked.slice(0,3).some(row=>row.url.split('#')[0]==='parameters.html'))throw Error(language+'/'+query+': missing preferred section');
  const keys=ranked.map(row=>row.url.startsWith('parameter-atlas.html#param-')?row.url:row.url.split('#')[0]);
  if(new Set(keys).size!==keys.length)throw Error('Duplicate search results');
  checks.push({language,query,count_after_deduplication:ranked.length,top_three:ranked.slice(0,3)});
 }
}
process.stdout.write(JSON.stringify(checks));
"""
    search_checks=json.loads(subprocess.check_output(['node','-e',search_code],cwd=ROOT,text=True))
    result=dict(audited_urls=len(baseline['urls']),declared_old_sections=len(sections),restored_public_aliases=restored,
        original_figure_ids=len(captions),new_course_pngs=len(courses),archive_inputs=len(manifest),
        canonical_navigation_records=navigation_checks,parameter_synonym_queries=search_checks,
        scope='Local files, declared source/known generated chapter anchors and download bytes; not remote deployment or arbitrary unrecorded historical fragments',failures=failures)
    (ROOT/'tools/site-audit/built-check.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    (ROOT/'tools/site-audit/section-migrations.json').write_text(json.dumps(sections,ensure_ascii=False,indent=2)+'\n')
    if failures:raise SystemExit('\n'.join(failures))
    print('PASS:',result['audited_urls'],'audited URLs,',len(sections),'declared old sections, 14 public aliases, figures/downloads,',navigation_checks,'navigation records and',len(search_checks),'synonym queries')
if __name__=='__main__':main()
