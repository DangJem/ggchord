#!/usr/bin/env python3
"""Check audited URLs, declared section destinations and actual download assets."""
import hashlib,json,re,subprocess,sys,zipfile
from pathlib import Path
from urllib.parse import unquote,urlsplit
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
from site_registry import canonical

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
    for row in courses:
        digest=rendered['figures'][row['id']]['sha256']
        require(hashlib.sha256((site/'assets'/Path(row['file']).name).read_bytes()).hexdigest()==digest,'Course figure: '+row['id'])
    result=dict(audited_urls=len(baseline['urls']),declared_old_sections=len(sections),restored_public_aliases=restored,
        original_figure_ids=len(captions),new_course_pngs=len(courses),archive_inputs=len(manifest),
        scope='Local files, declared source/known generated chapter anchors and download bytes; not remote deployment or arbitrary unrecorded historical fragments',failures=failures)
    (ROOT/'tools/site-audit/built-check.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    (ROOT/'tools/site-audit/section-migrations.json').write_text(json.dumps(sections,ensure_ascii=False,indent=2)+'\n')
    if failures:raise SystemExit('\n'.join(failures))
    print('PASS:',result['audited_urls'],'audited URLs,',len(sections),'declared old sections, 14 public aliases, figures and downloads')
if __name__=='__main__':main()
