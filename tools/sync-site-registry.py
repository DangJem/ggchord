#!/usr/bin/env python3
"""Apply reviewed canonical destinations without guessing page parents or API meaning."""
import json,re
from pathlib import Path
from documentation_lib import ROOT,WEB,snippets
from site_registry import canonical

def save(path,data):
    canonical=WEB/'content-zh'/path.relative_to(WEB)
    if canonical.exists():
        import hashlib
        translations=WEB/'translations/en/content'/path.relative_to(WEB)
        saved=json.loads(translations.read_text())
        # Registry updates may change shared destinations/metadata, never prose.
        def strings(value,pointer=()):
            if isinstance(value,dict):
                if set(value)=={'zh'}:return {'/'.join(map(str,pointer)):value['zh']}
                return {k:v for key,item in value.items() for k,v in strings(item,(*pointer,key)).items()}
            if isinstance(value,list):return {k:v for i,item in enumerate(value) for k,v in strings(item,(*pointer,i)).items()}
            return {}
        assert strings(json.loads(canonical.read_text()))==strings(data),'Registry must not change localized prose'
        canonical.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
        payload=json.dumps(data,ensure_ascii=False,sort_keys=True,separators=(',',':'))
        saved['source_sha256']=hashlib.sha256(payload.encode()).hexdigest()
        translations.write_text(json.dumps(saved,ensure_ascii=False,indent=2)+'\n')
    else:path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')

def main():
    pages=json.loads((WEB/'content-zh/pages.json').read_text())
    for p in pages:
        if p['kind']=='alias':continue
        for lang in ['zh']:
            path=WEB/lang/(p['slug']+'.md');text=path.read_text()
            def dest(url):
                target=canonical(url,p['slug'])
                return target if target!=p['slug']+'.html'+url else url
            text=re.sub(r'\]\(([^\s)]+)\)',lambda m:']('+dest(m[1])+')',text)
            text=re.sub(r'href="([^"]+)"',lambda m:'href="'+dest(m[1])+'"',text)
            path.write_text(text)
        text=(WEB/'zh'/(p['slug']+'.md')).read_text()
        keys=[k for v in re.findall(r'\{\{(?:code|prepare):([\w,]+)\}\}',text) for k in v.split(',')]
        recipes=json.loads((WEB/'showcase-recipes.json').read_text())
        for key in re.findall(r'\{\{recipe:(\w+)\}\}',text):keys.extend(recipes[key])
        p['recipe_ids']=list(dict.fromkeys(keys))
        code='\n'.join(snippets()[k] for k in keys)
        p['apis']=list(dict.fromkeys(re.findall(r'\b([A-Za-z][\w.]*)\s*\(',code)))
    save(WEB/'pages.json',pages)
    for filename in ['teaching-tasks.json','teaching-reviews.json','figure-lessons.json','parameter-comparisons.json','coverage.json']:
        path=WEB/filename;canonical_path=WEB/'content-zh'/filename
        data=json.loads((canonical_path if canonical_path.exists() else path).read_text())
        if filename=='teaching-reviews.json':
            for r in data.values():
                dest=canonical(r['page']+'.html#'+r['anchor']);r['page'],_,r['anchor']=dest.replace('.html','#',1).partition('#');r['anchor']=r['anchor'].lstrip('#')
                r['case']=canonical(r['case']+'.html').split('.html')[0]
        elif filename=='teaching-tasks.json':
            for r in data:r['page']=canonical(r['page'].replace('#','.html#',1) if '#' in r['page'] else r['page']+'.html').replace('.html','',1)
        elif filename=='figure-lessons.json':
            data={key:canonical(value.replace('#','.html#',1)).replace('.html','',1) for key,value in data.items()}
        elif filename=='parameter-comparisons.json':
            for r in data:r['page']=canonical(r['page']+'.html').split('.html')[0]
        elif filename=='coverage.json':
            for r in data:
                keys=set(r['snippets']);figures=set(r['figures'])
                # This locates actual shared material after splitting pages; it
                # does not infer or upgrade the editorial teaching depth.
                r['articles']=[p['slug'] for p in pages if p['kind']!='alias' and
                    (keys.intersection(p['recipe_ids']) or figures.intersection(p['figure_ids']))]
        save(path,data)
    import subprocess
    subprocess.run(["python3","tools/sync-website-languages.py"],cwd=ROOT,check=True)

if __name__=='__main__':main()
