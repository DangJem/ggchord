#!/usr/bin/env python3
"""Apply reviewed canonical destinations without guessing page parents or API meaning."""
import json,re
from pathlib import Path
from documentation_lib import ROOT,WEB,snippets
from site_registry import canonical

def save(path,data):path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')

def main():
    pages=json.loads((WEB/'pages.json').read_text())
    for p in pages:
        if p['kind']=='alias':continue
        for lang in ['en','zh']:
            path=WEB/lang/(p['slug']+'.md');text=path.read_text()
            def dest(url):
                target=canonical(url,p['slug'])
                return target if target!=p['slug']+'.html'+url else url
            text=re.sub(r'\]\(([^\s)]+)\)',lambda m:']('+dest(m[1])+')',text)
            text=re.sub(r'href="([^"]+)"',lambda m:'href="'+dest(m[1])+'"',text)
            path.write_text(text)
        text=(WEB/'en'/(p['slug']+'.md')).read_text()
        keys=[k for v in re.findall(r'\{\{(?:code|prepare):([\w,]+)\}\}',text) for k in v.split(',')]
        recipes=json.loads((WEB/'showcase-recipes.json').read_text())
        for key in re.findall(r'\{\{recipe:(\w+)\}\}',text):keys.extend(recipes[key])
        p['recipe_ids']=list(dict.fromkeys(keys))
        code='\n'.join(snippets()[k] for k in keys)
        p['apis']=list(dict.fromkeys(re.findall(r'\b([A-Za-z][\w.]*)\s*\(',code)))
    save(WEB/'pages.json',pages)
    for filename in ['teaching-tasks.json','teaching-reviews.json','figure-lessons.json','parameter-comparisons.json','coverage.json']:
        path=WEB/filename;data=json.loads(path.read_text())
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

if __name__=='__main__':main()
