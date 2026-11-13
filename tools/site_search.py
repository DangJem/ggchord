"""Structured local search index and complete example downloads from public lessons."""
import html,json,re,hashlib
from pathlib import Path
from urllib.parse import quote
from documentation_lib import ROOT,WEB,PAGES,RECIPES,snippets,expand
from site_registry import page_kind

def complete_examples(assets):
    folder=assets/'examples';folder.mkdir(exist_ok=True)
    chunks=snippets();records=[]
    for lang in ['en','zh']:
        for page in PAGES:
            slug=page['slug']
            if slug=='installation' or page.get('kind')=='alias':continue
            text=(WEB/lang/(slug+'.md')).read_text()
            # Inline the exact public preparation source, including recipe functions.
            for name,keys in RECIPES.items():text=text.replace('{{recipe:'+name+'}}','{{prepare:'+','.join(keys)+'}}')
            if not re.search(r'\{\{(?:code|prepare):|```r\n',text):continue
            expanded=expand(text,lang,figures='../assets')
            blocks=re.findall(r'```r\n(.*?)```',expanded,re.S)
            if not blocks:continue
            code='library(ggchord)\n\n'+'\n\n'.join(blocks)+'\n'
            name=slug+'-'+lang+'.R';(folder/name).write_text(code)
            records.append(dict(page=slug,language=lang,file='assets/examples/'+name,sha256=hashlib.sha256(code.encode()).hexdigest(),scope='All page R blocks and inlined public preparations; file examples require the input archive; preview/save calls retained'))
    (assets/'examples.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')

def search_records(out,lang,page_records):
    records=[]
    for page in page_records:
        kind=page_kind(page['url'])
        page['kind']='faq' if page['url']=='faq.html' else 'function' if kind=='api' else 'course' if kind in ['tutorial','module','case','index'] else 'tool'
        page['keywords']=page['url'].removesuffix('.html')
        page['summary']=page['text'][:240]
        records.append(page)
        # Hit-level context uses only public rendered prose; executable blocks are excluded.
        path=out/lang/page['url']
        raw=path.read_text()
        main=re.search(r'<main\b[^>]*>(.*?)</main>',raw,re.S)
        body=main[1] if main else ''
        headings=list(re.finditer(r'<h([23])\b[^>]*id="([^"]+)"[^>]*>(.*?)</h\1>',body,re.S))
        if page['url']!='parameter-atlas.html':
            for i,h in enumerate(headings):
                text=body[h.end():headings[i+1].start() if i+1<len(headings) else len(body)]
                text=re.sub(r'<pre\b[^>]*>.*?</pre>','',text,flags=re.S)
                text=html.unescape(re.sub(r'<[^>]+>',' ',text))
                text=re.sub(r'\s+',' ',text).strip()
                section=html.unescape(re.sub(r'<[^>]+>','',h[3]))
                records.append(dict(kind=page['kind'],url=page['url']+'#'+h[2],title=page['title']+' · '+section,
                    section=section,keywords=section,text=text[:6000],summary=text[:300]))

    contracts=json.loads((ROOT/'tools/teaching-audit/parameter-contracts.json').read_text())
    titles={p['slug']:p['title'][lang] for p in PAGES}
    decisions=json.loads((WEB/'teaching-reviews.json').read_text())
    for r in contracts:
        api=r['api'];parameter=r['parameter'];review=decisions[api]
        records.append(dict(kind='parameter',url='parameter-atlas.html#'+quote('param-'+api+'-'+parameter,safe='._-'),title=api+' · '+parameter,
            section=titles['parameter-atlas'],keywords=parameter+' '+api+' '+review['purpose'][lang],text=r['constraint'],summary=(('默认：' if lang=='zh' else 'Default: ')+(r['default'] or ('必需参数' if lang=='zh' else 'required'))+' · '+review['purpose'][lang])))
    for r in json.loads((WEB/'teaching-tasks.json').read_text()):
        slug,_,anchor=r['page'].partition('#')
        records.append(dict(kind='task',url=slug+'.html'+('#'+anchor if anchor else ''),title=r['task'][lang],section=titles[slug],keywords=r['api']+' '+r['parameters'],text=r['task'][lang]+' '+r['parameters'],summary=r['api']+' · '+r['parameters']+' · '+r['modes']))
    aliases={a:p.stem for p in (ROOT/'man').glob('*.Rd') for a in re.findall(r'\\alias\{([^}]+)\}',p.read_text())}
    for api,r in decisions.items():
        slug,_,anchor=(r['page']+'#'+r['anchor']).partition('#')
        records.append(dict(kind='function',url='reference/'+aliases[api]+'.html#ref-usage',title=api,section=titles[slug],keywords=api+' '+r['purpose'][lang],text=r['purpose'][lang],summary=r['purpose'][lang]))
    for intent in json.loads((WEB/'search-intents.json').read_text()):
        records.append(dict(kind='task',url=intent['url'],title=intent['title'][lang],queries=intent['queries'],keywords=' '.join(intent['queries']),text=intent['title'][lang],summary=intent['summary'][lang]))
    return records
