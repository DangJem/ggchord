#!/usr/bin/env python3
"""Current export/formal and expanded-page evidence. No automatic depth certification."""
import csv, json, re, subprocess, hashlib
from collections import Counter
from documentation_lib import ROOT, WEB, PAGES, snippets, expand

def main():
    revision = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    version = re.search(r'^Version: (.+)$',(ROOT/'DESCRIPTION').read_text(),re.M)[1]
    raw = subprocess.check_output(['Rscript','tools/teaching-inventory.R'],cwd=ROOT,text=True)
    entries = json.loads(raw)
    chunks = snippets()
    pages = {p['slug']:(WEB/'en'/(p['slug']+'.md')).read_text() for p in PAGES if p['slug'] not in ('task-index',)}
    rd = {}
    for path in (ROOT/'man').glob('*.Rd'):
        for alias in re.findall(r'\\alias\{([^}]+)\}',path.read_text()): rd[alias] = str(path.relative_to(ROOT))
    reviews=json.loads((WEB/'teaching-reviews.json').read_text())
    source_aliases={}
    for path in (ROOT/'R').glob('*.R'):
        source_aliases.update(re.findall(r'^(\w+)\s*<-\s*(\w+)\s*$',path.read_text(),re.M))
    expanded_pages={slug:expand(re.sub(r'\{\{figure:\w+\}\}', '', text),'en') for slug,text in pages.items()}
    rows=[]; params=[]
    for e in entries:
        name=e['name']; pattern=r'(?<![\w.])'+re.escape(name)+r'\s*\('
        evidence=[]; mentions=[]
        for slug, text in pages.items():
            expanded=expanded_pages[slug]
            if name in text: mentions.append(slug)
            code='\n'.join(re.findall(r'```r\n(.*?)```',expanded,re.S))
            # R AST inventory supplies signatures; call sites here are candidates for review.
            if re.search(pattern,code): evidence.append(slug)
        chunk_ids=[k for k,v in chunks.items() if re.search(pattern,v)]
        # Deliberately conservative: occurrence alone never proves L2 or execution.
        level='L1' if evidence or mentions else 'L0'
        row=dict(api=name,owner=e['owner'],alias=source_aliases.get(name,''),same_object_as=e['alias'],ownership='ggchord' if e['owner']=='ggchord' else 'upstream/object',family=name.split('_')[0],rd=rd.get(name,''),level=level,
            pages=';'.join(evidence),mentions=';'.join(mentions),snippets=';'.join(chunk_ids),
            validation='static; teaching depth and mode require review',
            gap='Needs semantic review of inputs, effects, constraints and independent execution',
            version=version,source_sha=revision)
        review=reviews.get(name)
        row['modes']=review['modes'] if review else 'not individually reviewed; consult Rd'
        row['review_reason']=review['reason'] if review else ''
        if review:
            row['level']=review['level']
            row['validation']='semantic review; runtime logs separate'
        rows.append(row)
        for parameter,default in (e['parameters'] or {}).items():
            params.append(dict(api=name,parameter=parameter,default=default,priority='standard ggplot2 interface' if parameter in ('mapping','data','position','show.legend','inherit.aes','...') else 'review by task',rd=row['rd'],pages=row['pages'],
                level=review['level'] if review and parameter in review['parameters'] else ('L1' if any(re.search(r'(?<![\w.])'+re.escape(parameter)+r'(?![\w.])', expanded_pages[p]) for p in evidence+mentions) else 'L0'),
                validation='formal inspected; parameter teaching not certified',version=version,source_sha=revision))
    out=ROOT/'tools/teaching-audit';out.mkdir(exist_ok=True)
    for filename, data in [('teaching-api-v0.16.7.csv',rows),('teaching-parameters-v0.16.7.csv',params)]:
        with (out/filename).open('w') as f:
            writer=csv.DictWriter(f,fieldnames=data[0],lineterminator='\n');writer.writeheader();writer.writerows(data)
    sources=[ROOT/'NAMESPACE',ROOT/'DESCRIPTION',ROOT/'tools/documentation-examples.R',*sorted((ROOT/'R').glob('*.R')),*sorted((ROOT/'man').glob('*.Rd')),*sorted((WEB/'en').glob('*.md')),*sorted((WEB/'zh').glob('*.md'))]
    summary={'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'source_sha':revision,'version':version,'exports':len(rows),'parameters':len(params),'owned':sum(r['ownership']=='ggchord' for r in rows),'upstream_or_object':sum(r['ownership']!='ggchord' for r in rows),
        'levels':dict(Counter(r['level'] for r in rows)),
        'limitation':'L2/L3 use explicit teaching-reviews.json judgments. Other L0/L1 rows are conservative evidence floors, not final depth certification. Runtime results recorded separately.'}
    (out/'teaching-summary-v0.16.7.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='source_sha256'},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
