#!/usr/bin/env python3
"""Generate English pages and paired content from one Han source and reviewed translations."""
import argparse
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WEB = ROOT / 'website'
CONFIG = WEB / 'content-policy.json'


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def blocks(text):
 lines=text.splitlines(keepends=True);out=[];i=0
 while i<len(lines):
  if not lines[i].strip():
   j=i+1
   while j<len(lines) and not lines[j].strip():j+=1
   out.append(('space',''.join(lines[i:j])));i=j;continue
  line=lines[i]
  if m:=re.match(r'^\s*(`{3,}|~{3,})',line):
   j=i+1
   while j<len(lines) and not re.match(r'^\s*'+re.escape(m[1][0])+r'{'+str(len(m[1]))+r',}\s*$',lines[j]):j+=1
   out.append(('code',''.join(lines[i:min(j+1,len(lines))])));i=j+1;continue
  if re.match(r'^\s*\{\{[^\n]+\}\}\s*$',line):kind='directive'
  elif re.match(r'^#{1,6} ',line):kind='heading'
  elif line.lstrip().startswith('|'):kind='table'
  elif re.match(r'^\s*(?:[-*+] |\d+[.)] )',line):kind='list'
  elif re.fullmatch(r'(?:\s*<[^>]+>)+\s*',line):kind='html'
  else:kind='paragraph'
  j=i+1
  if kind in ('paragraph','list'):
   while j<len(lines) and lines[j].strip() and not re.match(r'^(?:#{1,6} |\s*(?:[-*+] |\d+[.)] )|\s*\{\{|\s*\||\s*```|\s*<)',lines[j]):j+=1
  out.append((kind,''.join(lines[i:j])));i=j
 return out

# Protect executable snippets, shared directives, links, HTML and explicit IDs.
# Inline backticked prose in the generated parameter atlas remains translatable.
TOKENS = re.compile(r'`[^`\n]+`|\{#[^}\n]+\}|\{\{[^}\n]+\}\}|<[^>\n]+>|(?<=\]\()[^\)\n]*(?=\))')


def protected(text):
    values = []
    for match in TOKENS.finditer(text):
        value = match[0]
        if value.startswith('`') and re.search('[\u3400-\u9fff]', value):
            continue
        if value not in values:
            values.append(value)
    return values


def mask_target(source, target):
    for i, value in sorted(enumerate(protected(source)), key=lambda item: -len(item[1])):
        localized = value.replace('../assets/MIGRATION-Hans.md', '../assets/MIGRATION.md')
        if localized not in target:
            # Natural English may omit backticks around an inline identifier.
            if value.startswith('`') and localized.strip('`') in target:
                localized = localized.strip('`')
            else:
                raise ValueError('Translation lost shared token: ' + value)
        target = target.replace(localized, f'⟦shared:{i}⟧')
    return target


def unmask(source, target):
    for i, value in enumerate(protected(source)):
        marker = f'⟦shared:{i}⟧'
        if marker not in target:
            raise ValueError('Translation lost protected placeholder: ' + marker)
        value = value.replace('../assets/MIGRATION-Hans.md', '../assets/MIGRATION.md')
        target = target.replace(marker, value)
    if '⟦shared:' in target:
        raise ValueError('Unknown shared placeholder')
    return target


def validate_shape(source, target, kind):
    shape=[k for k,t in blocks(target) if k!='space']
    if shape!=[kind]:raise ValueError('Translation changes the shared block structure')
    if re.findall(r'\{#[^}]+\}',source)!=re.findall(r'\{#[^}]+\}',target):
        raise ValueError('Translated explicit anchors differ')
    if kind == 'heading':
        if re.match(r'^#+', source)[0] != re.match(r'^#+', target)[0]:
            raise ValueError('Translated heading level differs')
    if kind == 'table' and len(re.split(r'(?<!\\)\|', source)) != len(re.split(r'(?<!\\)\|', target)):
        raise ValueError('Translated table structure differs')
    # Block-level code/directives must always come from the canonical source.
    if any(token in target for token in ('```', '~~~', '{{code:', '{{prepare:')):
        raise ValueError('Translation contains executable structure')
    source_urls = re.findall(r'(?<=\]\()[^\)\n]*(?=\))', source)
    target_urls = re.findall(r'(?<=\]\()[^\)\n]*(?=\))', target)
    if sorted(u.replace('../assets/MIGRATION-Hans.md', '../assets/MIGRATION.md') for u in source_urls) != sorted(target_urls):
        raise ValueError('Translated links/images differ from source')


def translated_page(source, saved):
    if saved.get('source_sha256') != sha(source):
        raise ValueError('Han source changed; review translations before generating')
    entries = saved['blocks']
    used = set(); output = []; index = 0
    for kind, text in blocks(source):
        if kind == 'space':
            output.append(text); continue
        key = str(index); index += 1
        if kind in ('code', 'directive', 'html'):
            output.append(text); continue
        entry = entries.get(key)
        if entry is None:
            raise ValueError('Missing translation block ' + key)
        if entry.get('source_sha256') != sha(text) or entry.get('status') != 'reviewed':
            raise ValueError('Unreviewed/stale translation block ' + key)
        target = unmask(text, entry['translation'])
        validate_shape(text, target, kind)
        output.append(target); used.add(key)
    if used != set(entries):
        raise ValueError('Orphan translation blocks')
    return ''.join(output)


def resolve_content(source, saved, path=()):
    entries=saved['strings']; used=set()
    def visit(value, pointer):
        if isinstance(value, dict):
            if set(value)=={'zh'} and isinstance(value['zh'], str):
                key='/'.join(map(str,pointer)); entry=entries.get(key)
                if not entry or entry.get('source_sha256')!=sha(value['zh']) or entry.get('status')!='reviewed':
                    raise ValueError('Missing/stale structured translation: '+key)
                used.add(key)
                return dict(zh=value['zh'],en=entry['translation'])
            return {k:visit(v,(*pointer,k)) for k,v in value.items()}
        if isinstance(value,list):return [visit(v,(*pointer,i)) for i,v in enumerate(value)]
        return value
    canonical=json.dumps(source,ensure_ascii=False,sort_keys=True,separators=(',',':'))
    if saved.get('source_sha256')!=sha(canonical):raise ValueError('Structured source changed; translations need review')
    result=visit(source,path)
    if used!=set(entries):raise ValueError('Orphan structured translations')
    return result


def generated_outputs():
    policy=json.loads(CONFIG.read_text())
    if policy['source_language']!='zh' or policy['generated_language']!='en':raise ValueError('Unsupported content policy')
    outputs={}; counts={'pages':0,'extra_markdown':0,'structured_files':0,'translated_blocks':0,'translated_strings':0}
    for row in policy['markdown']:
        source=(WEB/row['source']).read_text()
        saved=json.loads((WEB/row['translations']).read_text())
        outputs[WEB/row['output']]=translated_page(source,saved)
        counts['pages' if row['source'].startswith('zh/') else 'extra_markdown']+=1
        counts['translated_blocks']+=len(saved['blocks'])
    for row in policy['structured']:
        source=json.loads((WEB/row['source']).read_text())
        saved=json.loads((WEB/row['translations']).read_text())
        result=resolve_content(source,saved)
        outputs[WEB/row['output']]=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
        counts['structured_files']+=1;counts['translated_strings']+=len(saved['strings'])
    return outputs,counts


def sync(check=True):
    outputs,counts=generated_outputs()
    stale=[str(path.relative_to(ROOT)) for path,text in outputs.items() if not path.exists() or path.read_text()!=text]
    if check and stale:raise ValueError('Generated language content was edited or is stale; regenerate from Han source and reviewed translations: '+', '.join(stale))
    if not check:
        for path,text in outputs.items():
            path.parent.mkdir(parents=True,exist_ok=True)
            if not path.exists() or path.read_text()!=text:path.write_text(text)
    return counts


def prepare():
    """Keep unchanged translations; new/changed source prose becomes an explicit draft."""
    policy=json.loads(CONFIG.read_text());drafts=0
    for row in policy['markdown']:
        source=(WEB/row['source']).read_text();path=WEB/row['translations']
        saved=json.loads(path.read_text()) if path.exists() else {'blocks':{}}
        previous=saved['blocks'];by_hash={v['source_sha256']:v for v in previous.values()}
        entries={};index=0
        for kind,text in blocks(source):
            if kind=='space':continue
            key=str(index);index+=1
            if kind in ('code','directive','html'):continue
            digest=sha(text)
            if previous.get(key,{}).get('source_sha256')==digest:entry=previous[key]
            elif digest in by_hash:entry=by_hash[digest]
            else:
                entry={'source_sha256':digest,'status':'needs_review',
                    'source_for_review':text,'translation':previous.get(key,{}).get('translation','')}
                drafts+=1
            entries[key]=entry
        fresh={'source_sha256':sha(source),'blocks':entries}
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(fresh,ensure_ascii=False,indent=2)+'\n')
    for row in policy['structured']:
        source=json.loads((WEB/row['source']).read_text());path=WEB/row['translations']
        saved=json.loads(path.read_text()) if path.exists() else {'strings':{}}
        previous=saved['strings'];entries={}
        def visit(value,pointer=()):
            nonlocal drafts
            if isinstance(value,dict):
                if set(value)=={'zh'} and isinstance(value['zh'],str):
                    key='/'.join(map(str,pointer));digest=sha(value['zh']);entry=previous.get(key,{})
                    if entry.get('source_sha256')!=digest:
                        entry={'source_sha256':digest,'status':'needs_review','source_for_review':value['zh'],
                            'translation':entry.get('translation','')};drafts+=1
                    entries[key]=entry
                else:
                    for key,item in value.items():visit(item,(*pointer,key))
            elif isinstance(value,list):
                for i,item in enumerate(value):visit(item,(*pointer,i))
        visit(source)
        payload=json.dumps(source,ensure_ascii=False,sort_keys=True,separators=(',',':'))
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps({'source_sha256':sha(payload),'strings':entries},ensure_ascii=False,indent=2)+'\n')
    return drafts


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--prepare',action='store_true',help='Mark changed Han prose as drafts; never accepts old translations automatically')
    args=parser.parse_args()
    if args.prepare:
        print("Prepared",prepare(),"translation drafts; review English and set status to reviewed before generation")
        return
    try:counts=sync(check=args.check)
    except (ValueError,KeyError) as error:raise SystemExit(str(error))
    print('PASS: Han source → reviewed English; '+', '.join(str(v)+' '+k for k,v in counts.items()))


if __name__=='__main__':main()
