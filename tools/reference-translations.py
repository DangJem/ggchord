#!/usr/bin/env python3
"""Versioned Han prose for canonical Rd; no code or API signature translation."""
import argparse, hashlib, json, re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'website/reference-zh'
PROSE={'title','description','details','value','format','note','source','references','author','seealso'}
PROTECTED={'code','verb','preformatted','link','url','eqn','deqn','kbd','samp','env','option','file','email','acronym','donttest','dontrun','testonly'}

def group(text, start):
    assert text[start]=='{'
    depth=1;i=start+1
    while depth:
        if text[i]=='\\': i+=2;continue
        depth += (text[i]=='{')-(text[i]=='}');i+=1
    return start+1,i-1,i

def commands(text):
    i=0
    while i<len(text):
        match=re.search(r'\\([A-Za-z]+)',text[i:])
        if not match:return
        a=i+match.start();name=match[1];j=i+match.end()
        if j<len(text) and text[j]=='[':
            end=text.index(']',j);j=end+1
        if j>=len(text) or text[j]!='{':i=j;continue
        b,c,end=group(text,j)
        yield name,a,b,c,end
        i=end

def fields(text):
    result=[]
    for name,a,b,c,end in commands(text):
        if name in PROSE:result.append((name,b,c))
        elif name=='arguments':
            for n,x,y,z,e in commands(text[b:c]):
                if n!='item':continue
                _,v,after=group(text[b:c],e)
                result.append(('argument:'+text[b+y:b+z],b+e+1,b+v))
        elif name=='section':
            result.append(('section-title',b,c))
            d,e,after=group(text,end)
            result.append(('section',d,e))
    return result

def protect(text):
    spans=[]
    for name,a,b,c,end in commands(text):
        if name in PROTECTED:spans.append((a,end))
        elif name=='href':
            # URL is fixed; visible labels are prose.
            spans.append((a,end))
        else:
            # Structural Rd macros may contain inline code/link nodes.
            for x,y in protected_spans(text[b:c]):spans.append((b+x,b+y))
    tokens=[];cursor=0;parts=[]
    for a,b in sorted(spans):
        parts.append(text[cursor:a]);parts.append('⟦'+str(len(tokens))+'⟧');tokens.append(text[a:b]);cursor=b
    parts.append(text[cursor:]);return re.sub(r'\s+',' ',''.join(parts)).strip(),tokens

def protected_spans(text):
    for name,a,b,c,end in commands(text):
        if name in PROTECTED or name=='href':yield a,end
        else:
            for x,y in protected_spans(text[b:c]):yield b+x,b+y

def inventory(path):
    text=path.read_text();entries=[]
    for label,b,c in fields(text):
        source,tokens=protect(text[b:c])
        # Pure names, formulae and links need no prose translation.
        residue=re.sub(r'⟦\d+⟧','',source)
        if not re.search(r'[A-Za-z]{2}',residue):continue
        entries.append(dict(key=str(len(entries)),field=label,source=source,tokens=tokens,start=b,end=c))
    return text,entries

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--export',action='store_true');ap.add_argument('--check',action='store_true');ap.add_argument('--output');ap.add_argument('--show',nargs='*');args=ap.parse_args()
    DEST.mkdir(exist_ok=True);missing=[];count=0
    for path in sorted((ROOT/'man').glob('*.Rd')):
        text,entries=inventory(path);target=DEST/(path.stem+'.json');digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if args.show is not None:
            if not args.show or path.stem in args.show:
                print('\n## '+path.stem)
                for e in entries:print(e['key']+' '+e['field']+': '+e['source'])
            continue
        if args.export:
            if not target.exists():target.write_text(json.dumps(dict(source_sha256=digest,translations={}),ensure_ascii=False,indent=2)+'\n')
            continue
        saved=json.loads(target.read_text()) if target.exists() else {}
        if saved.get('source_sha256')!=digest:missing.append(path.stem+': stale source hash')
        translations=saved.get('translations',{});replacements=[]
        for e in entries:
            value=translations.get(e['key'])
            if not value:missing.append(path.stem+': missing '+e['key']+' '+e['field']);continue
            if sorted(re.findall(r'⟦\d+⟧',value))!=sorted(re.findall(r'⟦\d+⟧',e['source'])):
                missing.append(path.stem+': changed inline tokens '+e['key']);continue
            if not re.search('[\u3400-\u9fff]',value):missing.append(path.stem+': no Han prose '+e['key'])
            value=re.sub(r'(?<!\\)%', r'\\%', value)
            for i,token in enumerate(e['tokens']):value=value.replace('⟦'+str(i)+'⟧',token)
            replacements.append((e['start'],e['end'],value));count+=1
        if set(translations)-{e['key'] for e in entries}:missing.append(path.stem+': orphan translations')
        if args.output:
            for b,c,value in reversed(replacements):text=text[:b]+value+text[c:]
            out=Path(args.output);out.mkdir(parents=True,exist_ok=True);(out/path.name).write_text(text)
    if missing:raise SystemExit('\n'.join(missing))
    if args.show is None:print(f'PASS: {count} translated prose blocks; code, identifiers and links preserved.')
if __name__=='__main__':main()
