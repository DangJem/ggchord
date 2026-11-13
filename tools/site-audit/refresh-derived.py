#!/usr/bin/env python3
"""Refresh copied presentation/search assets after a verified native build."""
import hashlib,json,re,shutil,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from site_search import search_records

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    site=ROOT/'.site';snapshot_path=site/'source-snapshot.json'
    snapshot=json.loads(snapshot_path.read_text());old_snapshot=sha(snapshot_path)
    assert snapshot['status']=='verified unchanged across rendering and finalization'
    allowed={'website/assets/site.css','website/assets/site.js','tools/site_search.py'}
    changed={name for name,value in snapshot['render_inputs'].items() if sha(ROOT/name)!=value}
    assert changed and changed<=allowed,'Only derived presentation/search inputs may change: '+str(changed)
    manifest_path=site/'build-manifest.json';manifest=json.loads(manifest_path.read_text())
    for name in sorted(changed):
        path=ROOT/name
        if name=='tools/site_search.py':
            before={str(p):sha(p) for p in site.rglob('*.html')}
            for lang in ['en','zh']:
                target=site/lang/'search.json';old=json.loads(target.read_text())
                pages=[r for r in old if '#' not in r['url'] and 'section' not in r and 'queries' not in r]
                assert pages and len({r['url'] for r in pages})==len(pages)
                target.write_text(json.dumps(search_records(site,lang,pages),ensure_ascii=False)+'\n')
            assert before=={str(p):sha(p) for p in site.rglob('*.html')}
        else:
            shutil.copy2(path,site/'assets'/path.name)
            digest=hashlib.md5(path.read_bytes()).hexdigest()
            for html in site.rglob('*.html'):
                text=html.read_text();fresh=re.sub(r'(assets/'+re.escape(path.name)+r'\?v=)[^"\s>]+',lambda m:m[1]+digest,text)
                if fresh!=text:html.write_text(fresh)
        snapshot['render_inputs'][name]=sha(path)
        manifest['source_sha256'][name]=sha(path)
    snapshot.setdefault('derived_refreshes',[]).append(dict(
        previous_snapshot_sha256=old_snapshot,files={n:sha(ROOT/n) for n in sorted(changed)},
        scope='Copied CSS/JS cache keys or local search index from unchanged canonical HTML; no native content/examples rerendered'))
    snapshot_path.write_text(json.dumps(snapshot,indent=2)+'\n')
    manifest['derived_refreshes']=snapshot['derived_refreshes'];manifest['dirty']=True
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('PASS: guarded derived-asset refresh',sorted(changed))
if __name__=='__main__':main()
