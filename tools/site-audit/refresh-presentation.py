#!/usr/bin/env python3
"""Refresh copied CSS/JS only after a verified native build, with strict input guards."""
import hashlib,json,re,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    site=ROOT/'.site';snapshot_path=site/'source-snapshot.json'
    snapshot=json.loads(snapshot_path.read_text());old_snapshot=sha(snapshot_path)
    assert snapshot['status']=='verified unchanged across rendering and finalization'
    allowed={'website/assets/site.css','website/assets/site.js'}
    changed={name for name,value in snapshot['render_inputs'].items() if sha(ROOT/name)!=value}
    assert changed and changed<=allowed,'Only copied presentation assets may change: '+str(changed)
    manifest_path=site/'build-manifest.json';manifest=json.loads(manifest_path.read_text())
    for name in changed:
        path=ROOT/name;target=site/'assets'/path.name;shutil.copy2(path,target)
        digest=hashlib.md5(path.read_bytes()).hexdigest()
        for html in site.rglob('*.html'):
            text=html.read_text();fresh=re.sub(r'(assets/'+re.escape(path.name)+r'\?v=)[^"\s>]+',lambda m:m[1]+digest,text)
            if fresh!=text:html.write_text(fresh)
        snapshot['render_inputs'][name]=sha(path)
        manifest['source_sha256'][name]=sha(path)
    snapshot.setdefault('presentation_refreshes',[]).append(dict(
        previous_snapshot_sha256=old_snapshot,files={n:sha(ROOT/n) for n in sorted(changed)},
        scope='Copied presentation assets and their cache query strings; native HTML content and examples not rerendered'))
    snapshot_path.write_text(json.dumps(snapshot,indent=2)+'\n')
    manifest['presentation_refreshes']=snapshot['presentation_refreshes']
    manifest['dirty']=True
    manifest_path.write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print('PASS: presentation-only refresh',sorted(changed))

if __name__=='__main__':main()
