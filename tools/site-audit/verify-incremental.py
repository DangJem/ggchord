#!/usr/bin/env python3
"""Explicit scoped incremental-build acceptance; not a default package test."""
import copy, importlib.util, json, shutil, subprocess, sys, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
from site_cache import Cache,key,reference_inputs

def run(*args):subprocess.run(args,cwd=ROOT,check=True)
def main():
    findings={}
    with tempfile.TemporaryDirectory(prefix='ggchord-cache-verification-') as folder:
        folder=Path(folder);cache=Cache(folder);entries=list((cache.root/'reference').glob('*/entry.json'))
        assert entries,'A successfully verified cold build is required'
        current=reference_inputs();receipt=next(json.loads(p.read_text()) for p in entries if json.loads(p.read_text())['identity']['inputs']==current)
        identity=receipt['identity'];assert cache.restore('reference',identity,folder/'baseline')
        source=ROOT/'R/geom-residue.R'; original=source.read_bytes()
        try:
            source.write_bytes(original+b'\n# Incremental source invalidation probe.\n')
            changed=dict(identity,inputs=reference_inputs());assert changed!=identity
            assert not cache.restore('reference',changed,folder/'source-miss')
        finally:source.write_bytes(original)
        findings['actual_R_source_change']='miss'
        source=ROOT/'tools/website-dependencies.R';original=source.read_bytes()
        try:
            source.write_bytes(original+b'\n# Dependency configuration invalidation probe.\n')
            assert not cache.restore('reference',dict(identity,inputs=reference_inputs()),folder/'dependency-miss')
        finally:source.write_bytes(original)
        findings['actual_dependency_configuration_change']='miss'
        changed=copy.deepcopy(identity);changed['shell']['runtime']=key({'changed_installed_dependency_bytes':True})
        assert not cache.restore('reference',changed,folder/'runtime-miss')
        findings['simulated_installed_dependency_fingerprint_change']='miss'
        isolated=folder/'isolated';target=isolated/'reference'/key(identity)
        shutil.copytree(cache.root/'reference'/key(identity),target)
        damaged=next((target/'files').rglob('*.html'));damaged.write_bytes(b'damaged')
        cache.root=isolated
        assert not cache.restore('reference',identity,folder/'corrupt')
        findings['corrupt_output']='miss; no bytes restored'
    spec=importlib.util.spec_from_file_location('language_sync',ROOT/'tools/sync-website-languages.py')
    sync=importlib.util.module_from_spec(spec);spec.loader.exec_module(sync)
    canonical=ROOT/'website/zh/faq.md';translations=ROOT/'website/translations/en/faq.json';english=ROOT/'website/en/faq.md'
    originals={p:p.read_bytes() for p in (canonical,translations,english)}
    try:
        paragraph='This paragraph verifies incremental website rebuilding.\n'
        canonical.write_text(canonical.read_text().rstrip()+'\n\n增量构建验证：本段用于确认正文更新生效。\n')
        # Prepare only this page's reviewed translation, preserving every shared token.
        saved=json.loads(translations.read_text());by_hash={v['source_sha256']:v for v in saved['blocks'].values()}
        entries={};i=0
        for kind,text in sync.blocks(canonical.read_text()):
            if kind=='space':continue
            index=str(i);i+=1
            if kind in ('code','directive','html'):continue
            old=by_hash.get(sync.sha(text))
            entries[index]=old or dict(source_sha256=sync.sha(text),status='reviewed',translation=sync.mask_target(text,paragraph))
        translations.write_text(json.dumps(dict(source_sha256=sync.sha(canonical.read_text()),blocks=entries),ensure_ascii=False,indent=2)+'\n')
        sync.sync(check=False)
        run('python3','tools/build-website.py','--render-jobs','2','--render-batch-size','12')
        report=json.loads((ROOT/'.site/cache-report.json').read_text())
        refs=[r for r in report if r['kind']=='reference'];articles=[r for r in report if r['kind']=='article']
        assert len(refs)==2 and all(r['status']=='hit' for r in refs)
        assert sum(r['status']=='miss' for r in articles)==2,articles
        assert paragraph.strip() in (ROOT/'.site/en/faq.html').read_text()
        assert '增量构建验证' in (ROOT/'.site/zh/faq.html').read_text()
        run('python3','tools/check-documentation.py')
        findings['real_bilingual_prose_edit']=dict(API_groups='2 hits',article_misses=2,other_article_hits=len(articles)-2,updated_text_visible=True)
    finally:
        for p,data in originals.items():p.write_bytes(data)
    # Reconstruct the actual checkout for deployment; never publish the probe text.
    run('python3','tools/build-website.py','--render-jobs','2','--render-batch-size','12','--publication')
    run('python3','tools/site_cache.py','--verify-and-publish')
    assert 'This paragraph verifies incremental website rebuilding.' not in (ROOT/'.site/en/faq.html').read_text()
    (ROOT/'.site/incremental-verification.json').write_text(json.dumps(dict(status='passed',checks=findings),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(findings,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
