#!/usr/bin/env python3
"""Explicit scoped incremental-build acceptance; not a default package test."""
import copy, importlib.util, json, shutil, subprocess, sys, tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
from site_cache import Cache,key,reference_inputs

def run(*args):subprocess.run(args,cwd=ROOT,check=True)
def verify_reuse(report, pages, edited):
    languages = ('en', 'zh')
    references = [row for row in report if row['kind'] == 'reference']
    assert sorted((row['language'], row['status']) for row in references) == [
        (lang, 'hit') for lang in languages], 'Both native reference groups must be reused'
    rendered = {p['slug'] for p in pages if p['kind'] != 'alias' and p['slug'] != 'reference'}
    assert edited in rendered, f'{edited} is not a rendered handbook article'
    expected = sorted((lang, slug, 'miss' if slug == edited else 'hit')
                      for lang in languages for slug in rendered)
    for kind in ('markdown', 'article'):
        actual = sorted((r['language'], r['page'], r['status']) for r in report if r['kind'] == kind)
        assert actual == expected, f'{kind} reuse differs from the current bilingual article inventory'
    return len(expected) - len(languages)


def main():
    findings={}
    out = ROOT / '.site'
    metadata = json.loads((out / 'build-metadata.json').read_text())
    mode = metadata['deployment']
    if mode not in ('local preview; not deployed', 'prepared for authorized deployment; builder does not deploy'):
        raise ValueError('Unknown original website build mode: ' + mode)
    build = ['python3', 'tools/build-website.py', '--render-jobs', '2', '--render-batch-size', '12']
    if mode.startswith('prepared for authorized deployment'):
        build.append('--publication')
    pages = json.loads((ROOT / 'website/pages.json').read_text())
    baseline = json.loads((out / 'cache-report.json').read_text())
    entry = next((row for row in baseline if row['kind'] == 'reference' and row['language'] == 'en'), None)
    if not entry:
        raise ValueError('A verified cache-enabled reference build is required')
    with tempfile.TemporaryDirectory(prefix='ggchord-cache-verification-') as folder:
        folder = Path(folder)
        cache = Cache(folder)
        receipt = json.loads((cache.root / 'reference' / entry['key'] / 'entry.json').read_text())
        identity = receipt['identity']
        assert identity['inputs'] == reference_inputs()
        assert cache.restore('reference', identity, folder/'baseline')
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
    slug = 'installation'
    canonical = ROOT/'website/en'/f'{slug}.md'
    translations = ROOT/'website/translations/zh'/f'{slug}.json'
    han = ROOT/'website/zh'/f'{slug}.md'
    originals={p:p.read_bytes() for p in (canonical,translations,han)}
    try:
        paragraph='This paragraph verifies incremental website rebuilding.\n'
        source = canonical.read_text().rstrip() + '\n\n' + paragraph
        reviewed = han.read_text().rstrip() + '\n\n增量构建验证：本段用于确认正文更新生效。\n'
        record = sync.reviewed_page(source, reviewed)
        canonical.write_text(source)
        translations.write_text(sync.json_text(record))
        sync.sync(check=False)
        run(*build)
        report=json.loads((out/'cache-report.json').read_text())
        other_articles = verify_reuse(report, pages, slug)
        assert paragraph.strip() in (ROOT/'.site/en'/f'{slug}.html').read_text()
        assert '增量构建验证' in (ROOT/'.site/zh'/f'{slug}.html').read_text()
        run('python3','tools/check-documentation.py')
        findings['real_bilingual_prose_edit']=dict(API_groups='2 hits',article_misses=2,other_article_hits=other_articles,updated_text_visible=True)
    finally:
        for p,data in originals.items():p.write_bytes(data)
    # Reconstruct the actual checkout for deployment; never publish the probe text.
    run(*build)
    run('python3','tools/site_cache.py','--verify-and-publish')
    assert 'This paragraph verifies incremental website rebuilding.' not in (ROOT/'.site/en'/f'{slug}.html').read_text()
    (ROOT/'.site/incremental-verification.json').write_text(json.dumps(dict(status='passed',checks=findings),ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(findings,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
