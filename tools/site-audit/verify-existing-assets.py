#!/usr/bin/env python3
"""Prove old PNG dependencies remain equivalent after documentation-only R edits.

Does not stamp generation hashes or claim a new render. Any executable change,
shared recipe change, source-data change or changed PNG fails this comparison.
"""
import csv,hashlib,importlib.util,json,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
spec=importlib.util.spec_from_file_location('fingerprints',ROOT/'tools/documentation-fingerprints.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
BASE='0f5dd96da57876da74cd1234c7de9cf5579c80b8'
def blob(name):return subprocess.check_output(['git','show',BASE+':'+name],cwd=ROOT)
def main():
    changed=subprocess.check_output(['git','diff',BASE,'--name-only','--','R','data'],cwd=ROOT,text=True).splitlines()
    if set(changed)-{'R/api-contract.R'}:raise SystemExit('Runtime/data changes require new figure verification: '+str(changed))
    old=blob('R/api-contract.R')
    with tempfile.TemporaryDirectory(prefix='ggchord-doc-equivalence-') as folder:
        file=Path(folder)/'api-contract.R';file.write_bytes(old)
        subprocess.run(['Rscript','--vanilla','-e','args<-commandArgs(TRUE);stopifnot(identical(parse(args[1],keep.source=FALSE),parse("R/api-contract.R",keep.source=FALSE)));cat("PASS: identical executable expressions; only documentation changed\\n")',str(file)],cwd=ROOT,check=True)
    actual=(ROOT/'tools/documentation-examples.R').read_bytes()
    assert actual==blob('tools/documentation-examples.R'),'Original shared plot recipes changed'
    overrides={'R/api-contract.R':hashlib.md5(old).hexdigest()}
    expected=module.fingerprints(runtime_hashes=overrides)['assets']
    rows=list(csv.DictReader((ROOT/'website/assets/plots.csv').open()))
    for row in rows:
        assert expected[row['id']]==row['dependency_md5'],row['id']+': changed dependencies'
        assert hashlib.md5((ROOT/row['file']).read_bytes()).hexdigest()==row['md5'],row['id']+': changed PNG'
    result=dict(status='equivalent dependencies; old renders retained',baseline=BASE,
        method='R parse(keep.source=FALSE) equality; unchanged original recipes/data/runtime and every stored dependency/PNG hash',
        current_api_document_sha256=hashlib.sha256((ROOT/'R/api-contract.R').read_bytes()).hexdigest(),
        figure_count=len(rows),figures={r['id']:r['md5'] for r in rows},
        scope='Not a fresh render, device matrix, package test or general behavior certification')
    (ROOT/'tools/site-audit/existing-asset-equivalence.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS:',len(rows),'existing PNGs retained with proved equivalent dependencies')
if __name__=='__main__':main()
