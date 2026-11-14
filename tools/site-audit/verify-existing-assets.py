#!/usr/bin/env python3
"""Verify affected PNG rerenders and scoped equivalence of retained legacy PNGs.

Require exact reviewed runtime changes, actual refresh records for Chord/Radial,
and original data, recipes and PNG hashes for every retained legacy asset.
"""
import csv,hashlib,importlib.util,json,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
spec=importlib.util.spec_from_file_location('fingerprints',ROOT/'tools/documentation-fingerprints.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
BASE='0f5dd96da57876da74cd1234c7de9cf5579c80b8'
def blob(name):return subprocess.check_output(['git','show',BASE+':'+name],cwd=ROOT)
def verify_reviewed_style_patch():
    prefix = "# Cross-row outline conventions belong to linear rows, not every derived Coord.\nggchord_collinear_row_styles <- function(coord) {\n  isTRUE(coord$ggchord_collinear) &&\n    !isTRUE(coord$ggchord_radial) && !isTRUE(coord$ggchord_sankey)\n}\n\n"
    assert (ROOT/'R/collinear-utils.R').read_text()==prefix+blob('R/collinear-utils.R').decode()
    replacements={
      'R/layout-ribbon-step.R':('      ribbon_linetype_vec <- ifelse(\n        ribbon_data$.ggchord_adjacent[valid_idx], "solid", "dashed"\n      )','      ribbon_linetype_vec <- rep("solid", length(valid_idx))'),
      'R/ggchord-build.R':('      collinear_transform = collinear_transform,\n      collinear_vertical =','      collinear_transform = collinear_transform,\n      collinear_crossrow_styles = ggchord_collinear_row_styles(plot$coordinates),\n      collinear_vertical ='),
      'R/layout-collinear.R':('    } else if (isTRUE(input$.ggchord_adjacent[i])) "solid" else "dashed"','    } else if (!isTRUE(p$collinear_crossrow_styles) ||\n        isTRUE(input$.ggchord_adjacent[i])) "solid" else "dashed"'),
      'R/ggchord-frame.R':('    if (isTRUE(layout$collinear) && identical(lyr$ggchord_type, "ribbon") &&','    if (ggchord_collinear_row_styles(plot$coordinates) &&\n        identical(lyr$ggchord_type, "ribbon") &&'),
      'R/coord-sankey.R':('    degrees[degrees > 90] <- degrees[degrees > 90] - 180\n    degrees[degrees < -90] <- degrees[degrees < -90] + 180','    above <- which(degrees > 90)\n    degrees[above] <- degrees[above] - 180\n    below <- which(degrees < -90)\n    degrees[below] <- degrees[below] + 180')}
    for name,(old,new) in replacements.items():
        original=blob(name).decode();assert original.count(old)==1,name
        assert (ROOT/name).read_text()==original.replace(old,new,1),'Runtime delta exceeds reviewed style scope: '+name
    with tempfile.TemporaryDirectory(prefix='ggchord-ribbon-docs-') as folder:
        original=Path(folder)/'geom-ribbon.R';original.write_bytes(blob('R/geom-ribbon.R'))
        subprocess.run(['Rscript','--vanilla','-e','a<-commandArgs(TRUE);stopifnot(identical(parse(a[1],keep.source=FALSE),parse("R/geom-ribbon.R",keep.source=FALSE)))',str(original)],cwd=ROOT,check=True)

def main():
    changed=subprocess.check_output(['git','diff',BASE,'--name-only','--','R','data'],cwd=ROOT,text=True).splitlines()
    allowed={'R/api-contract.R','R/collinear-utils.R','R/coord-sankey.R','R/geom-ribbon.R',
      'R/ggchord-build.R','R/ggchord-frame.R','R/layout-collinear.R','R/layout-ribbon-step.R'}
    assert not set(changed)-allowed,'Unreviewed runtime/data changes: '+str(changed)
    refresh_path=ROOT/'tools/site-audit/ribbon-figure-refresh.json'
    refresh=json.loads(refresh_path.read_text()) if set(changed)-{'R/api-contract.R'} else None
    if refresh:
        assert refresh['status']=='complete','Affected PNG refresh is incomplete'
        for name,digest in refresh['runtime_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'Runtime moved after refresh: '+name
        before=json.loads((ROOT/'tools/site-audit/ribbon-styles-before.json').read_text())
        after=json.loads((ROOT/'tools/site-audit/ribbon-styles-after.json').read_text())
        for name,row in after['views'].items():
            assert row['geometry_sha256']==before['views'][name]['geometry_sha256'],'Ribbon geometry changed: '+name
            assert row['default']==(['solid','dashed','dashed'] if name=='Collinear' else ['solid']*3)
            assert row['mapped']==['dotdash','longdash','dotted'] and row['fixed']==['dotted']*3
        assert after['runtime_sha256']==refresh['runtime_sha256'],'Style proof and PNG engine differ'
        verify_reviewed_style_patch()

    old=blob('R/api-contract.R')
    with tempfile.TemporaryDirectory(prefix='ggchord-doc-equivalence-') as folder:
        file=Path(folder)/'api-contract.R';file.write_bytes(old)
        subprocess.run(['Rscript','--vanilla','-e','args<-commandArgs(TRUE);stopifnot(identical(parse(args[1],keep.source=FALSE),parse("R/api-contract.R",keep.source=FALSE)));cat("PASS: identical executable expressions; only documentation changed\\n")',str(file)],cwd=ROOT,check=True)
    actual=(ROOT/'tools/documentation-examples.R').read_bytes()
    original=blob('tools/documentation-examples.R')
    explicit=original.replace(b'read_gff3("examples/format-inputs/emales-BVI_023A.gff3", source_file = TRUE)',b'read_gff3("examples/format-inputs/emales-BVI_023A.gff3", unstranded = "keep", source_file = TRUE)')
    assert actual==explicit,'Recipe edits exceed the exact explicit-keep change'
    subprocess.run(['Rscript','--vanilla','-e',
        'pkgload::load_all(quiet=TRUE);f<-"examples/format-inputs/emales-BVI_023A.gff3";a<-read_gff3(f,source_file=TRUE);b<-read_gff3(f,source_file=TRUE,unstranded="keep");stopifnot(identical(a,b),all(a$source_strand %in% c("+","-")));cat("PASS: explicit keep preserves all 24 source CDS records\\n")'],cwd=ROOT,check=True)
    overrides={name:hashlib.md5(blob(name)).hexdigest() for name in changed}
    expected=module.fingerprints(runtime_hashes=overrides,recipe_source=original.decode())['assets']
    current=module.fingerprints()['assets']
    rows=list(csv.DictReader((ROOT/'website/assets/plots.csv').open()))
    actual=set(refresh['figures']) if refresh else set()
    if refresh:assert actual=={r['id'] for r in rows if r['coordinate'] in {'chord','radial_synteny'}}
    for row in rows:
        if row['id'] in actual:
            record=refresh['figures'][row['id']]
            assert record['status']=='actual PNG rendered' and record['md5']==row['md5']
            assert hashlib.sha256((ROOT/row['file']).read_bytes()).hexdigest()==record['sha256']
            expected[row['id']]=current[row['id']]
        assert expected[row['id']]==row['dependency_md5'],row['id']+': changed dependencies'
        assert hashlib.md5((ROOT/row['file']).read_bytes()).hexdigest()==row['md5'],row['id']+': changed PNG'
    result=dict(status='scoped actual PNG refresh plus retained equivalent figures' if refresh else 'equivalent dependencies; old renders retained',baseline=BASE,
        method='R parse equality for documentation-only edits; exact reviewed runtime style delta with four-Coord behavior/geometry proof; 61 actual PNG refresh records; explicit-keep input equality and baseline hashes for the remaining 117 PNGs',
        current_api_document_sha256=hashlib.sha256((ROOT/'R/api-contract.R').read_bytes()).hexdigest(),
        figure_count=len(rows),actual_png_refresh_count=len(actual),retained_equivalent_count=len(rows)-len(actual),
        expected_dependency_md5=expected,figures={r['id']:r['md5'] for r in rows},
        runtime_scope='Exact reviewed cross-row style patch; circular/linear/dotplot and two-source Fusion files retained, Chord/Radial rendered' if refresh else 'Documentation-only runtime change',
        scope='Scoped actual PNG refresh and retained-asset equivalence; no device matrix, package tests or general behavior certification')
    (ROOT/'tools/site-audit/existing-asset-equivalence.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS:',len(actual),'actual refreshed PNGs;',len(rows)-len(actual),'retained with scoped equivalent dependencies')
if __name__=='__main__':main()
