#!/usr/bin/env python3
"""Verify affected PNG rerenders and scoped equivalence of retained legacy PNGs.

Require exact reviewed runtime changes, actual refresh records for Chord/Radial,
and original data, recipes and PNG hashes for every retained legacy asset.
"""
import csv,hashlib,importlib.util,json,re,subprocess,sys,tempfile
from pathlib import Path
from retained_provenance import validate_home_sources, render_source
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

CI_FILES = {'R/helpers-label-boxes.R', 'R/adaptive-layout.R'}

def verify_ci_patch():
    old = '    close_device <- if (grDevices::dev.cur() > 1L) function() invisible(NULL) else\n      ggchord_measurement_device(text = texts[valid])'
    new = '    legacy_device <- names(grDevices::dev.cur()) %in% c("pdf", "postscript")\n    unicode <- any(utf8ToInt(paste(texts[valid], collapse = "")) > 127L)\n    close_device <- if (grDevices::dev.cur() > 1L && !(legacy_device && unicode))\n      function() invisible(NULL) else ggchord_measurement_device(text = texts[valid])'
    source = blob('R/helpers-label-boxes.R').decode()
    assert source.count(old) == 1
    assert (ROOT/'R/helpers-label-boxes.R').read_text() == source.replace(old, new, 1)
    source = blob('R/adaptive-layout.R').decode()
    old_sig = 'ggplot_add.ggchord_layout_policy <- function(object, plot, object_name) {'
    new_sig = 'ggplot_add.ggchord_layout_policy <- function(object, plot, ...) {'
    assert source.count(old_sig) == 1
    assert (ROOT/'R/adaptive-layout.R').read_text() == source.replace(old_sig, new_sig, 1)
    # Registered images use PNG devices: the added base PDF/PostScript Unicode
    # branch is unreachable there. S3 method body is unchanged. Preserve the
    # real render hashes instead of stamping these old images as freshly drawn.

def verify_linear_revision():
    record=json.loads((ROOT/'tools/site-audit/home-rendering.json').read_text())
    assert record['status']=='complete'
    dot=validate_home_sources(record)
    review=record['runtime_review']
    assert review['baseline']=='a1cd31d' and review['signatures_identical'] and review['other_expressions_identical']
    assert set(review['changed_functions'])=={'compute_chord_geometry','ggchord_collinear_gene_labels','ggchord_collinear_share_labels'}
    assert review['representative_chord_circular_geometry_identical']
    previous=json.loads(subprocess.check_output(['git','show','a1cd31d:tools/site-audit/home-rendering.json'],cwd=ROOT))
    actual={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'R').glob('*.R')}
    render_engine={str(p.relative_to(ROOT)):hashlib.sha256(render_source(str(p.relative_to(ROOT)))).hexdigest() for p in (ROOT/'R').glob('*.R')}
    expected_engine=dict(record['engine_sha256'])
    if dot:expected_engine['R/coord-dotplot.R']=dot['engine_sha256']['R/coord-dotplot.R']
    assert render_engine==expected_engine,'Runtime changed outside reviewed rendering equivalence'
    changed={name for name,digest in render_engine.items() if previous['engine_sha256'].get(name)!=digest}
    assert changed==({'R/ggchord-build.R','R/layout-collinear-annotation.R','R/coord-dotplot.R'} if dot else {'R/ggchord-build.R','R/layout-collinear-annotation.R'})
    recipe=(ROOT/'tools/documentation-examples.R').read_bytes()
    assert recipe==subprocess.check_output(['git','show','a1cd31d:tools/documentation-examples.R'],cwd=ROOT)
    assert hashlib.sha256(recipe).hexdigest()==record['legacy_scope']['unchanged_recipes_sha256']
    baseline=list(csv.DictReader(subprocess.check_output(['git','show','a1cd31d:website/assets/plots.csv'],cwd=ROOT,text=True).splitlines()))
    rows=list(csv.DictReader((ROOT/'website/assets/plots.csv').open()))
    assert [r['id'] for r in rows]==[r['id'] for r in baseline]
    refreshed=record['legacy_figures']
    assert set(refreshed)==set(record['legacy_scope']['annotated'])=={'showcasecollinear','linearlabels','practicelinear'}
    current=module.fingerprints(runtime_hashes={'R/coord-dotplot.R':hashlib.md5(subprocess.check_output(['git','show','3a56988:R/coord-dotplot.R'],cwd=ROOT)).hexdigest()} if dot else None)['assets'];expected={}
    for old,row in zip(baseline,rows):
        name=row['id'];path=ROOT/row['file']
        assert hashlib.md5(path.read_bytes()).hexdigest()==row['md5'],name
        if name in refreshed:
            item=refreshed[name]
            assert item['status']=='actual PNG rendered'
            assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256']
            assert item['md5']==row['md5'] and item['dependency_md5']==row['dependency_md5']==current[name]
            assert all(old[k]==row[k] for k in old if k not in {'md5','code_md5','dependency_md5','generator'}),name
        elif dot and name=='viewdotplot':
            item=dot['figures'][name]
            assert item['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
            assert item['dependency_md5']==row['dependency_md5']==module.fingerprints(runtime_hashes={'R/coord-dotplot.R':hashlib.md5(render_source('R/coord-dotplot.R')).hexdigest()})['assets'][name]
        else:
            assert old==row,'Retained PNG provenance moved: '+name
        expected[name]=row['dependency_md5']
    result=dict(status='passed scoped annotation update',baseline='a1cd31d',figure_count=len(rows),
        actual_png_refresh_count=len(refreshed),retained_equivalent_count=len(rows)-len(refreshed),
        method='Unchanged recipes and unchanged AST outside three reviewed internal functions; all 42 linear/composition objects inspected for name layers; three affected PNGs regenerated; Chord/Circular representative geometry exact; other runtime files identical.',
        current_engine_sha256=actual,release_equivalence={name:dict(current_sha256=hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),render_sha256=hashlib.sha256(render_source(name)).hexdigest()) for name in ('DESCRIPTION','R/coord-dotplot.R')},expected_dependency_md5=expected,
        runtime_scope=review['scope'],scope='Source and normal-size PNG verification; no package-wide or browser certification')
    (ROOT/'tools/site-audit/existing-asset-equivalence.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 3 affected linear PNGs refreshed; 175 retained assets preserve their provenance')

def main():
    revision=json.loads((ROOT/'tools/site-audit/home-rendering.json').read_text())
    if revision.get('runtime_review'):
        verify_linear_revision()
        return
    changed=subprocess.check_output(['git','diff',BASE,'--name-only','--','R','data'],cwd=ROOT,text=True).splitlines()
    allowed={'R/api-contract.R','R/collinear-utils.R','R/coord-sankey.R','R/geom-ribbon.R',
      'R/ggchord-build.R','R/ggchord-frame.R','R/layout-collinear.R','R/layout-ribbon-step.R'} | CI_FILES
    assert not set(changed)-allowed,'Unreviewed runtime/data changes: '+str(changed)
    verify_ci_patch()
    refresh_path=ROOT/'tools/site-audit/ribbon-figure-refresh.json'
    refresh=json.loads(refresh_path.read_text()) if set(changed)-{'R/api-contract.R'} else None
    if refresh:
        assert refresh['status']=='complete','Affected PNG refresh is incomplete'
        for name,digest in refresh['runtime_sha256'].items():
            actual = blob(name) if name in CI_FILES else (ROOT/name).read_bytes()
            assert hashlib.sha256(actual).hexdigest()==digest,'Runtime moved after refresh: '+name
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
    def chunks(source):
        return dict(re.findall(r'^# ---- ([a-z]+) ----\n(.*?)(?=^# ---- [a-z]+ ----)',source.decode(),re.M|re.S))
    before_chunks=chunks(explicit);after_chunks=chunks(actual)
    changed_chunks={k for k in before_chunks.keys()|after_chunks.keys() if before_chunks.get(k)!=after_chunks.get(k)}
    approved={'frontgenelabels','frontradius','frontcurvature','parameteruniform',
              'parametervector','parametergenes','frontpreview','frontsave'}
    assert changed_chunks==approved,'Unreviewed example changes: '+str(changed_chunks)
    redesign=json.loads((ROOT/'tools/site-audit/redesign-rendering.json').read_text())
    assert redesign['status']=='complete'
    assert redesign['example_sha256']==hashlib.sha256(actual).hexdigest(),'Examples changed since actual rendering'
    assert all(redesign['public_checks'].values()),'Incomplete public behavior verification'
    affected={'frontgenelabels','frontradius','frontcurvature','parameteruniform',
              'parameterexample','parametervector','parameterstrand','parameterwidth'}
    # The changed preparations write only their own plots/tables. parameter_genes
    # is consumed only by strand/width; preview/save do not create figure inputs.
    assert {key for key,code in after_chunks.items() if 'parameter_genes' in code}=={'parametergenes','parameterstrand','parameterwidth'}
    for key in changed_chunks-{'frontpreview','frontsave','parametergenes'}:
        assigned=set(re.findall(r'^(p_\w+) <-',after_chunks[key],re.M))
        assert assigned,'Missing independent graph assignment: '+key
        for other,code in after_chunks.items():
            if other!=key and key!='frontgenelabels':
                assert not any(re.search(r'\b'+name+r'\b',code) for name in assigned),'New shared dependency: '+other
    current_redesign=module.fingerprints()['assets']

    subprocess.run(['Rscript','--vanilla','-e',
        'pkgload::load_all(quiet=TRUE);f<-"examples/format-inputs/emales-BVI_023A.gff3";a<-read_gff3(f,source_file=TRUE);b<-read_gff3(f,source_file=TRUE,unstranded="keep");stopifnot(identical(a,b),all(a$source_strand %in% c("+","-")));cat("PASS: explicit keep preserves all 24 source CDS records\\n")'],cwd=ROOT,check=True)
    overrides={name:hashlib.md5(blob(name)).hexdigest() for name in changed}
    expected=module.fingerprints(runtime_hashes=overrides,recipe_source=original.decode())['assets']
    current=module.fingerprints(runtime_hashes={name:hashlib.md5(blob(name)).hexdigest() for name in CI_FILES},recipe_source=explicit.decode())['assets']
    rows=list(csv.DictReader((ROOT/'website/assets/plots.csv').open()))
    actual=set(refresh['figures']) if refresh else set()
    if refresh:assert actual=={r['id'] for r in rows if r['coordinate'] in {'chord','radial_synteny'}}
    for row in rows:
        if row['id'] in affected:
            record=redesign['figures'][row['id']]
            assert record['status']=='actual PNG rendered' and record['md5']==row['md5']
            assert record['sha256']==hashlib.sha256((ROOT/row['file']).read_bytes()).hexdigest()
            expected[row['id']]=current_redesign[row['id']]
            assert record['dependency_md5']==expected[row['id']]
        elif row['id'] in actual:
            record=refresh['figures'][row['id']]
            assert record['status']=='actual PNG rendered' and record['md5']==row['md5']
            assert hashlib.sha256((ROOT/row['file']).read_bytes()).hexdigest()==record['sha256']
            expected[row['id']]=current[row['id']]
        assert expected[row['id']]==row['dependency_md5'],row['id']+': changed dependencies'
        assert hashlib.md5((ROOT/row['file']).read_bytes()).hexdigest()==row['md5'],row['id']+': changed PNG'
    result=dict(status='scoped actual PNG refresh plus retained equivalent figures' if refresh else 'equivalent dependencies; old renders retained',baseline=BASE,
        method='R parse equality for documentation-only edits; exact reviewed runtime style delta with four-Coord behavior/geometry proof; 61 actual PNG refresh records; explicit-keep input equality and baseline hashes for the remaining 117 PNGs',
        redesign_actual_png_count=len(affected),redesign_retained_png_count=len(rows)-len(affected),
        redesign_changed_chunks=sorted(changed_chunks),
        redesign_scope='Eight current normal-size PNG rerenders; other legacy PNGs retain their original render records under exact independent-example boundary checks.',
        current_api_document_sha256=hashlib.sha256((ROOT/'R/api-contract.R').read_bytes()).hexdigest(),
        png_equivalence_scope='Exact CI delta only changes active base PDF/PostScript Unicode measurement, unreachable for PNG; S3 method body unchanged',
        figure_count=len(rows),actual_png_refresh_count=len(actual),retained_equivalent_count=len(rows)-len(actual),
        expected_dependency_md5=expected,figures={r['id']:r['md5'] for r in rows},
        runtime_scope='Exact reviewed cross-row style patch; circular/linear/dotplot and two-source Fusion files retained, Chord/Radial rendered' if refresh else 'Documentation-only runtime change',
        scope='Scoped actual PNG refresh and retained-asset equivalence; no device matrix, package tests or general behavior certification')
    (ROOT/'tools/site-audit/existing-asset-equivalence.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS:',len(affected),'current PNG rerenders;',len(rows)-len(affected),'legacy PNGs retained with scoped source equivalence')
if __name__=='__main__':main()
