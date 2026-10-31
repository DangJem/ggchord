# Selected plot/layer/cache/render ownership risks. No release checks.
pkgload::load_all(quiet=TRUE)
ns <- asNamespace('ggchord')
private <- function(name) get(name,ns)
results <- list()
record <- function(case,detail) {
  results[[length(results)+1L]] <<- data.frame(case=case,result='PASS',detail=detail)
  cat(case, ': PASS\n')
}
fail <- function(expr,pattern) {
  e<-tryCatch(force(expr),error=identity)
  stopifnot(inherits(e,'error'),grepl(pattern,conditionMessage(e),ignore.case=TRUE))
}
local({
  path<-tempfile(fileext='.pdf');grDevices::pdf(path,width=8,height=6)
  on.exit({grDevices::dev.off();unlink(path)},add=TRUE)
  s<-data.frame(accver=c('A','B'),length=1000)
  g<-data.frame(accver='A',start=100,end=300,strand='+',anno='one')
  snapshot<-function(p) lapply(p$layers,function(l) list(data=l$data,mapping=l$mapping,
    params=l$ggchord_params,input=l$ggchord_input_data,input_mapping=l$ggchord_input_mapping,
    fixed=l$aes_params,geom=l$geom_params,legend=l$show.legend,id=l$ggchord_layer_id))
  shared<-geom_gene(aes(gene_fill=anno),data=g)
  p<-ggchord(s)+geom_seq()+shared+shared+scale_gene_fill_manual(values=c(one='#123456'))
  before<-snapshot(p)
  invisible(ggplot2::ggplot_build(p))
  out<-export_ggchord_layout(p,include=c('seq','gene'),original_data=TRUE)
  stopifnot(nrow(out$seq)>0,nrow(out$gene)>0,
    length(unique(out$gene$layer_id))==2L,identical(snapshot(p),before))
  again<-export_ggchord_layout(p,include=c('seq','gene'),original_data=TRUE)
  stopifnot(identical(out,again))
  record('F48_build_export_and_reused_layer','Build then export retains rows; repeated layer objects get independent IDs without source mutation')
  sibling<-p+coord_collinear()
  copy<-private('ggchord_render_plot_copy')(p)
  invisible(ggplot2::ggplot_build(sibling));invisible(ggplot2::ggplot_build(copy))
  stopifnot(identical(out,export_ggchord_layout(p,include=c('seq','gene'),original_data=TRUE)),
    identical(snapshot(p),before))
  record('F05_sibling_and_render_copy','Interleaved sibling/copy builds preserve base geometry, source data, mappings and styles')
  sites<-data.frame(accver='A',position=c(100,300,600),enzyme=c('e1','e2','e3'))
  circular<-ggchord(s[1,])+geom_seq()+geom_restriction_site(data=sites)+coord_circular()
  original<-snapshot(circular)
  invisible(get_chord_layout(circular))
  stopifnot(identical(original,snapshot(circular)))
  record('F48_derived_parameters','Restriction backbone offsets and circular styles remain build-local')
  ctx<-private('ggchord_render_context');with_ctx<-private('ggchord_with_render_context')
  initial<-ctx()
  with_ctx(list(status='outer'),{
    fail(with_ctx(list(status='inner'),stop('intentional inner error')),'intentional')
    stopifnot(identical(ctx()$status,'outer'))
  })
  stopifnot(identical(ctx(),initial))
  record('F26_nested_error_context','Nested errors restore outer context and final exit restores the initial context')
  values<-data.frame(accver='A',start=100,end=200,value=5)
  physical<-ggchord(s[1,])+geom_seq()+geom_track_bar(aes(track_value=value),
    data=values,position=position_track('signal'))+coord_collinear()
  original<-snapshot(physical)
  baseline<-export_ggchord_layout(physical,include='track')
  cache<-get_chord_layout(physical,build=FALSE)
  for (size in list(c(6,4),c(9,5))) {
    d<-diagnose_ggchord(physical,width=size[1],height=size[2],device='png')
    stopifnot(identical(d$measurement$status,'measured'),
      all(is.finite(d$measurement$panel_mm)),nrow(d$track_registry)>0,
      all(d$track_registry$measurement_status==d$measurement$geometry_metric_status))
    if (isTRUE(d$measurement$converged)) {
      mm<-with(d$track_registry,(normal_max-normal_min)/d$measurement$units_per_mm[['y']])
      stopifnot(all(abs(mm-4)<.15))
    }
  }
  stopifnot(identical(snapshot(physical),original),identical(ctx(),initial),
    identical(get_chord_layout(physical,build=FALSE),cache),
    identical(export_ggchord_layout(physical,include='track'),baseline))
  record('F22_F26_target_size','Two PNG target sizes measure physical bands consistently; original layer inputs/cache and dynamic context stay unchanged')
  preview<-view_ggchord(physical,width=6,height=4,device='png',viewer='none')
  stopifnot(is.character(preview),file.exists(preview),identical(snapshot(physical),original),
    identical(get_chord_layout(physical,build=FALSE),cache))
  unlink(c(preview,sub('\\.png$','.html',preview)))
  record('F22_preview','Explicit-size static preview returns an image path and preserves source/cache')
  previous<-grDevices::dev.cur()
  fail(diagnose_ggchord(physical+scale_track_value_continuous('signal',limits=c(1,10)),
    width=6,height=4,device='png'),'baseline')
  stopifnot(grDevices::dev.cur()==previous,identical(ctx(),initial))
  record('F26_failed_device_build','Failed diagnostic build restores the active device/context')
  # Even ordinary fixed labels use device metrics, not only physical tracks.
  labels<-ggchord(s)+geom_seq()+geom_gene_label(data=g)
  invisible(get_chord_layout(labels))
  second<-tempfile(fileext='.pdf');grDevices::pdf(second,width=5,height=4)
  tryCatch({
    fail(get_chord_layout(labels,build=FALSE),'not found')
    stopifnot(nrow(get_chord_layout(labels)$seq_arcs[[1]])>0)
  },finally={grDevices::dev.off();unlink(second)})
  record('F26_device_invalidation','Changing device/size invalidates label-layout estimates before cache-only access')
})
out<-do.call(rbind,results)
write.csv(out,'design-roadmap/audit-generated/engine-state-validation-v0.16.3.csv',row.names=FALSE)
print(out,row.names=FALSE)
