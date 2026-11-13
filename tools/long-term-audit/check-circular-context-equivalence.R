# Compare the unchanged circular branch after a non-circular/empty-input fix.
source('tools/audit-source-snapshot.R')
audit<-audit_source_snapshot(additional='tools/long-term-audit/check-circular-context-equivalence.R')
pkgload::load_all(quiet=TRUE)
previous<-readRDS('/tmp/ggchord-pre-compatibility-functions.rds')
reference<-jsonlite::read_json('tools/long-term-audit/v0.27.0-compact-final-results.json')
stopifnot(all(vapply(names(reference$runtime_sha256),function(path)
 identical(previous$engine_sha256[[path]],reference$runtime_sha256[[path]]),logical(1))))
paths<-setNames(paste0('/tmp/v0.27.0-compact-final-',names(reference$cases),'.rds'),names(reference$cases))
for(id in c('pTRE_Tight_BI','pCAMBIA1300','pETDuet_1')) for(style in c('individual','bracket')) {
 name<-paste(id,style,sep='-')
 paths[name]<-paste0('/tmp/ggchord-compact-',name,'.rds')
}
ggchord_diagnostic_device('/tmp/ggchord-circular-equivalence.png','png',10,8,110)
observations<-list()
for(name in names(paths)) {
 l<-readRDS(paths[[name]])
 count<-0L
 for(owner in names(l$layer_geometry)) {
  geometry<-l$layer_geometry[[owner]]$gene_label_repel
  if(!is.data.frame(geometry) || !nrow(geometry))next
  sub<-l$layer_layouts[[owner]]
  stopifnot(isTRUE(sub$circular))
  before<-previous$repel(sub)
  after<-ggchord_repel_geometry(sub)
  stopifnot(identical(before,after))
  count<-count+1L
 }
 stopifnot(count>0L)
 observations[[name]]<-list(equal=TRUE,label_owners=count)
}
grDevices::dev.off()
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(scope='Exact circular label geometry equality; not a new full device rendering',
 cases=observations,rendered_baseline_sha256=reference$runtime_sha256,
 current_runtime_sha256=as.list(audit$engine)),
 'tools/long-term-audit/v0.27.0-circular-context-equivalence-results.json',pretty=TRUE,auto_unbox=TRUE)
cat('Circular label geometry unchanged across the 13 full maps and six style controls\n')
