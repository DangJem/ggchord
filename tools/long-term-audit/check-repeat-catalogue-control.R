# Second real positive grouping control, with the same complete source profile.
source('tools/audit-source-snapshot.R')
source('tools/long-term-audit/feature-association-metrics.R')
audit<-audit_source_snapshot(additional=c('tools/long-term-audit/check-repeat-catalogue-control.R',
 'tools/long-term-audit/feature-association-metrics.R'))
pkgload::load_all(quiet=TRUE)
s<-plasmid_example_pTRIPZ
f<-find_common_features(s)
all_sites<-find_restriction_sites(s)
sites<-filter_restriction_sites(all_sites,set='reference',parent_set='commercial_nonredundant')
pos<-position_feature_stack(spacing=.085,base_position=position_plasmid())
p<-ggchord(s[c('accver','length')],feature_data=f)+geom_seq(seq_style='double')+
 geom_feature_plasmid(position=pos)+geom_feature_label_repel(position=pos,label_repeat='bracket',
 label_fit='shrink',label_min_size=2.5)+geom_seq_center_label()+
 geom_restriction_site(data=sites,label_size=2.5)+coord_circular()+theme_ggchord_plasmid()
path<-'tools/long-term-audit/v0.27.0-repeat-pTRIPZ-bracket.png'
ggchord_diagnostic_device(path,'png',10,8,110)
conditions<-character()
table<-withCallingHandlers(ggplot2::ggplotGrob(p),warning=function(w) {
 conditions<<-c(conditions,conditionMessage(w));invokeRestart('muffleWarning')})
grid::grid.newpage();grid::grid.draw(table)
l<-attr(table,'ggchord_render_layout')
visibility<-ggchord_label_visibility(l)
labels<-l$gene_labels
stopifnot(sum(labels$repeat_count>1L)==1L,
 labels$repeat_count[labels$repeat_count>1L]==6L,
 sum(visibility$component=='gene_label_repel' &
   visibility$source_text=='tet operator' &
   visibility$reason=='represented_by_shared_label')==6L,all(visibility$visible))
reference<-readRDS('/tmp/v0.27.0-compact-final-pTRIPZ.rds')
for(component in c('gene_poly','restriction_site')) {
 collect<-function(layout) {
  parts<-lapply(layout$layer_inputs,function(owner)owner[[component]])
  ggchord_rbind_fill(parts) %||% data.frame()
 }
 a<-collect(reference);b<-collect(l)
 fields<-intersect(c('accver','start','end','strand','anno','position','enzyme','segments'),names(a))
 stopifnot(identical(a[fields],b[fields]))
}
result<-list(file=path,source_features=nrow(f),selected_sites=nrow(sites),
 visibility=visibility,converged=attr(table,'ggchord_render_context')$converged,
 warnings=conditions,quality=association_text_quality(table)$summary)
grDevices::dev.off()
audit_assert_source_snapshot(audit)
jsonlite::write_json(c(result,list(runtime_sha256=as.list(audit$engine))),
 'tools/long-term-audit/v0.27.0-repeat-catalogue-control-results.json',pretty=TRUE,auto_unbox=TRUE)
cat('pTRIPZ six consecutive operators share one label with all source rows retained\n')
