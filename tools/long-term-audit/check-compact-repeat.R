# Same catalogue inputs and physical canvas; explicit repeat style comparison.
source('tools/audit-source-snapshot.R')
audit<-audit_source_snapshot(additional='tools/long-term-audit/check-compact-repeat.R')
pkgload::load_all(quiet=TRUE)
results<-list()
ids <- commandArgs(trailingOnly = TRUE)
if (!length(ids)) ids <- c('pTRE_Tight_BI','pCAMBIA1300','pETDuet_1')
for(id in ids) for(style in c('individual','bracket')) {
 e<-new.env();data(list=paste0('plasmid_example_',id),package='ggchord',envir=e)
 s<-e[[paste0('plasmid_example_',id)]];f<-find_common_features(s)
 pos<-position_feature_stack(spacing=.085,base_position=position_plasmid())
 p<-ggchord(s[c('accver','length')],feature_data=f)+geom_seq(seq_style='double')+
  geom_feature_plasmid(position=pos)+geom_feature_label_repel(position=pos,
   label_repeat=style,label_fit='shrink',label_min_size=2.5)+
  geom_seq_center_label()+coord_circular()+theme_ggchord_plasmid()
 name<-paste(id,style,sep='-');cat(name,'\n');flush.console()
 path<-paste0('tools/long-term-audit/v0.27.0-compact-',name,'.png')
 ggchord_diagnostic_device(path,'png',10,8,110)
 warnings<-character()
 table<-withCallingHandlers(ggplot2::ggplotGrob(p),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 grid::grid.newpage();grid::grid.draw(table)
 l<-attr(table,'ggchord_render_layout')
 saveRDS(l,paste0('/tmp/ggchord-compact-',name,'.rds'))
 results[[name]]<-list(file=path,visibility=ggchord_label_visibility(l),
  converged=attr(table,'ggchord_render_context')$converged,warnings=warnings,
  fit=ggchord_label_fit_diagnostics(l),
  polygons=length(unique(l$gene_polys$group)),source_features=nrow(f))
 if (any(l$gene_labels$feature_label_mode %in% 'external')) {
   built <- attr(table,'ggchord_render_built')
   drawn <- ggchord_rbind_fill(lapply(built$data,function(d) {
     if (!all(c('.component','.feature_internal_exit') %in% names(d)))return(NULL)
     d[d$.component %in% 'segment' & d$.feature_internal_exit %in% TRUE,,drop=FALSE]
   }))
   stopifnot(is.data.frame(drawn),nrow(drawn)>0L)
 }
 stopifnot(all(results[[name]]$visibility$visible),
   length(unique(l$gene_polys$source_row))==nrow(f))
 if(id=='pTRE_Tight_BI' && style=='bracket') {
   labels <- l$gene_labels
   stopifnot(sum(labels$repeat_count>1L)==1L,
     labels$repeat_count[labels$repeat_count>1L]==7L,
     sum(results[[name]]$visibility$reason=='represented_by_shared_label')==7L,
     sum(results[[name]]$fit$label_fit_reason=='represented_by_shared_label')==7L)
   rendered <- l$layer_geometry[[names(l$layer_geometry)[3L]]]$gene_label_repel
   # Remaining leaders must resolve to a visible stable text/source group.
   segments <- rendered[rendered$.component=='segment',,drop=FALSE]
   texts <- rendered[rendered$.component=='text',,drop=FALSE]
   stopifnot(all(segments$group %in% texts$group),
     all(segments$source_row %in% texts$source_row))
 }
 grDevices::dev.off()
}
if ('pTRE_Tight_BI' %in% ids) {
 e<-new.env();data('plasmid_example_pTRE_Tight_BI',package='ggchord',envir=e)
 s<-e$plasmid_example_pTRE_Tight_BI;f<-find_common_features(s)
 pos<-position_feature_stack(spacing=.085,base_position=position_plasmid())
 p<-ggchord(s[c('accver','length')],feature_data=f)+geom_seq()+
  geom_feature_plasmid(position=pos)+geom_feature_label_repel(position=pos,
    label_repeat='bracket',external=FALSE,label_fit='shrink',label_min_size=2.5)+
  coord_circular(rotation=48,direction='counterclockwise')+theme_ggchord_plasmid()
 mirror<-get_chord_layout(p)
 stopifnot(sum(mirror$gene_labels$repeat_count>1L)==1L,
   mirror$gene_labels$repeat_count[mirror$gene_labels$repeat_count>1L]==7L,
   length(unique(mirror$gene_polys$source_row))==nrow(f))
}
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(cases=results,runtime_sha256=as.list(audit$engine),input_sha256=as.list(audit$inputs)),
 'tools/long-term-audit/v0.27.0-compact-repeat-results.json',pretty=TRUE,auto_unbox=TRUE)
