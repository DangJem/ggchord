# Explicit reference-profile comparison, independent of the sparse gallery.
source('tools/audit-source-snapshot.R')
source('tools/long-term-audit/feature-association-metrics.R')
audit <- audit_source_snapshot(additional=c('tools/long-term-audit/render-reference-catalogue.R',
 'tools/long-term-audit/feature-association-metrics.R'))
pkgload::load_all(quiet=TRUE)
ids<-c('pUC19','pBR322','pSB1C3','pBluescript_II_SK_plus','pET_28a_plus',
 'pETDuet_1','pcDNA3_1_plus','pDONR221','pCAMBIA1300','pTRE_Tight_BI',
 'pEarleyGate_201','pTRIPZ','pSpCas9_BB_2A_GFP_PX458')
# Optional subset is useful for revisiting one observed failure.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)) {stopifnot(all(args %in% ids));ids<-args}
output_prefix <- Sys.getenv('GGCHORD_CATALOGUE_PREFIX', 'v0.27.0-reference-catalogue')
stopifnot(grepl('^[A-Za-z0-9._-]+$',output_prefix))
results<-list()
for(id in ids) {
 cat('Reference catalogue:',id,'\n');flush.console()
 e<-new.env();data(list=paste0('plasmid_example_',id),package='ggchord',envir=e)
 sequence<-e[[paste0('plasmid_example_',id)]]
 f<-find_common_features(sequence)
 all_sites<-find_restriction_sites(sequence)
 profile<-attr(all_sites,'reference_enzyme_profiles')
 sites<-filter_restriction_sites(all_sites,set='reference',parent_set='commercial_nonredundant')
 policy<-'saved reference profile; current commercial nonredundant representatives'
 if(nrow(profile) && profile$profile_type=='none') {
  # The old PNG visibly contains cutters although the .dna saves None.
  # This is an explicit comparison policy, not a recovered export setting.
  sites<-filter_restriction_sites(all_sites,set='unique_6plus',parent_set='commercial_nonredundant')
  policy<-'explicit unique 6+ commercial nonredundant comparison; saved None disagrees with PNG'
 }
 primers<-find_primer_bindings(sequence,set='reference')
 position<-position_feature_stack(spacing=.085,base_position=position_plasmid())
 p<-ggchord(sequence[c('accver','length')],feature_data=f)+geom_seq(seq_style='double')+
  geom_feature_plasmid(position=position)+
  geom_feature_label_repel(position=position,label_fit='shrink',label_min_size=2.5)+
  geom_seq_center_label()+coord_circular()+theme_ggchord_plasmid()
 if(nrow(sites))p<-p+geom_restriction_site(data=sites,label_size=2.5)
 if(nrow(primers))p<-p+geom_primer(data=primers)+geom_primer_label_repel(data=primers,size=2.5)
 path<-paste0('tools/long-term-audit/',output_prefix,'-',id,'.png')
 ggchord_diagnostic_device(path,'png',10,8,110)
 warnings<-character()
 table<-withCallingHandlers(ggplot2::ggplotGrob(p),warning=function(w){
  warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 grid::grid.newpage();grid::grid.draw(table)
 l<-attr(table,'ggchord_render_layout')
 axis <- ggchord_axis_geometry(l)
 axis_rows <- axis$.component %in% 'text' & !is.na(axis$.axis_label)
 axis_labels <- do.call(rbind,lapply(unique(axis$.axis_label[axis_rows]),function(value) {
   at <- axis_rows & axis$.axis_label %in% value
   data.frame(label=value,visible=any(!is.na(axis$label[at]) & nzchar(trimws(axis$label[at]))),
     endpoint=any(axis$is_origin[at] %in% TRUE | axis$is_endpoint[at] %in% TRUE |
       axis$.axis_seam_combined[at] %in% TRUE),
     status=paste(unique(axis$.axis_clearance_status[at]),collapse=';'))
 }))
 stopifnot(is.data.frame(axis_labels),any(axis_labels$endpoint),
   all(axis_labels$visible[axis_labels$endpoint]))
 results[[id]]<-list(file=path,site_policy=policy,selected_sites=nrow(sites),
  selected_positions=length(unique(sites$position)),source_features=nrow(f),source_primers=nrow(primers),
  visibility=ggchord_label_visibility(l),converged=attr(table,'ggchord_render_context')$converged,
  warnings=warnings,axis_labels=axis_labels,quality=association_text_quality(table)$summary,
  association=association_metrics(l)$summary,
  feature_exterior_from_inside=sum(vapply(l$layer_geometry,function(owner) {
    d<-owner$gene_label_repel
    if (!is.data.frame(d) || !'.feature_internal_exit' %in% names(d))return(0L)
    sum(d$.component %in% 'text' & d$.feature_internal_exit %in% TRUE)
  },integer(1))),
  built_internal_exit_segments=sum(vapply(attr(table,'ggchord_render_built')$data,
    function(d)if (all(c('.component','.feature_internal_exit') %in% names(d)))
      sum(d$.component %in% 'segment' & d$.feature_internal_exit %in% TRUE) else 0L,integer(1))))
 if (results[[id]]$feature_exterior_from_inside>0L)
   stopifnot(results[[id]]$built_internal_exit_segments>0L)
 saveRDS(l,paste0('/tmp/',output_prefix,'-',id,'.rds'))
 grDevices::dev.off()
 jsonlite::write_json(list(cases=results,runtime_sha256=as.list(audit$engine),
   input_sha256=as.list(audit$inputs),scope='10x8in/110dpi source-profile comparison; not exact original PNG export settings'),
  paste0('tools/long-term-audit/',output_prefix,'-results.json'),pretty=TRUE,auto_unbox=TRUE)
}
audit_assert_source_snapshot(audit)
stopifnot(all(vapply(results,function(x) x$association$anchors_off_own_feature==0L,logical(1))))
cat('Reference catalogue completed\n')
