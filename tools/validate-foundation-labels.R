# Local controls isolate content, size and the old count-based solver switch.
# No package tests, vignette knitting or release device matrix.
pkgload::load_all('.',quiet=TRUE,export_all=TRUE)
ns <- asNamespace('ggchord')
old_env <- new.env(parent=ns)
eval(parse(text=system2('git',c('show','f32e23a:R/layout-external-annotation.R'),stdout=TRUE)),old_env)
original <- get('ggchord_share_external_annotations',ns)
select_solver <- function(old=FALSE) assignInNamespace('ggchord_share_external_annotations',
  if(old) old_env$ggchord_share_external_annotations else original,ns='ggchord')
metrics <- function(p, l = NULL) {
  if (is.null(l)) l <- get_chord_layout(p,build=TRUE)
  unit <- l$render_context$units_per_inch %||% l$text_units_per_inch
  tasks <- ggchord_perimeter_tasks(l$layer_geometry,unit)
  overlaps <- 0L
  if(nrow(tasks)>1) for(i in seq_len(nrow(tasks)-1L)) for(j in seq.int(i+1L,nrow(tasks)))
    overlaps <- overlaps + as.integer(abs(tasks$text_x[i]-tasks$text_x[j]) < (tasks$w[i]+tasks$w[j])/2 &&
      abs(tasks$text_y[i]-tasks$text_y[j]) < (tasks$h[i]+tasks$h[j])/2)
  visibility <- ggchord_label_visibility(l)
  list(external_labels=nrow(tasks),overlapping_boxes=overlaps,
    requested=nrow(visibility),visible=sum(visibility$visible),
    hidden_reasons=if(nrow(visibility)) as.list(table(visibility$reason[!visibility$visible])) else list(),
    minimum_size=if(nrow(tasks)) min(tasks$size) else NULL)
}
synthetic <- function(n,size=2.5) {
  seqs<-data.frame(accver='synthetic',length=10000)
  # One fixed, explicitly synthetic 58-label population; no duplicated records.
  sites<-data.frame(accver='synthetic',position=round(seq(80,9800,length.out=58)),
    enzyme=paste0('SyntheticEnzyme',seq_len(58)))
  ggchord(seqs)+geom_seq()+geom_restriction_site(data=sites[seq_len(n),],label_size=size)+coord_circular()
}
if ("--flagship-sizes" %in% commandArgs(trailingOnly=TRUE)) {
  data(plasmid_example_pETDuet_1)
  source <- plasmid_example_pETDuet_1
  features <- find_common_features(source)
  primers <- find_primer_bindings(source,set="reference")
  sites <- filter_restriction_sites(find_restriction_sites(source),set="unique_6plus",parent_set="commercial_all")
  position <- position_feature_stack(spacing=.085,base_position=position_plasmid())
  p <- ggchord(source[c("accver","length")],feature_data=features)+geom_seq(seq_style="double")+
    geom_feature_plasmid(position=position)+geom_feature_label_repel(position=position,external=TRUE,
      label_fit="shrink",label_min_size=2.5,max_overlaps=0)+geom_primer(data=primers)+
    geom_primer_label_repel(data=primers,max_overlaps=0)+geom_restriction_site(data=sites)+
    geom_seq_center_label()+coord_circular()+theme_ggchord_plasmid()
  results <- list(input=list(features=nrow(features),primers=nrow(primers),sites=nrow(sites)),sizes=list())
  widths <- sub("^--width=","",grep("^--width=",commandArgs(trailingOnly=TRUE),value=TRUE))
  if (!length(widths)) widths <- c(10,16.5,8) else widths <- as.numeric(widths)
  for (width in widths) {
    cat("flagship size",width,"x 8\n");flush.console()
    grDevices::png(paste0("/private/tmp/ggchord-flagship-",width,"x8.png"),width=width,height=8,units="in",res=110)
    table <- ggplot2::ggplotGrob(p)
    grid::grid.newpage();grid::grid.draw(table)
    layout <- attr(table,"ggchord_render_layout")
    results$sizes[[paste0(width,"x8")]] <- c(metrics(p,layout),list(measurement=attr(table,"ggchord_render_context")))
    grDevices::dev.off()
  }
  jsonlite::write_json(results,"tools/long-term-audit/foundation-size-results-v0.17.0.json",pretty=TRUE,auto_unbox=TRUE)
  quit(save="no",status=0)
}
results <- list(controls=list(),catalog=list())
for(n in c(0,1,3,10,29,30,31,58)) {
  cat("label control",n,"\n"); flush.console()
  select_solver(TRUE);before<-metrics(synthetic(n))
  select_solver();after<-metrics(synthetic(n))
  stopifnot(before$external_labels==after$external_labels)
  results$controls[[paste0('n',n)]]<-list(before=before,after=after)
}
results$controls$size_2_6<-metrics(synthetic(30,2.6))
select_solver()
objects<-c('pSB1C3','pUC19','pBR322','pBluescript_II_SK_plus','pET_28a_plus',
 'pETDuet_1','pDONR221','pcDNA3_1_plus','pTRE_Tight_BI','pTRIPZ','pCAMBIA1300','pEarleyGate_201','pSpCas9_BB_2A_GFP_PX458')
for(id in objects) {
  cat('catalog',id,'\n');flush.console()
  env<-new.env();data(list=paste0('plasmid_example_',id),package='ggchord',envir=env)
  source<-env[[paste0('plasmid_example_',id)]]
  features<-find_common_features(source);sites<-find_restriction_sites(source)
  nearest<-if(nrow(sites)) unique(vapply(source$length*c(1,3,5)/6,function(x) which.min(abs(sites$position-x)),integer(1))) else integer()
  position<-position_feature_stack(spacing=.085,base_position=position_plasmid())
  p<-ggchord(source[c('accver','length')],feature_data=features)+geom_seq(seq_style='double')+
    geom_feature_plasmid(position=position)+geom_feature_label_repel(position=position,external=TRUE,
      label_fit='shrink',label_min_size=2.5)+geom_restriction_site(data=sites[nearest,],label_size=2.5)+
    coord_circular()+theme_ggchord_plasmid()
  results$catalog[[id]]<-metrics(p)
  if(id=='pETDuet_1') {
    primers<-find_primer_bindings(source,set='reference')
    flagship<-p+geom_primer(data=primers)+geom_primer_label_repel(data=primers)
    results$flagship<-metrics(flagship)
    ggplot2::ggsave('/private/tmp/ggchord-foundation-plasmid.png',flagship,width=10,height=8,dpi=120)
  }
}
jsonlite::write_json(results,'tools/long-term-audit/foundation-label-results-v0.17.0.json',pretty=TRUE,auto_unbox=TRUE)
cat('Local label controls recorded; collision counts are reported, not assumed zero.\n')
