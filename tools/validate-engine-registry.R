# Selected shared-annotation and physical-band link risks, not full visual acceptance.
pkgload::load_all(quiet=TRUE)
local({
  file<-tempfile(fileext='.pdf');grDevices::pdf(file,width=8,height=6)
  on.exit({grDevices::dev.off();unlink(file)},add=TRUE)
  results<-list()
  record<-function(case,detail) results[[length(results)+1L]]<<-
    data.frame(case=case,result='PASS',detail=detail)
  s<-data.frame(accver='A',length=1000)
  a<-data.frame(accver='A',start=100,end=280,strand='+',type='CDS',label='first')
  b<-data.frame(accver='A',start=600,end=820,strand='-',type='misc',label='second')
  position<-position_plasmid(-.1)
  entities<-ggchord(s)+geom_seq()+geom_feature(data=a,position=position)+
    geom_feature(data=b,position=position)+
    geom_feature_label_repel(data=a,position=position)+
    geom_feature_label_repel(data=b,position=position)
  circular<-export_ggchord_layout(entities+coord_circular(direction='counterclockwise',gap=15),
    include=c('feature','labels'),original_data=TRUE)
  registry<-circular$annotation_registry
  labels<-circular$labels
  text<-labels[labels$.component %in% c('text','arc_text'),]
  stopifnot(length(unique(circular$feature$layer_id))==2L,
    length(unique(text$layer_id))==2L,all(text$source_row==1L),
    all(c('first','second') %in% unique(text$label)),nrow(registry)>0,
    all(c('track','slot','owner','component','bbox_xmin','leader_corridor') %in% names(registry)),
    all(registry$owner %in% c('.coord',unique(c(circular$feature$layer_id,labels$layer_id)))))
  record('F21_circular_ownership','Same source-row numbers across independent feature/label owners retain text and circular resource provenance')
  linear<-export_ggchord_layout(entities+coord_collinear(),include=c('feature','labels'))
  d<-linear$label_diagnostics
  stopifnot(nrow(d)==2L,length(unique(d$layer_id))==2L,all(d$source_row==1L),
    all(d$status %in% c('placed','hidden')),all(!is.na(d$reason[d$status=='hidden'])),
    all(d$leader_crossings[d$status=='placed']==0))
  record('F21_linear_ownership','Linear shared solver retains each owner and reports placement/hiding and leader evidence')
  pair<-data.frame(accver=c('A','B'),length=1000)
  genes<-data.frame(accver='A',start=100,end=300,strand='+',anno='gene')
  r<-data.frame(qaccver='A',saccver='B',qstart=100,qend=250,sstart=400,send=600)
  bands<-data.frame(accver='A',start=100,end=300,value=5)
  base<-ggchord(pair,r)+geom_seq(seq_gap=.5)+geom_gene(data=genes)+
    geom_track_bar(aes(track_value=value),data=bands,
      position=position_track('inner',side='inside',width=.03,gap=.01))+
    coord_collinear()
  automatic<-export_ggchord_layout(base+geom_link_ribbon(avoid='uniform'),include='ribbon')
  explicit<-export_ggchord_layout(base+geom_link_ribbon(avoid='uniform',link_gap=.02),include='ribbon')
  stopifnot(nrow(automatic$ribbon)>0,nrow(explicit$ribbon)>0,
    all(automatic$ribbon$q_gap>explicit$ribbon$q_gap),
    all(explicit$ribbon$q_gap==.02),
    setequal(automatic$ribbon$source_row,1L))
  record('F18_F21_link_band_avoidance','Requested automatic avoidance consumes registered bands; explicit endpoint gap remains authoritative')
  # Branching consumes mapped visible styles and must retain source membership.
  three<-data.frame(accver=c('A','B','C'),length=1000)
  blocks<-data.frame(qaccver='A',saccver=c('B','C'),qstart=100,qend=300,sstart=400,send=600)
  shared<-ggchord(three,blocks)+geom_seq()+geom_link_ribbon(fill='orange',
    link_branch='query',link_ctrl_point=c(4,4))
  invisible(ggplot2::ggplot_build(shared))
  branch<-export_ggchord_layout(shared,include='ribbon',original_data=TRUE)
  stopifnot(any(branch$ribbon$.component=='trunk'),
    all(c(1L,2L) %in% unlist(branch$ribbon$source_rows)),
    all(branch$ribbon$layer_id %in% names(branch$original_data)))
  links<-data.frame(qaccver='A',saccver=c('B','C'),qpos=200,spos=500)
  styled<-ggplot2::ggplot_build(ggchord(three)+geom_seq()+geom_link_line(
    aes(link_alpha=after_scale(seq_along(x)/length(x))),data=links,link_branch='query'))$data[[2]]
  stopifnot(!any(styled$.component=='trunk'),all(is.finite(styled$link_alpha)))
  record('F05_F22_branch_membership','Post-scale ribbon trunks retain source sets; varying staged path styles stay separate')
  out<-do.call(rbind,results)
  write.csv(out,'design-roadmap/audit-generated/engine-registry-validation-v0.16.3.csv',row.names=FALSE)
  print(out,row.names=FALSE)
})
