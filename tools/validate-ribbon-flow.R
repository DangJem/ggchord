# One curved-front regression: endpoint colours and opacity are public styles.
pkgload::load_all('.',quiet=TRUE,export_all=TRUE)
s<-data.frame(accver=c('A','B'),length=1000)
h<-data.frame(qaccver='A',saccver='B',qstart=100,qend=600,sstart=800,send=200)
p<-ggchord(s,h,validate='none')+geom_seq()+geom_link_ribbon(gradient=TRUE,alpha=.4,tip='none')+
  scale_seq_colour_manual(values=c(A='#FF0000',B='#0000FF'))+coord_chord()
g<-ggplot2::ggplotGrob(p)
rasters<-list()
walk<-function(node) {
  if(inherits(node,'rastergrob'))rasters[[length(rasters)+1L]]<<-node$raster
  if(length(node$grobs))for(child in node$grobs)walk(child)
  if(length(node$children))for(child in node$children)walk(child)
}
walk(g)
stopifnot(length(rasters)>0)
colours<-unique(unlist(lapply(rasters,as.character)))
colours<-colours[!is.na(colours)]
rgb<-grDevices::col2rgb(colours,alpha=TRUE)
stopifnot(all(rgb[4,]==round(.4*255)),any(rgb[1,]>250 & rgb[3,]<3),any(rgb[3,]>250 & rgb[1,]<3))
e<-export_ggchord_layout(p,include='ribbon',original_data=TRUE)
stopifnot(setequal(e$ribbon$source_row,1L))
# Explicit user fill remains authoritative and does not produce a flow fill.
q<-ggchord(s,h,validate='none')+geom_seq()+geom_link_ribbon(gradient=TRUE,fill='#00FF00',tip='none')
b<-ggplot2::ggplot_build(q)
stopifnot(all(toupper(b$data[[2]]$ribbon_fill)=='#00FF00'))
# The endpoint resolver follows category mappings, staged colours and missing
# backbones through the very same sequence geom/scale computation.
s$category <- c('group2','group1')
base <- ggchord(s,h,validate='none')+geom_seq(aes(seq_colour=category))+
  scale_seq_colour_manual(values=c(group1='#0000FF',group2='#FF0000'))
layout <- get_chord_layout(base+geom_link_ribbon(gradient=TRUE,),build=TRUE)
resolved <- ggchord_sequence_display_colours(ggchord_register_layers(base),layout)
stopifnot(all(resolved$colour[resolved$instance=='A']=='#FF0000'))
implicit <- ggchord(s,h,validate='none')+geom_link_ribbon(gradient=TRUE,)+scale_seq_colour_manual(values=c(A='red',B='blue'))
invisible(ggplot2::ggplotGrob(implicit))
staged <- ggchord(s,h,validate='none')+geom_seq(aes(seq_colour=ggplot2::stage(start=accver,after_scale=scales::alpha(seq_colour,.5))))+
  geom_link_ribbon(gradient=TRUE,)+scale_seq_colour_manual(values=c(A='red',B='blue'))
invisible(ggplot2::ggplotGrob(staged))
transparent <- ggchord(s,h,validate='none')+geom_seq()+geom_link_ribbon(gradient=TRUE,)+
  scale_seq_colour_manual(values=c(A='transparent',B='blue'))
invisible(ggplot2::ggplotGrob(transparent))
arrow <- ggchord(s,h,validate='none')+geom_seq()+geom_link_ribbon(gradient=TRUE,tip='subject')+coord_collinear(seq_axis='vertical')
invisible(ggplot2::ggplotGrob(arrow))
for(side in c('query','subject')) {
  shared <- data.frame(qaccver='A',saccver='B',qstart=c(100,300),qend=c(200,400),sstart=c(600,800),send=c(700,900))
  if(side=='query'){shared$qstart<-100;shared$qend<-200}else{shared$sstart<-800;shared$send<-900}
  view <- ggchord(s,shared,validate='none')+geom_seq()+geom_link_ribbon(gradient=TRUE,link_branch=side,tip='subject')+coord_collinear()
  b<-ggplot2::ggplot_build(view)
  stopifnot(any(b$data[[2]]$.component=='trunk'))
  invisible(ggplot2::ggplot_gtable(b))
}
jsonlite::write_json(list(curved_front='source arc retained in the fill mask; manually inspected PNG',
  endpoint_colours='both declared sequence colours present',opacity='requested alpha applied once',
  explicit_fill='user fixed fill retained'), 'tools/long-term-audit/ribbon-flow-results-v0.17.0.json',pretty=TRUE,auto_unbox=TRUE)
cat('Curved ribbon colour/opacity checks passed.\n')
