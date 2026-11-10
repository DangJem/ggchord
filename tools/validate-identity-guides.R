# Confirmed regression: generic 5-mm legend keys collapsed the whole ruler.
pkgload::load_all('.',quiet=TRUE,export_all=TRUE)
s<-data.frame(accver=c('A','B'),length=1000)
h<-data.frame(qaccver='A',saccver='B',qstart=seq(10,60,10),qend=seq(100,150,10),
 sstart=seq(200,250,10),send=seq(300,350,10),pident=c(0,50,80,90,95,100))
p<-ggchord(s,h)+geom_seq()+geom_link_ribbon()+coord_collinear()+
 theme_ggchord_publication()+theme(legend.text=element_text(size=11),
  legend.key.height=unit(5,'mm'),legend.key.width=unit(9,'mm'))
results<-list()
for(position in c('right','bottom')) {
 d<-diagnose_ggchord(p+theme(legend.position=position),width=10,height=7,device='png')
 ruler<-d$guide_measurements[d$guide_measurements$role=='ribbon',,drop=FALSE]
 stopifnot(nrow(ruler)==1L)
 extent<-if(position=='right')ruler$requested_key_height_mm else ruler$requested_key_width_mm
 text<-if(position=='right')ruler$max_label_height_mm else ruler$max_label_width_mm
 stopifnot(extent*.05>text)
 results[[position]]<-list(guides=d$guides,measurements=ruler,
  closest_boundary_distance_mm=extent*.05,label_extent_mm=text)
}
manual <- ggchord(s,h)+geom_seq(show.legend=FALSE)+geom_link_ribbon()+
  coord_collinear()+theme_ggchord_publication()+theme(legend.position="bottom",
    legend.text=element_text(size=11),legend.key.height=unit(5,"mm"),legend.key.width=unit(9,"mm"))+
  guides(ribbon_fill=guide_ggchord_colourbar())
md <- diagnose_ggchord(manual,width=10,height=7,device="png")
box <- md$guides[md$guides$guide=="guide-box-bottom",,drop=FALSE]
stopifnot(nrow(box)==1L,box$width_mm>2*box$height_mm)
results$wrapper_bottom <- list(guides=md$guides,measurements=md$guide_measurements,
  evidence="late bottom placement retains a horizontal ruler")
b<-ggplot2::ggplot_build(p)
stopifnot(!ggchord_ribbon_uses_flow(p$layers[[2]],p),
 length(unique(b$data[[2]]$ribbon_fill))>1L)
jsonlite::write_json(list(scope='synthetic 0/50/80/90/95/100 Identity boundaries with small generic keys',
 default_gradient=FALSE,results=results),
 'tools/long-term-audit/identity-guide-results-v0.17.0.json',pretty=TRUE,auto_unbox=TRUE)
cat('Identity labels retain physical separation in both directions\n')
