# Synthetic style/layout fixture; no biological measurements are implied.
pkgload::load_all('.',quiet=TRUE,export_all=TRUE)
s<-data.frame(accver=paste0('sequence-',seq_len(12)),length=1000)
f<-data.frame(accver=s$accver,start=100,end=400,strand='+',type=rep(paste0('class-',1:6),2))
p<-ggchord(s,feature_data=f)+geom_seq()+geom_feature(aes(feature_fill=type),shape='block')+
  coord_collinear()+theme(legend.position='right',legend.text=element_text(size=9),legend.key.height=grid::unit(6,'mm'),axis.text=element_blank())
d<-diagnose_ggchord(p,width=7,height=5,device='png')
stopifnot(identical(d$measurement$status,'measured'),setequal(d$guide_measurements$n_items,c(12L,6L)))
jsonlite::write_json(list(scope='two categorical roles; shared side budget and final guide/panel measurement',
  guides=d$guides,guide_measurements=d$guide_measurements,checks=d$checks,measurement=d$measurement),
  'tools/long-term-audit/foundation-guide-results-v0.17.0.json',pretty=TRUE,auto_unbox=TRUE)
cat('Guide measurements recorded.\n')
