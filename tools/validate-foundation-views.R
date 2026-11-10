# Direct Foundation behaviour checks, not package or release acceptance.
pkgload::load_all('.',quiet=TRUE,export_all=TRUE)
results <- list()
s <- data.frame(accver=c('A','B','C'),length=c(1000,600,400))
h <- data.frame(qaccver=c('A','B'),saccver=c('B','A'),qstart=c(100,200),qend=c(400,400),
  sstart=c(400,500),send=c(100,800),pident=c(95,80),source_hsp_id=c('h1','h2'))
f <- data.frame(accver=c('A','B','C'),start=50,end=100,type=c('CDS','misc_feature','tRNA'),strand=c('+','.','-'))
p <- suppressWarnings(ggchord(s,h,feature_data=f)) +
  geom_seq(seq_order=c(A_high='A',B='B',A_low='A'),seq_orientation=c(A=1,B=1,A_low=-1)) +
  geom_feature(shape='block') +
  geom_link_ribbon(data=h[1,],between=c('A_high','B')) +
  geom_link_ribbon(data=h[2,],between=c(2,3))
for (coord in list(coord_chord(),coord_collinear(),coord_collinear(seq_axis='vertical'))) {
  view <- p+coord
  e <- export_ggchord_layout(view,include=c('seq','ribbon','feature'),original_data=TRUE)
  stopifnot(nrow(e$display_instances)==3L,nrow(e$instance_routing)==2L,
    identical(e$instance_routing$source_hsp_id,c('h1','h2')),
    setequal(e$instance_routing$qaccver,c('A','B')))
  ggplot2::ggplotGrob(view)
}
results$instances <- 'A–B–A routes and source identities persist across Chord and both Collinear axes; omitted source annotation is permitted'
for (coord in list(coord_chord(),coord_collinear(),coord_collinear(seq_axis='vertical'))) {
  q <- ggchord(s)+geom_seq(length_mode='equal')+coord
  lay <- get_chord_layout(q,build=TRUE)
  stopifnot(identical(unname(export_ggchord_layout(q,include="seq")$sequence_display$source_length),s$length))
  if (isTRUE(coord$ggchord_collinear)) {
    m <- map_ggchord_position(q,s$accver,s$length)
    stopifnot(all(abs(map_ggchord_position(q,s$accver,if(identical(coord$collinear_seq_axis,'vertical'))m$y else m$x,inverse=TRUE)$position-s$length)<1e-7))
  }
  ggplot2::ggplotGrob(q)
}
results$equal <- 'source lengths retained; forward and inverse positions round-trip in both axes'
r <- rbind(h,transform(h[1,],qaccver='A',saccver='A',source_hsp_id='self'))
q <- suppressWarnings(ggchord(s,r))+geom_seq()+geom_link_ribbon()+coord_sankey()
e <- export_ggchord_layout(q,include=c('seq','ribbon'))
stopifnot(nrow(e$display_instances)==6,nrow(e$instance_routing)==2,
  all(grepl('^left:',e$instance_routing$q_instance)),all(grepl('^right:',e$instance_routing$s_instance)))
m <- map_ggchord_position(q,'right:B',300)
stopifnot(abs(map_ggchord_position(q,'right:B',m$y,inverse=TRUE)$position-300)<1e-7)
invisible(ggplot2::ggplotGrob(q))
results$sankey <- 'real query-left/subject-right endpoints; reciprocals retained, same entity excluded; inverse mapping round-trip'
q <- ggchord(s)+geom_seq()+coord_collinear()
lay <- get_chord_layout(q,build=TRUE)
for (id in s$accver) stopifnot(all(c('1',as.character(s$length[s$accver==id])) %in% lay$axis_ticks$label[lay$axis_ticks$accver==id]))
q <- q+scale_seq_position_continuous(breaks=c(100,200),minor_breaks=NULL)
stopifnot(setequal(na.omit(get_chord_layout(q,build=TRUE)$axis_ticks$label),c('100','200')))
results$axis <- 'default real endpoints and explicit-break precedence'
loci <- data.frame(accver='A',start=25,end=200,locus_id='a')
stopifnot(nrow(focus_ggchord_data(s,feature_data=f,loci=loci)$feature_data)==1,
  nrow(segment_ggchord_data(s,feature_data=f,windows=loci)$feature_data)==1,
  nrow(clean_ggchord_data(s,feature_data=f)$feature_data)==3,
  nrow(prepare_ggchord_data(s,feature_data=f,steps=list(clean=list()))$feature_data)==3)
results$workflow <- 'unified mixed/unknown-strand features survive focus, segment, clean, prepare'
jsonlite::write_json(results,'tools/long-term-audit/foundation-views-results-v0.17.0.json',pretty=TRUE,auto_unbox=TRUE)
cat('Foundation view checks passed.\n')
