# One circular-axis glyph regression; no package suite or device matrix.
pkgload::load_all('.',quiet=TRUE,export_all=TRUE)
p<-ggchord(data.frame(accver='synthetic-axis',length=10000))+geom_seq(seq_style='double')+
  coord_circular()+theme_ggchord_plasmid()
g<-ggplot2::ggplotGrob(p)
l<-attr(g,'ggchord_render_layout');if(is.null(l))l<-get_chord_layout(p,build=TRUE)
ticks<-l$axis_ticks;ticks<-ticks[!is.na(ticks$label),,drop=FALSE]
glyphs<-ggchord_axis_label_geometry(l,ticks)
stopifnot(nrow(glyphs)>nrow(ticks))
for(label in unique(glyphs$.axis_label)) {
  rows<-which(glyphs$.axis_label==label)
  stopifnot(paste0(glyphs$label[rows],collapse='')==label)
  if(length(rows)>1L)stopifnot(diff(range(glyphs$label_angle[rows]))>0)
}
horizontal<-p+theme_ggchord(axis.text.orientation='horizontal')
h<-get_chord_layout(horizontal,build=TRUE)
hs<-h$axis_ticks;hs<-hs[!is.na(hs$label),,drop=FALSE]
stopifnot(nrow(ggchord_axis_label_geometry(h,hs))==nrow(hs))
jsonlite::write_json(list(scope='synthetic circular axis and explicit horizontal override',
  source_labels=as.character(ticks$label),glyph_count=nrow(glyphs),
  label_text_preserved=TRUE,tangents_follow_arc=TRUE,horizontal_stays_straight=TRUE),
  'tools/long-term-audit/plasmid-axis-results-v0.17.0.json',pretty=TRUE,auto_unbox=TRUE)
cat('Circular axis glyph and horizontal override checks passed\n')
