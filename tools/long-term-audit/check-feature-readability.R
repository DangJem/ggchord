# Local geometry and drawing controls; not testthat or a release device matrix.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional='tools/long-term-audit/check-feature-readability.R')
pkgload::load_all(quiet=TRUE)
s <- plasmid_example_pBluescript_II_SK_plus
f <- find_common_features(s)
pos <- position_feature_stack(spacing=.085,base_position=position_plasmid())
base <- ggchord(s[c('accver','length')],feature_data=f)+geom_seq(seq_style='double')+
 geom_feature_plasmid(position=pos)+coord_circular()+theme_ggchord_plasmid()
p <- base+geom_feature_label_repel(position=pos,label_fit='shrink',label_min_size=2.5)
l <- get_chord_layout(p)
# Every visible boundary fragment lies inside its own drawable contour.
polygon <- l$layer_geometry[['layer-0002']]$gene_poly
marks <- split(polygon[polygon$.component %in% 'boundary',],
 polygon$group[polygon$.component %in% 'boundary'])
stopifnot(length(marks)>0L)
for(mark in marks) {
 own <- polygon[polygon$.component %in% 'polygon' & polygon$source_row==mark$source_row[1L],]
 contours <- split(own,own$group)
 for(t in c(.1,.5,.9)) {
  point <- unlist(mark[1L,c('x','y')]+t*(mark[2L,c('x','y')]-mark[1L,c('x','y')]))
  stopifnot(any(vapply(contours,function(poly) ggchord_point_in_polygon(point[1],point[2],poly),logical(1))))
 }
}
# Curved glyph baselines use the radius selected for their complete label.
angles <- radii <- numeric()
for(owner in l$layer_geometry) {
 d<-owner$gene_label_repel
 if(!is.data.frame(d))next
 glyph<-d[d$.component %in% 'arc_text',]
 text<-d[d$.component %in% 'text',]
 if(!nrow(glyph))next
 tangent<-atan2(glyph$y,glyph$x)*180/pi+90
 angles<-c(angles,abs(((glyph$angle-tangent+90)%%180)-90))
 expected<-sqrt(text$x^2+text$y^2)[match(glyph$source_row,text$source_row)]
 radii<-c(radii,abs(sqrt(glyph$x^2+glyph$y^2)-expected))
}
stopifnot(length(angles)>0,max(angles)<1,max(radii)<.05*l$text_units_per_inch/25.4)
# Wrapped internal text uses one concentric baseline for each line.
wrapped<-data.frame(text="first line\nsecond line",text_x=.4,text_y=.4,text_angle=-45,
 size=2.5,family="",fontface=1,lineheight=1.2,group=1L,hjust=.5,vjust=.5)
arc<-ggchord_arc_text_layout(wrapped)
stopifnot(setequal(unique(arc$glyphs$.arc_line),1:2))
for(line in 1:2) {
 glyph<-arc$glyphs[arc$glyphs$.arc_line==line,]
 radius<-sqrt(glyph$text_x^2+glyph$text_y^2)
 error<-abs(((glyph$text_angle-atan2(glyph$text_y,glyph$text_x)*180/pi-90+90)%%180)-90)
 stopifnot(diff(range(radius))<.05*.35/25.4,max(error)<1)
}
# The draw component respects automatic styling, theme and explicit params.
segment<-data.frame(x=.1,y=.1,xend=.2,yend=.2,.component='segment',
 feature_label_mode='external',.feature_internal_exit=TRUE,colour='black',alpha=1,linewidth=.25,linetype='solid',
 label=NA_character_,size=2.5,angle=0,hjust=.5,vjust=.5,family='',fontface=1,lineheight=1.2)
built<-ggplot2::ggplot_build(ggplot2::ggplot()+ggplot2::geom_blank()+ggplot2::coord_fixed())
draw<-function(style=list(colour='auto',linetype='auto'),params=list(.automatic_feature_segment=TRUE))
 GeomChordGeneLabelRepel$draw_panel(segment,built$layout$panel_params[[1]],built$plot$coordinates,
 feature_segment_style=style,segment_params=params)$children[[1]]$gp
standard<-draw();explicit<-draw(list(colour='red',linetype='dotdash'))
segment$.feature_internal_exit <- FALSE
inside<-draw()
segment$.feature_internal_exit <- TRUE
custom<-draw(params=list(.automatic_feature_segment=FALSE,colour='blue',linetype='dotted'))
stopifnot(standard$col=='#8A8A8AFF',standard$lty=='dashed',
 inside$col=='#000000FF',inside$lty=='solid',
 explicit$col=='#FF0000FF',explicit$lty=='dotdash',custom$col=='#0000FFFF',custom$lty=='dotted')
# Isolated nearby names omit a spur; displaced names still need association.
row<-data.frame(feature_label_mode='adjacent',repeat_count=1L,x=0,y=.85,
 anchor_x=0,anchor_y=.86,text_x=0,text_y=.85,.attachment_id='own')
box<-data.frame(cx=0,cy=.85,ow=.08,oh=.02,angle=0,xmin=-.04,xmax=.04,ymin=.84,ymax=.86)
stopifnot(ggchord_feature_short_attachment(row,box,list(),.35))
row$x<-.1;row$text_x<-.1
stopifnot(!ggchord_feature_short_attachment(row,box,list(),.35))
# Actual-size wrapped control; source feature records and coordinates unchanged.
wrapped_plot <- base+geom_feature_label_repel(position=pos,label_wrap=10,
 label_fit='wrap',label_max_lines=3)
path <- 'tools/long-term-audit/v0.27.0-wrapped-feature-labels.png'
ggchord_diagnostic_device(path,'png',10,8,110)
wrapped_table <- ggplot2::ggplotGrob(wrapped_plot)
grid::grid.newpage();grid::grid.draw(wrapped_table)
wrapped_layout <- attr(wrapped_table,'ggchord_render_layout')
wrapped_text <- wrapped_layout$gene_labels
wrapped_internal <- wrapped_text$feature_label_mode %in% c('inside','adjacent') &
 grepl('\n',wrapped_text$text,fixed=TRUE)
stopifnot(any(wrapped_internal))
wrapped_geometry <- wrapped_layout$layer_geometry[['layer-0003']]$gene_label_repel
wrapped_sources <- wrapped_text$source_row[wrapped_internal]
stopifnot(all(vapply(wrapped_sources,function(source) {
 d<-wrapped_geometry[wrapped_geometry$.component %in% 'arc_text' &
   wrapped_geometry$source_row==source,,drop=FALSE]
 nrow(d)>0L && length(unique(d$.arc_line))>1L
},logical(1))))
grDevices::dev.off()
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(boundary_fragments=length(marks),
 max_glyph_tangent_error_deg=max(angles),max_glyph_radius_error_native=max(radii),
 automatic_inner_to_exterior_grey_dashed=TRUE,explicit_and_theme_style_priority=TRUE,
 isolated_spur_suppression=TRUE,wrapped_internal_labels=sum(wrapped_internal),
 wrapped_converged=attr(wrapped_table,'ggchord_render_context')$converged,
 runtime_sha256=as.list(audit$engine)),
 'tools/long-term-audit/v0.27.0-feature-readability-results.json',pretty=TRUE,auto_unbox=TRUE)
cat('Boundary, circular glyph, short attachment and drawing style controls passed\n')
