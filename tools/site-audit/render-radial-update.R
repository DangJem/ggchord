source('tools/audit-source-snapshot.R')
snapshot <- audit_source_snapshot(additional=c('tools/site-flagships.R','tools/home-previews.R','tools/documentation-examples.R','tools/site-audit/render-radial-update.R','website/home-previews.json','tools/site-course-examples.R','website/course-figures.json'))
pkgload::load_all(quiet=TRUE)
library(ggplot2)
library(dplyr)
library(grid)
blocks <- function(file) {
 x <- readLines(file,warn=FALSE); m <- grep('^# ---- [a-z]+ ----$',x)
 setNames(lapply(seq_len(length(m)-1),function(i) paste(x[seq.int(m[i]+1,m[i+1]-1)],collapse='\n')),sub('^# ---- ([a-z]+) ----$','\\1',x[m[-length(m)]]))
}
b <- c(blocks('tools/documentation-examples.R'),blocks('tools/site-course-examples.R'),blocks('tools/site-flagships.R'))
for(k in c('sequencepalette','docstyle','showcasedata','chordshowcasegenes','showcasechord','showcasecollinear','fusiondata','flagshipdata','flagshipmulti','fusionoverview','flagshipdot','flagshipfusion')) eval(parse(text=b[[k]]))
source('tools/home-previews.R')
layout <- get_chord_layout(flag_radial)
ends <- export_ggchord_layout(flag_radial)$radial_endpoints
stopifnot(nrow(ends)==12L, sum(!ends$adjacent)==4L,
  all(radial_orientation[c("MT108731.1","OR222515.1")]==1))
input <- Filter(function(x) !is.null(x$ribbon),layout$layer_inputs)[[1]]$ribbon
fields <- c("qaccver","saccver","qstart","qend","sstart","send","pident")
stopifnot(identical(as.data.frame(input[fields]),as.data.frame(show_links[fields])),nrow(chord_genes)==20L)
polys <- layout$ribbon_polys
frames <- layout$radial_frames
# Intersections of each nonadjacent ribbon polygon with other full backbone rays.
clearances <- numeric()
for(i in ends$source_row[!ends$adjacent]) {
  e <- ends[ends$source_row==i, ]
  xy <- as.matrix(polys[polys$source_row==i,c("x","y")]);xy <- rbind(xy,xy[1,])
  for(j in which(!frames$key %in% c(e$qaccver,e$saccver))) {
    u <- c(frames$tx[j],frames$ty[j]);v <- c(-u[2],u[1])
    r <- as.vector(xy%*%u);z <- as.vector(xy%*%v)
    k <- which(head(z,-1)*tail(z,-1)<0)
    if(!length(k)) next
    crossing <- r[k]+(r[k+1]-r[k])*(-z[k])/(z[k+1]-z[k])
    crossing <- crossing[crossing>0]
    tip <- frames$inner[j]+layout$sequence_layout$xmax[j]-layout$sequence_layout$xmin[j]
    if(length(crossing)) clearances <- c(clearances,crossing-tip)
  }
}
stopifnot(length(clearances)>0,all(clearances>.10))
built <- ggplot2::ggplot_build(flag_radial)$data[[1]]
stopifnot(all(built$linewidth[built$linetype=="dashed"]==.4),
  setequal(unique(built$linetype),c("solid","dashed")))
fixed <- flag_radial
fixed$layers[[1]] <- geom_link_ribbon(adjacent_only=FALSE,linetype="dotted",nonadjacent_linewidth=0)
stopifnot(all(ggplot2::ggplot_build(fixed)$data[[1]]$linetype=="dotted"))
fusion_export <- export_ggchord_fusion(flag_fusion, width=16, height=12)
stopifnot(nrow(flagship_hsp)==32L,nrow(flagship_features)==266L,
  nrow(fusion_export$connections)==3L,all(fusion_export$connections$retained),
  nrow(fusion_display_hsp)==6L,nrow(fusion_display_features)==12L,
  nrow(fusion_detail_features)<=6L,
  fusion_export$metadata$guide_layout$collected$shared==3L,
  all(vapply(fusion_export$diagnostics,function(d) length(d$guide_sizes_mm)==0L,logical(1))),
  all(fusion_export$source_layouts$center$sequence_reference$orientation==1),
  setequal(paste(fusion_export$connections$qpanel,fusion_export$connections$spanel),
    c("chord linear","chord radial","radial center")),
  nrow(fusion_export$panels)==5L,nrow(fusion_export$focus_regions)==2L,
  all(fusion_export$focus_regions$length_ratio>1),
  all(vapply(fusion_export$diagnostics,function(d) isTRUE(d$render_context$converged),logical(1))))
# Detail selection comes from the same full-source records shown in context.
stopifnot(nrow(fusion_detail_source_features)==4L,
  all(fusion_detail_features$feature_id %in% fusion_context_features$feature_id),
  identical(fusion_detail_source_features,
    flagship_features |> semi_join(fusion_detail_features,by=c("accver","feature_id"))))
in_window <- fusion_context_features |>
  inner_join(flagship_windows,by="accver",suffix=c("", ".window")) |>
  filter(end >= start.window, start <= end.window)
stopifnot(setequal(in_window$feature_id,fusion_detail_features$feature_id))
check_focus <- function(exported) {
  stopifnot(nrow(exported$leaders)==2L,all(exported$leaders$obstacle_collisions==0),
    all(exported$leaders$route_status=="direct_clear"),
    all(table(exported$leader_paths$boundary)==2L),
    all(exported$focus_regions$connector=="envelope"),
    all(exported$focus_regions$highlight_alpha==.35),
    all(exported$focus_regions$source_highlight_outline),
    all(exported$focus_regions$source_highlight_padding_mm==0),
    identical(as.data.frame(exported$focus_region_data[c("accver","start","end")]),
      as.data.frame(flagship_windows)))
  source<-exported$focus_fronts[exported$focus_fronts$panel_id=="center",]
  bands<-lapply(split(source,source$region_row),function(f) {
    # Include the visible outline stroke; no genomic padding is permitted.
    row<-exported$focus_regions[exported$focus_regions$region_row==f$region_row[1],]
    half<-row$highlight_width_mm/2
    polygon<-data.frame(x=c(f$x+half*f$normal_x,rev(f$x-half*f$normal_x)),
      y=c(f$y+half*f$normal_y,rev(f$y-half*f$normal_y)))
    attr(polygon,"stroke")<-row$leader_linewidth/2
    polygon
  })
  gate_hits<-vapply(bands,function(b) vapply(seq_len(nrow(exported$leaders)),function(i)
    {
      point<-c(exported$leaders$x[i],exported$leaders$y[i])
      a<-as.matrix(b);z<-a[c(2:nrow(a),1),];delta<-z-a
      along<-pmin(1,pmax(0,rowSums(sweep(a,2,point,"-")*(-delta))/pmax(rowSums(delta^2),.Machine$double.eps)))
      distance<-sqrt(rowSums(sweep(a+delta*along,2,point,"-")^2))
      ggchord_point_in_polygon(point[1],point[2],b) || min(distance)<=attr(b,"stroke")
    },logical(1)),logical(2))
  stopifnot(all(rowSums(gate_hits)==1L),all(colSums(gate_hits)==1L))
  d<-exported$routing_obstacles
  text<-d[d$panel_id=="zoom" & d$kind=="text",]
  heading<-d[d$panel_id=="zoom" & d$kind=="decoration",]
  contacts<-sum(vapply(seq_len(nrow(text)),function(i) sum(
    pmin(text$xmax[i],heading$xmax)>pmax(text$xmin[i],heading$xmin) &
    pmin(text$ymax[i],heading$ymax)>pmax(text$ymin[i],heading$ymin)),integer(1)))
  stopifnot(contacts==0L)
  frame_top<-max(exported$leaders$yend)
  frame_heading_gap<-min(heading$ymin)-frame_top
  detail_fronts<-exported$focus_fronts[exported$focus_fronts$panel_id=="zoom",]
  sequence_heading_gap<-min(heading$ymin)-max(detail_fronts$y)
  stopifnot(frame_heading_gap>0)
  headings<-d[d$panel_id %in% c("chord","radial") & d$kind=="decoration",]
  path<-exported$paths[exported$paths$layer_id %in% c("fusion-ribbon-1","fusion-ribbon-2"),]
  a_contacts<-sum(vapply(split(path,path$layer_id),function(p)
    sum(ggchord_route_polygon_hits(as.matrix(p[c("x","y")]),headings)),integer(1)))
  stopifnot(a_contacts==0L)
  list(leaders=nrow(exported$leaders),leader_bends=0L,leader_collisions=sum(exported$leaders$obstacle_collisions),
    heading_contacts=contacts,overview_ribbon_heading_contacts=a_contacts,
    detail_frame_excludes_heading=TRUE,frame_heading_gap_mm=frame_heading_gap,
    sequence_heading_gap_mm=sequence_heading_gap,
    source_gates_on_true_highlights=sum(gate_hits),source_interval_padding_mm=0,
    shared_guide_height_mm=exported$metadata$guide_bounds_mm$height_mm)
}
focus_check <- check_focus(fusion_export)
preview_export <- export_ggchord_fusion(home_fusion,width=10,height=8)
preview_focus_check <- check_focus(preview_export)
probe <- ggchord(flagship_sequences) + geom_seq(aes(seq_colour=accver)) +
  scale_seq_colour_manual(values=flagship_colours,name="Sequence") +
  theme_ggchord_publication() + theme(legend.position="bottom",legend.text=element_text(size=7)) +
  guides(seq_colour=guide_ggchord_legend(nrow=2))
shared <- ggchord_fusion(probe,guides="collect",center_box=c(0,.48,0,1)) + fusion_panel(probe,"other",box=c(.52,1,0,1))
shared_export <- export_ggchord_fusion(shared,width=12,height=8)
stopifnot(shared_export$metadata$guide_layout$collected$shared==1L)
changed <- flagship_colours;changed[1] <- "#000000"
conflict <- ggchord_fusion(probe,guides="collect",center_box=c(0,.48,0,1)) +
  fusion_panel(probe+scale_seq_colour_manual(values=changed,name="Sequence"),"other",box=c(.52,1,0,1))
stopifnot(export_ggchord_fusion(conflict,width=12,height=8)$metadata$guide_layout$collected$shared==2L)
kept <- export_ggchord_fusion(ggchord_fusion(probe,center_box=c(0,.48,0,1))+fusion_panel(probe,"other",box=c(.52,1,0,1)),width=12,height=8)
stopifnot(all(vapply(kept$diagnostics,function(d) length(d$guide_sizes_mm)>0L,logical(1))))
# Explicit mapped styles still win over automatic cyclic line types/widths.
manual <- ggchord(show_sequences, show_links) + geom_seq() +
  geom_link_ribbon(aes(ribbon_linetype=qaccver),adjacent_only=FALSE) +
  scale_ribbon_linetype_manual(values=setNames(rep("dotdash",4),show_sequences$accver)) +
  coord_radial_synteny()
stopifnot(all(ggplot2::ggplot_build(manual)$data[[2]]$ribbon_linetype=="dotdash"))
# Small physical-front regression: a grazing arrival must not acquire an S bend.
old_env <- new.env(parent=asNamespace('ggchord'))
previous <- parse(text=system2('git',c('show','7bcc1ad:R/fusion.R'),stdout=TRUE))
for(e in previous) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&
  identical(e[[2]],as.name('ggchord_fusion_path'))) eval(e,old_env)
previous_heading <- parse(text=system2('git',c('show','17db4c1:R/fusion-routing.R'),stdout=TRUE))
for(e in previous_heading) if(is.call(e)&&identical(e[[1]],as.name('<-'))&&
  identical(e[[2]],as.name('ggchord_fusion_guide_obstacles'))) eval(e,old_env)
export_with_old_headings <- function(plot,width,height) {
  current <- get('ggchord_fusion_guide_obstacles',asNamespace('ggchord'))
  on.exit(assignInNamespace('ggchord_fusion_guide_obstacles',current,ns='ggchord'))
  assignInNamespace('ggchord_fusion_guide_obstacles',old_env$ggchord_fusion_guide_obstacles,ns='ggchord')
  export_ggchord_fusion(plot,width=width,height=height)
}
q <- data.frame(x=seq(0,1,length.out=40),y=seq(0,-1,length.out=40))
z <- data.frame(x=seq(20,28,length.out=40),y=-.4)
changes <- function(p) {
 d<-p[-1,,drop=FALSE]-p[-nrow(p),,drop=FALSE]
 cross<-head(d[,1],-1)*tail(d[,2],-1)-head(d[,2],-1)*tail(d[,1],-1)
 sign<-sign(cross[abs(cross)>1e-8]);sum(head(sign,-1)!=tail(sign,-1))
}
old_path<-old_env$ggchord_fusion_path(q,z,.2);new_path<-ggchord_fusion_path(q,z,.2)
stopifnot(changes(old_path[41:100,])==1L,changes(new_path[41:100,])==0L,
  identical(old_path[1:40,],new_path[1:40,]),identical(old_path[101:140,],new_path[101:140,]))
# Check the existing composition illustrations at their declared dimensions.
# Only drawings with changed control normals are regenerated.
for(key in c('fusionpanels','fusionlinear','fusionview','fusionfour')) invisible(eval(parse(text=b[[key]])))
legacy <- list();course_equivalence <- list()
manifest <- read.csv('website/assets/plots.csv',stringsAsFactors=FALSE,check.names=FALSE)
fingerprints <- jsonlite::fromJSON(paste(system2('python3','tools/documentation-fingerprints.py',stdout=TRUE),collapse='\n'))
course_rows <- jsonlite::fromJSON('website/course-figures.json',simplifyVector=FALSE)
candidates <- c(lapply(which(manifest$id %in% c('fusionview','fusionfour')),function(i) as.list(manifest[i,])),
  Filter(function(r) 'fusiondata' %in% r$snippets,course_rows))
courses <- jsonlite::fromJSON('tools/site-audit/course-rendering.json',simplifyVector=FALSE)
prior <- if(file.exists('tools/site-audit/radial-update.json'))
  jsonlite::fromJSON('tools/site-audit/radial-update.json',simplifyVector=FALSE) else NULL
legacy_inputs <- c(names(snapshot$engine),'tools/documentation-examples.R',
  'tools/site-course-examples.R','website/course-figures.json')
reuse <- !is.null(prior$legacy_geometry) && all(vapply(legacy_inputs,function(name)
  identical(prior$source_sha256[[name]],snapshot$inputs[[name]]),logical(1)))
if(reuse) {
  legacy <- prior$figures[setdiff(names(prior$figures),c('flagradial','homeradial','flagfusion','homefusion'))]
  reuse <- all(vapply(legacy,function(item) identical(item$sha256,
    digest::digest(file=item$file,algo='sha256')),logical(1)))
}
if(reuse) {
  course_equivalence <- prior$legacy_geometry
  for(name in names(legacy)) legacy[[name]]$retained_generation <- TRUE
} else {
for(row in candidates) {
 if(!is.null(row$snippets)) for(key in row$snippets) invisible(eval(parse(text=b[[key]])))
 plot <- get(row$object)
 if(!inherits(plot,'ggchord_fusion')) next
 routes<-vapply(plot$ribbons,function(r) r$route %||% "direct",character(1))
 exported <- export_ggchord_fusion(plot,width=row$width,height=row$height)
 changed_routes <- 0L
 if(any(routes!="direct") || (!is.null(plot$focus) && !identical(plot$focus$route %||% "direct","direct"))) {
   before <- export_with_old_headings(plot,row$width,row$height)
   changed_routes <- as.integer(!identical(exported$paths,before$paths) ||
     !identical(exported$leader_paths,before$leader_paths))
 }
 changed <- 0L
 for(i in seq_len(nrow(exported$connections))) {
  connection<-exported$connections[i,];if(!isTRUE(connection$retained)) next
  f<-exported$fronts[exported$fronts$layer_id==connection$layer_id & exported$fronts$source_row==connection$source_row,]
  q<-f[f$endpoint=='q' & f$panel_id==connection$qpanel,];z<-f[f$endpoint=='s' & f$panel_id==connection$spanel,]
  old<-old_env$ggchord_fusion_path(q,z,plot$ribbons[[match(connection$layer_id,paste0('fusion-ribbon-',seq_along(plot$ribbons)))]]$curvature %||% .3)
  now<-ggchord_fusion_path(q,z,plot$ribbons[[match(connection$layer_id,paste0('fusion-ribbon-',seq_along(plot$ribbons)))]]$curvature %||% .3)
  if(!identical(attr(old,'endpoint_normals'),attr(now,'endpoint_normals'))) changed<-changed+1L
 }
 course_equivalence[[row$id]] <- list(connections=nrow(exported$connections),changed_controls=changed,
   changed_obstacle_routes=changed_routes,
   retained_png_sha256=if(changed==0L && changed_routes==0L) digest::digest(file=row$file,algo='sha256') else NULL)
 if(changed==0L && changed_routes==0L) next
 ggsave(row$file,plot,width=row$width,height=row$height,dpi=as.numeric(row$dpi),bg='white')
 item<-c(row,list(sha256=digest::digest(file=row$file,algo='sha256'),status='actual PNG rendered'))
 if(row$id %in% manifest$id) {
  i<-match(row$id,manifest$id);manifest$md5[i]<-digest::digest(file=row$file,algo='md5')
  manifest$dependency_md5[i]<-fingerprints$assets[[row$id]]
  manifest$generator[i]<-'tools/site-audit/render-radial-update.R'
  item$md5<-manifest$md5[i];item$dependency_md5<-manifest$dependency_md5[i]
 } else {
  if(is.null(courses$figures[[row$id]]$previous_render))
    courses$figures[[row$id]]$previous_render<-courses$figures[[row$id]][c('sha256','rendered_source_sha256')]
  courses$figures[[row$id]]$sha256<-item$sha256
  courses$figures[[row$id]]$rendered_source_sha256<-as.list(snapshot$inputs)
 }
 legacy[[row$id]]<-item
}
write.csv(manifest,'website/assets/plots.csv',row.names=FALSE,na='')
jsonlite::write_json(courses,'tools/site-audit/course-rendering.json',auto_unbox=TRUE,pretty=TRUE)
}
rows <- c(jsonlite::fromJSON('website/manual-figures.json',simplifyVector=FALSE),
  jsonlite::fromJSON('website/home-previews.json',simplifyVector=FALSE))
figures <- legacy
for(row in rows) if(row$id %in% c('flagradial','homeradial','flagfusion','homefusion')) {
  ggsave(row$file,get(row$object),width=row$width,height=row$height,dpi=row$dpi,bg='#FFFCF5')
  figures[[row$id]] <- c(row,list(sha256=digest::digest(file=row$file,algo='sha256'),status='actual PNG rendered'))
}
audit_assert_source_snapshot(snapshot)
record <- list(status='complete',baseline='4f68edf',source_sha256=as.list(snapshot$inputs),
  engine_sha256=as.list(snapshot$engine),figures=figures,
  checks=list(HSP=12,CDS=20,nonadjacent=4,backbone_clearance_min=min(clearances),
    source_fields_identical=TRUE,fixed_linetype_priority=TRUE,nonadjacent_linewidth=.4,fusion_connections=3,fusion_HSP=32,fusion_CDS=266,
    detail_HSP_length=flagship_selected$length,mapped_linetype_priority=TRUE,
    fusion_display_HSP=6,fusion_display_CDS=12,shared_child_guides=3,
    child_guides_removed=TRUE,distinct_guide_contents_retained=TRUE,keep_default_retained=TRUE,detail_forward=TRUE,zoom_panels=1,focus_windows=2,
    magnification_ratio_min=min(fusion_export$focus_regions$length_ratio),
    all_five_panels_converged=TRUE,grazing_double_bend_fixed=TRUE,
    context_detail_feature_ids_identical=TRUE,detail_source_features=4),
  focus_layout=focus_check,preview_focus_layout=preview_focus_check,
  legacy_geometry_reused=reuse,legacy_geometry=course_equivalence,
  route_bulge=ends[,c('qaccver','saccver','source_row','route_bulge')],
  scope='Radial normal-size PNG update; no browser or package-wide release checks')
jsonlite::write_json(record,'tools/site-audit/radial-update.json',auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('PASS: source records retained; backbone clearance and fixed style priority; four current PNGs rendered\n')
