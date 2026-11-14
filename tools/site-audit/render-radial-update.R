source('tools/audit-source-snapshot.R')
snapshot <- audit_source_snapshot(additional=c('tools/site-flagships.R','tools/home-previews.R','tools/documentation-examples.R','tools/site-audit/render-radial-update.R'))
pkgload::load_all(quiet=TRUE)
library(ggplot2)
library(dplyr)
library(grid)
blocks <- function(file) {
 x <- readLines(file,warn=FALSE); m <- grep('^# ---- [a-z]+ ----$',x)
 setNames(lapply(seq_len(length(m)-1),function(i) paste(x[seq.int(m[i]+1,m[i+1]-1)],collapse='\n')),sub('^# ---- ([a-z]+) ----$','\\1',x[m[-length(m)]]))
}
b <- c(blocks('tools/documentation-examples.R'),blocks('tools/site-flagships.R'))
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
  nrow(flagship_selected)==1L,flagship_selected$length<3000)
# Explicit mapped styles still win over automatic cyclic line types/widths.
manual <- ggchord(show_sequences, show_links) + geom_seq() +
  geom_link_ribbon(aes(ribbon_linetype=qaccver),adjacent_only=FALSE) +
  scale_ribbon_linetype_manual(values=setNames(rep("dotdash",4),show_sequences$accver)) +
  coord_radial_synteny()
stopifnot(all(ggplot2::ggplot_build(manual)$data[[2]]$ribbon_linetype=="dotdash"))
rows <- c(jsonlite::fromJSON('website/manual-figures.json',simplifyVector=FALSE),
  jsonlite::fromJSON('website/home-previews.json',simplifyVector=FALSE))
figures <- list()
for(row in rows) if(row$id %in% c('flagradial','homeradial','flagfusion','homefusion')) {
  ggsave(row$file,get(row$object),width=row$width,height=row$height,dpi=row$dpi,bg='#FFFCF5')
  figures[[row$id]] <- c(row,list(sha256=digest::digest(file=row$file,algo='sha256'),status='actual PNG rendered'))
}
audit_assert_source_snapshot(snapshot)
record <- list(status='complete',baseline='4f68edf',source_sha256=as.list(snapshot$inputs),
  engine_sha256=as.list(snapshot$engine),figures=figures,
  checks=list(HSP=12,CDS=20,nonadjacent=4,backbone_clearance_min=min(clearances),
    source_fields_identical=TRUE,fixed_linetype_priority=TRUE,nonadjacent_linewidth=.4,fusion_connections=3,fusion_HSP=32,fusion_CDS=266,
    detail_HSP_length=flagship_selected$length,mapped_linetype_priority=TRUE),
  route_bulge=ends[,c('qaccver','saccver','source_row','route_bulge')],
  scope='Radial normal-size PNG update; no browser or package-wide release checks')
jsonlite::write_json(record,'tools/site-audit/radial-update.json',auto_unbox=TRUE,pretty=TRUE,digits=16)
cat('PASS: source records retained; backbone clearance and fixed style priority; four current PNGs rendered\n')
