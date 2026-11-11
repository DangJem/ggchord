# Direct Composition checks, not package/release acceptance.
pkgload::load_all(quiet = TRUE)
hits <- read_paf("examples/format-inputs/emales-pair.paf")
seqs <- unique(rbind(data.frame(accver = hits$qaccver, length = hits$qlen),
  data.frame(accver = hits$saccver, length = hits$slen)))
center <- ggchord(seqs[1, , drop = FALSE], validate = "none") +
  geom_seq(show.legend = FALSE) + coord_circular()
linear <- ggchord(seqs, hits, validate = "none") +
  geom_seq(show.legend = FALSE) + geom_link_ribbon(show.legend = FALSE) + coord_collinear()
vertical <- linear + coord_collinear(seq_axis = "vertical")
f0 <- ggchord_fusion(center) + fusion_panel(vertical, "detail") +
  fusion_ribbon(hits, qpanel = "center", spanel = "detail")
cache_before <- serialize(center$ggchord$ref, NULL)
e <- export_ggchord_fusion(f0, 10, 8)
stopifnot(identical(cache_before, serialize(center$ggchord$ref, NULL)),
  nrow(e$connections) == nrow(hits), all(e$connections$retained),
  identical(e$relationship_data[[1]], hits))
cols <- c("qstart", "qend", "sstart", "send")
stopifnot(isTRUE(all.equal(e$connections[cols], hits[cols], check.attributes = FALSE)))
check_fronts <- function(e) {
  for (rows in split(seq_len(nrow(e$fronts)), interaction(e$fronts$layer_id,
      e$fronts$source_row, e$fronts$endpoint, drop = TRUE))) {
    f <- e$fronts[rows, ]
    c <- e$connections[e$connections$layer_id == f$layer_id[1] & e$connections$source_row == f$source_row[1], ]
    prefix <- f$endpoint[1]; a <- c[[paste0(prefix, "start")]]; b <- c[[paste0(prefix, "end")]]
    expected <- c(min(a,b)-1, max(a,b)); if(b<a) expected <- rev(expected)
    stopifnot(isTRUE(all.equal(f$position_boundary[c(1,nrow(f))],expected,check.attributes=FALSE)))
    # Check final-coordinate projection against actual physical panel ranges.
    panel <- e$panels[e$panels$panel_id == f$panel_id[1], ]
    context <- e$diagnostics[[f$panel_id[1]]]$render_context
    projected_x <- panel$panel_xmin + (f$local_x-context$x_range[1])/diff(context$x_range) *
      (panel$panel_xmax-panel$panel_xmin)
    projected_y <- panel$panel_ymin + (f$local_y-context$y_range[1])/diff(context$y_range) *
      (panel$panel_ymax-panel$panel_ymin)
    stopifnot(isTRUE(all.equal(f$x,unname(projected_x),check.attributes=FALSE)),
      isTRUE(all.equal(f$y,unname(projected_y),check.attributes=FALSE)))
  }
}
check_fronts(e)
resized <- export_ggchord_fusion(f0, 12, 6)
check_fronts(resized)
stopifnot(identical(resized$connections[cols],e$connections[cols]),
  !isTRUE(all.equal(resized$fronts$x,e$fronts$x)))
# Four physical directions and independent data layers. Repeated source IDs
# are opt-in; row-local identity does not guess global identity by value.
four <- ggchord_fusion(center, layout="adaptive") +
  fusion_panel(vertical,"east","right") + fusion_panel(vertical,"west","left") +
  fusion_panel(linear,"north","top") + fusion_panel(linear,"south","bottom")
for (id in c("east","west","north","south")) four <- four +
  fusion_ribbon(hits[2,,drop=FALSE],qpanel="center",spanel=id,fill="#26738c",show.legend=FALSE)
ef <- export_ggchord_fusion(four,14,11)
stopifnot(nrow(ef$panels)==5L,nrow(ef$connections)==4L)
check_fronts(ef)
# Chord central frames, single positions, backwards orientation and explicit
# mappings/scales. Synthetic intervals are labelled as such, not new evidence.
chord <- ggchord(seqs,validate="none") + geom_seq(show.legend=FALSE) + coord_chord(rotation=20)
single <- hits[2,,drop=FALSE]; single$qend <- single$qstart; single$send <- single$sstart
single$direction <- NULL; single$relative_strand <- NULL
f1 <- ggchord_fusion(chord) + fusion_panel(linear,"detail","bottom") +
  fusion_ribbon(single,ggplot2::aes(fill=pident),qpanel="center",spanel="detail") +
  ggplot2::scale_fill_gradient(low="black",high="white",limits=c(0,100))
es <- export_ggchord_fusion(f1,10,8)
check_fronts(es)
stopifnot(all(vapply(split(es$fronts,es$fronts$endpoint),function(f)
  diff(range(f$x))+diff(range(f$y))>0,logical(1))))
# Explicit instance namespace within a repeated-source peripheral panel.
repeat_plot <- ggchord(seqs,validate="none") +
  geom_seq(seq_order=stats::setNames(c(seqs$accver[1],seqs$accver[2],seqs$accver[1]),c("A1","B","A2")),show.legend=FALSE) + coord_collinear()
route <- single; route$saccver <- seqs$accver[1]; route$sstart <- 100; route$send <- 200
route$s_instance <- "A2"
fr <- ggchord_fusion(center) + fusion_panel(repeat_plot,"repeat") +
  fusion_ribbon(route,qpanel="center",spanel="repeat",fill="grey70")
er <- export_ggchord_fusion(fr)
stopifnot(er$connections$s_instance=="A2")
# Empty relationships stay empty and carry no phantom endpoints.
empty <- ggchord_fusion(center) + fusion_panel(linear,"detail") +
  fusion_ribbon(hits[FALSE,],qpanel="center",spanel="detail")
stopifnot(nrow(export_ggchord_fusion(empty)$connections)==0L)
expect_error <- function(expr, pattern) {
  err <- tryCatch({force(expr); NULL},error=identity)
  stopifnot(inherits(err,"error"),grepl(pattern,conditionMessage(err)))
}
expect_error(export_ggchord_fusion(ggchord_fusion(chord)+fusion_panel(linear,"detail")+
  fusion_ribbon(hits[1,,drop=FALSE])),"ambiguous")
expect_error(ggchord_fusion(center)+fusion_panel(linear,"a")+fusion_panel(linear,"b"),"one panel per side")
expect_error(export_ggchord_fusion(ggchord_fusion(center)+fusion_panel(linear,"detail")+
  fusion_ribbon(hits[1,,drop=FALSE],qpanel="bad",spanel="detail")),"unknown panel")
expect_error(fusion_panel(linear+coord_sankey(),"no"),"Collinear")
expect_error(fusion_panel(linear+coord_radial_synteny(),"no"),"Collinear")
expect_error(export_ggchord_fusion(f0,2,2),"canvas|guides")
# Optional PNG output is local visual evidence, not a release device matrix.
output <- commandArgs(trailingOnly=TRUE)
if(length(output)) {
  ggplot2::ggsave(file.path(output[1],"fusion-f0.png"),f0,width=11,height=8,dpi=120,bg="white")
  ggplot2::ggsave(file.path(output[1],"fusion-f1.png"),four,width=14,height=11,dpi=120,bg="white")
  gradient <- ggchord_fusion(center)+fusion_panel(vertical,"detail")+
    fusion_ribbon(hits[2,,drop=FALSE],qpanel="center",spanel="detail",gradient=TRUE,tip="subject")
  ggplot2::ggsave(file.path(output[1],"fusion-gradient.png"),gradient,width=10,height=7,dpi=110,bg="white")
}
message("Fusion: real HSP fronts, physical resize, four directions, Chord/Circular, instances, styles and failures passed")
# Panel-local source statistics, labels and explicit theme/scale stay local.
annotations <- data.frame(accver=seqs$accver[1],start=1000,end=1500,strand="+",type="CDS",anno="source-backed interval demo")
# Above annotation is synthetic and is used only to check the layer contract.
events <- data.frame(accver=hits$qaccver,pos=hits$qstart)
annotated <- ggchord(seqs[1,,drop=FALSE],feature_data=annotations,validate="none")+
  geom_seq(show.legend=FALSE)+geom_feature(shape="arrow")+geom_feature_label()+
  geom_track(data=events,geom="bar",track=position_track("events",width=grid::unit(2,"mm")))+coord_circular()
fa <- ggchord_fusion(annotated)+fusion_panel(vertical,"detail")+
  fusion_ribbon(hits[3,,drop=FALSE],qpanel="center",spanel="detail",fill="#0072B2")+
  ggplot2::guides(fill="none")+ggplot2::labs(fill="Cross-panel identity")+ggplot2::theme(legend.position="bottom")
ea <- export_ggchord_fusion(fa)
stopifnot(nrow(ea$source_layouts$center$track_registry)>0,
  nrow(ea$diagnostics$center$label_visibility)>0)
# An unknown orientation is explicitly omitted and retained in the report.
u <- hits[1,,drop=FALSE];u$direction <- "unknown"
eu <- suppressWarnings(export_ggchord_fusion(ggchord_fusion(center)+fusion_panel(linear,"detail")+
  fusion_ribbon(u,qpanel="center",spanel="detail")))
stopifnot(!eu$connections$retained,eu$connections$reason=="unknown_direction",nrow(eu$fronts)==0)
# A declared HSP cannot be duplicated across independent layers by accident.
id <- hits[3,,drop=FALSE];id$source_hsp_id <- "real-paf:3"
duplicate <- ggchord_fusion(center)+fusion_panel(linear,"detail")+
  fusion_ribbon(id,qpanel="center",spanel="detail",fill="grey")+
  fusion_ribbon(id,qpanel="center",spanel="detail",fill="grey")
expect_error(export_ggchord_fusion(duplicate),"allow_repeated")
# Shared preview is the same grid.draw path, not a different composition engine.
preview <- view_ggchord(fa,width=10,height=NULL,viewer="none",dpi=72)
stopifnot(file.exists(preview))
# Existing single-plot build and label export remain available.
old <- export_ggchord_layout(annotated,original_data=TRUE)
stopifnot(nrow(old$feature)>0,nrow(old$label_visibility)>0)
message("Fusion: source Track/labels, unknown reports, duplicate IDs, native guides/themes and preview passed")
# Source accession and source coordinates survive explicitly renamed windows.
loci <- data.frame(accver=seqs$accver[2],start=3500,end=5000,locus_id="subject-window")
focused <- focus_ggchord_data(seqs,loci=loci)
window_plot <- ggchord(focused$seq_data,validate="none")+geom_seq(show.legend=FALSE)+coord_collinear()
fw <- ggchord_fusion(center)+fusion_panel(window_plot,"window")+
  fusion_ribbon(hits[9,,drop=FALSE],qpanel="center",spanel="window",fill="grey70")
ew <- export_ggchord_fusion(fw)
stopifnot(ew$connections$s_instance=="subject-window",
  min(ew$fronts$position_boundary[ew$fronts$endpoint=="s"])==min(hits$sstart[9],hits$send[9])-1)
check_fronts(ew)
message("Fusion: explicitly renamed source windows passed")
