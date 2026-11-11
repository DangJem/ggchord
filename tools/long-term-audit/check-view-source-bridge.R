# One real normalized PAF table across all supported relationship views.
# Compare biological source rows/endpoints, never screen-coordinate equality.
pkgload::load_all(quiet = TRUE)
hits <- read_paf("examples/format-inputs/emales-pair.paf")
seqs <- unique(rbind(data.frame(accver=hits$qaccver,length=hits$qlen),
  data.frame(accver=hits$saccver,length=hits$slen)))
coords <- list(chord=coord_chord(),collinear=coord_collinear(),sankey=coord_sankey(),
  dotplot=coord_dotplot(),radial=coord_radial_synteny())
for (mode in names(coords)) {
  p <- ggchord(seqs,hits,validate="none")
  p <- if (mode == "dotplot") p + geom_alignment_segment() else
    p + geom_seq() + geom_link_ribbon()
  p <- p + coords[[mode]]
  e <- export_ggchord_layout(p,original_data=TRUE)
  shown <- if (mode=="dotplot") e$alignment else e$ribbon
  stopifnot(setequal(unique(shown$source_row),seq_len(nrow(hits))),
    identical(p$ggchord$sources$relationship,hits))
  input <- Filter(function(x)!is.null(x$alignment %||% x$ribbon),e$original_data)[[1]]
  input <- input$alignment %||% input$ribbon
  cols <- c("qstart","qend","sstart","send")
  stopifnot(isTRUE(all.equal(input[cols],hits[cols],check.attributes=FALSE)))
  invisible(ggplot2::ggplotGrob(p))
  message(mode, ": all ",nrow(hits)," real source HSPs and four endpoints retained")
}
