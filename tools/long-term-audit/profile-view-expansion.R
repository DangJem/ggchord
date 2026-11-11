# Bounded synthetic replication benchmark; not additional biological evidence.
# Run explicitly; no time threshold, no release output/device matrix.
pkgload::load_all(quiet = TRUE)
library(ggplot2)
h <- read_paf("examples/format-inputs/emales-pair.paf")
s <- unique(rbind(data.frame(accver=h$qaccver,length=h$qlen),
  data.frame(accver=h$saccver,length=h$slen)))
rows <- list()
for (mode in c("dotplot", "radial")) {
  n <- if (mode == "dotplot") 10000L else 1000L
  input <- h[rep(seq_len(nrow(h)), length.out=n), ]
  p <- ggchord(s,input,validate="none") + if (mode == "dotplot")
    list(geom_alignment_segment(),coord_dotplot()) else
    list(geom_seq(),geom_link_ribbon(),coord_radial_synteny())
  elapsed <- system.time({b <- ggplot_build(p); g <- ggplot_gtable(b)})[["elapsed"]]
  layout <- attr(g,"ggchord_render_layout")
  shown <- if (mode == "dotplot") nrow(layout$alignment[[1]]) else nrow(layout$radial_endpoints)
  stopifnot(shown == n)
  rows[[mode]] <- list(mode=mode,input_rows=n,displayed_rows=shown,elapsed_seconds=elapsed,
    built_object_MiB=as.numeric(object.size(b))/1024^2,
    note="Explicit synthetic repetition of real PAF rows; object size is not RSS; no fixed time pass criterion")
}
jsonlite::write_json(rows,"tools/long-term-audit/view-expansion-profile.json",auto_unbox=TRUE,pretty=TRUE)
