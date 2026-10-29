# Explicit local development validation, not a testthat or release run.
pkgload::load_all(quiet = TRUE)
seqs <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
windows <- data.frame(accver = c("A", "A", "B", "B"),
  start = c(100, 700, 100, 700), end = c(300, 900, 300, 900),
  locus_id = c("A1", "A2", "B1", "B2"), row_id = c("A", "A", "B", "B"),
  orientation = c(1, -1, 1, 1))
links <- data.frame(qaccver = "A", saccver = "B", qstart = 100,
  qend = 900, sstart = 100, send = 900, pident = 95, length = 801)
features <- data.frame(accver = "A", start = 150, end = 850,
  strand = "+", feature_id = "spanning", anno = "spanning feature")
points <- data.frame(qaccver = "A", saccver = "B", qpos = c(200, 500), spos = c(200, 500))
view <- segment_ggchord_data(seqs, links, features, windows, point_data = points)
stopifnot(nrow(view$ribbon_data) == 2L, nrow(view$point_data) == 1L,
  all(view$ribbon_data$.alignment_clipped),
  identical(view$report$omitted_point_rows, 2L),
  all(view$gene_data$feature_id == "spanning"))
p <- ggchord(view$seq_data, view$ribbon_data, validate = "none") +
  geom_seq() + geom_link_ribbon() + geom_feature(data = view$gene_data) +
  geom_feature_label_repel(data = view$gene_data) +
  geom_link_line(data = view$point_data) +
  coord_collinear(segments = view$segments, segment_gap = 50)
out <- export_ggchord_layout(p, include = c("seq", "ribbon", "feature", "link"))
layout <- out$sequence_layout
stopifnot(layout$y[1] == layout$y[2], layout$y[3] == layout$y[4], layout$y[1] != layout$y[3],
  layout$xmax[1] < layout$xmin[2], layout$orientation[2] == -1,
  any(out$seq$.decoration %in% "segment_break"),
  nrow(out$ribbon) > 0L, nrow(out$link) > 0L,
  all(out$connection_report$retained))
mapped <- map_ggchord_position(p, "A2", 800, source = TRUE)
inverse <- map_ggchord_position(p, "A2", mapped$x, source = TRUE, inverse = TRUE)
stopifnot(abs(inverse$position - 800) < 1e-7)
gap <- tryCatch(map_ggchord_position(p, "A1", mean(c(layout$xmax[1], layout$xmin[2])),
  inverse = TRUE), error = identity)
stopifnot(inherits(gap, "error"))
dropped <- segment_ggchord_data(seqs, links, windows = windows, cross_gap = "drop")
stopifnot(nrow(dropped$ribbon_data) == 0L, identical(dropped$report$omitted_ribbon_rows, 1L))
# Count all pairs that focus will allocate even when a whole-link policy
# subsequently filters partial windows.
overlapping <- windows[rep(1L, 3L), ]
overlapping$locus_id <- paste0("overlap", 1:3)
overlapping$start <- c(100, 110, 120)
guard <- tryCatch(segment_ggchord_data(seqs,
  data.frame(qaccver = "A", saccver = "A", qstart = 100, qend = 200, sstart = 100, send = 200),
  windows = overlapping, cross_gap = "drop", max_links = 1), error = identity)
stopifnot(inherits(guard, "error"), grepl("max_links", conditionMessage(guard)))
# Repeated focusing must preserve the earliest point provenance.
again <- segment_ggchord_data(view$seq_data, windows = data.frame(
  accver = c("A1", "B1"), start = 20, end = 150,
  locus_id = c("againA", "againB")), point_data = view$point_data)
stopifnot(again$point_data$.source_qaccver == "A", again$point_data$.source_qpos == 200)
stopifnot(inherits(ggplot2::ggplot_build(p), "ggplot_built"))
cat("Segmented windows, local reversals, gap routing and inverse mapping passed.\n")
