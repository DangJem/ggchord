# Direct coordinate checks for focused composite annotations.
# Run from the repository root: Rscript tools/validate-focus-coordinates.R
pkgload::load_all(quiet = TRUE)

sequences <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
features <- data.frame(accver = "A", start = 400, end = 700,
  strand = "+", type = "CDS", label = "segmented")
features$segments <- I(list(data.frame(segment_index = 1:2,
  segment_type = "standard", start = c(400, 601), end = c(500, 700),
  length_bp = c(101, 100), color = "#4488AA")))
features$cleavage_arrows <- I(list(c(425, 475, 625, 675)))
# Cached render fields are deliberately stale; the source segments govern.
features$.feature_runs <- I(list(data.frame(start = 400, end = 700,
  segment_index = 1L, color = "#4488AA", draw_head = TRUE,
  draw_start_head = FALSE, boundaries = I(list(475)),
  boundary_styles = I(list("dotted")))))
loci <- data.frame(accver = c("A", "B"), start = c(450, 200), end = c(650, 800))
links <- data.frame(qaccver = "A", saccver = "B", qstart = 400, qend = 700,
  sstart = 900, send = 100, length = 801, pident = 95)

focused <- focus_ggchord_data(sequences, links, features, loci = loci)
segment <- focused$gene_data$segments[[1L]]
stopifnot(all(segment$start >= 1), all(segment$end <= focused$seq_data$length[1L]),
  identical(as.numeric(segment$start), c(1, 152)),
  identical(as.numeric(segment$end), c(51, 201)),
  identical(as.numeric(focused$gene_data$cleavage_arrows[[1L]]), c(26, 176)),
  !".feature_runs" %in% names(focused$gene_data),
  focused$gene_data$.truncated_left, focused$gene_data$.truncated_right,
  focused$gene_data$.source_start == features$start,
  focused$gene_data$.source_end == features$end,
  focused$seq_data$.source_origin[1L] == loci$start[1L],
  all(focused$ribbon_data$.alignment_projection == "linear"),
  all(focused$ribbon_data$.length_estimated),
  identical(focused$ribbon_data$pident, links$pident),
  all(focused$ribbon_data$.source_pident == links$pident))

# A second focus composes source coordinates rather than overwriting origin.
again <- focus_ggchord_data(focused$seq_data, focused$ribbon_data,
  focused$gene_data, loci = data.frame(accver = c("A", "B"),
    start = c(11, 11), end = c(191, 591)))
stopifnot(again$seq_data$.source_origin[1L] == loci$start[1L] + 10,
  again$gene_data$.source_start == features$start,
  again$gene_data$.source_end == features$end,
  all(again$gene_data$segments[[1L]]$end <= again$seq_data$length[1L]),
  all(again$ribbon_data$.source_qstart == links$qstart),
  all(again$ribbon_data$.source_sstart == links$sstart))

# The focused public rendering stays within its biological sequence interval.
plot <- ggchord(focused$seq_data, gene_data = focused$gene_data, validate = "none") +
  geom_seq() + geom_feature() + coord_collinear()
geometry <- export_ggchord_layout(plot, include = "feature")
sequence_width <- geometry$sequence_layout$xmax[
  geometry$sequence_layout$accver == focused$gene_data$accver[1L]]
stopifnot(all(is.finite(geometry$feature$x)),
  all(geometry$feature$x >= 0 & geometry$feature$x <= sequence_width))

# Already-derived interval runs are handled even without raw segments.
run_features <- features
run_features$segments <- NULL
run_focused <- focus_ggchord_data(sequences, gene_data = run_features, loci = loci)
run <- run_focused$gene_data$.feature_runs[[1L]]
stopifnot(run$start == 1, run$end == 201,
  identical(as.numeric(run$boundaries[[1L]]), 26))

# Cropping a gap-only part does not turn it into a solid false interval.
gap_focus <- focus_ggchord_data(sequences, gene_data = features,
  loci = data.frame(accver = "A", start = 550, end = 580))
stopifnot(nrow(gap_focus$gene_data) == 0L,
  gap_focus$report$dropped_input_rows[gap_focus$report$component == "gene"] == 1L)

dropped <- focus_ggchord_data(sequences, links, features, loci = loci, boundary = "drop")
stopifnot(nrow(dropped$gene_data) == 0L, nrow(dropped$ribbon_data) == 0L)
cat("Focused composite coordinates and provenance passed.\n")
