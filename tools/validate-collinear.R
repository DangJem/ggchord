# Focused development checks for linear sequence layouts.
# Run from the repository root:
#   Rscript tools/validate-collinear.R [optional-preview.png]
# This script does not run package tests, package checks or a device matrix.
pkgload::load_all(quiet = TRUE)
source(file.path("examples", "collinear.R"), local = TRUE)

assert_finite_geometry <- function(exported) {
  geometry <- Filter(function(x) is.data.frame(x) &&
    all(c("x", "y") %in% names(x)) && nrow(x) > 0L, exported)
  stopifnot(length(geometry) > 0L)
  for (part in geometry) {
    # Mixed text/leader exports retain their component-specific anchors.
    valid <- is.finite(part$x) & is.finite(part$y)
    if (all(c("text_x", "text_y") %in% names(part))) {
      valid <- valid | (is.finite(part$text_x) & is.finite(part$text_y))
    }
    stopifnot(all(valid))
  }
  invisible(exported)
}
expect_conflict <- function(expr, parameter) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"),
    grepl(parameter, conditionMessage(error), fixed = TRUE))
}

# One representative composed build exercises the public rendering path.
built <- ggplot2::ggplot_build(collinear_example)
stopifnot(inherits(built, "ggplot_built"))
out <- export_ggchord_layout(collinear_example,
  include = c("seq", "ribbon", "gene", "feature", "link", "labels", "axis"),
  original_data = TRUE)
assert_finite_geometry(out)
stopifnot(identical(out$metadata$coordinate, "collinear"),
  identical(as.character(out$sequence_layout$accver),
    c("cluster_C", "cluster_A", "cluster_B")),
  all(diff(out$sequence_layout$y) < 0),
  !"annotation_registry" %in% names(out),
  !any(c("radius", "track") %in% names(out$sequence_layout)))

# Source rows and inversion survive adjacent and cross-sequence alignments.
stopifnot(setequal(unique(out$ribbon$source_row), seq_len(nrow(collinear_alignments))),
  setequal(unique(out$gene$source_row), seq_len(nrow(collinear_genes))))
gene <- unique(out$gene[c("source_row", "strand", "normal_offset")])
stopifnot(all(gene$normal_offset[gene$strand == "+"] > 0),
  all(gene$normal_offset[gene$strand == "-"] < 0))
point_paths <- split(out$link, out$link$source_row)
stopifnot(length(point_paths) == nrow(collinear_points),
  point_paths[[1L]]$x[1L] < point_paths[[2L]]$x[1L],
  tail(point_paths[[1L]]$x, 1L) > tail(point_paths[[2L]]$x, 1L))
stacked <- out$feature[out$feature$component == "gene_poly", , drop = FALSE]
stopifnot(length(unique(stacked$lane)) > 1L,
  all(stacked$normal_offset > 0))
axis_labels <- out$axis$label[!is.na(out$axis$label)]
stopifnot(any(grepl("bp", axis_labels, fixed = TRUE)))

# Shared collision handling retains every label in this small representative map.
# Layer/source-row identities distinguish the three repeated "alpha" names.
fixed_text <- out$labels[out$labels$component == "gene_text" &
  !is.na(out$labels$text) & nzchar(out$labels$text),
  c("layer_id", "source_row", "text"), drop = FALSE]
feature_text <- out$labels[out$labels$component == "gene_label_repel" &
  out$labels$.component %in% "text" &
  !is.na(out$labels$label) & nzchar(out$labels$label),
  c("layer_id", "source_row", "label"), drop = FALSE]
fixed_text <- unique(fixed_text)
feature_text <- unique(feature_text)
stopifnot(length(unique(fixed_text$layer_id)) == 1L,
  nrow(fixed_text) == nrow(collinear_genes),
  setequal(fixed_text$source_row, seq_len(nrow(collinear_genes))),
  identical(as.character(fixed_text$text[match(seq_len(nrow(collinear_genes)),
    fixed_text$source_row)]), as.character(collinear_genes$anno)),
  length(unique(feature_text$layer_id)) == 1L,
  nrow(feature_text) == nrow(collinear_features),
  setequal(feature_text$source_row, seq_len(nrow(collinear_features))),
  identical(as.character(feature_text$label[match(seq_len(nrow(collinear_features)),
    feature_text$source_row)]), as.character(collinear_features$label)))

# Additional layers retain their own data functions, mapping and source rows.
separate <- ggchord(collinear_sequences, validate = "none") +
  geom_seq() +
  geom_feature(data = collinear_features[1L, ], feature_shape = "block") +
  geom_feature(data = function(x) collinear_features[2L, ],
    feature_shape = "chevron") + coord_collinear()
independent <- export_ggchord_layout(separate, include = "feature",
  original_data = TRUE)
stopifnot(length(unique(independent$feature$layer_id)) == 2L,
  setequal(independent$feature$feature_shape, c("block", "chevron")))

# The same source tables also remain valid in the default chord coordinate.
chord <- ggchord(collinear_sequences, collinear_alignments, collinear_genes,
  validate = "none") + geom_seq() + geom_link_ribbon() + geom_gene()
assert_finite_geometry(export_ggchord_layout(chord,
  include = c("seq", "ribbon", "gene")))

# Focus returns ordinary relocated tables; linear rendering needs no new input API.
focused <- focus_ggchord_data(collinear_sequences, collinear_alignments,
  collinear_genes, loci = data.frame(accver = collinear_sequences$accver,
    start = c(80, 90, 100), end = c(800, 800, 900)))
focused_plot <- ggchord(focused$seq_data, focused$ribbon_data,
  focused$gene_data, validate = "none") + geom_seq() +
  geom_link_ribbon() + geom_gene() + coord_collinear()
focused_out <- export_ggchord_layout(focused_plot,
  include = c("seq", "ribbon", "gene"), original_data = TRUE)
assert_finite_geometry(focused_out)
stopifnot(setequal(focused_out$sequence_layout$accver, focused$seq_data$accver),
  nrow(focused_out$gene) > 0L)

# Both composition orders must reject chord-specific geometry immediately.
base <- ggchord(collinear_sequences, validate = "none")
expect_conflict(base + geom_seq(seq_radius = 1) + coord_collinear(), "seq_radius")
expect_conflict(base + coord_collinear() + geom_seq(seq_radius = 1), "seq_radius")
expect_conflict(base + geom_seq(seq_curvature = 0) + coord_collinear(), "seq_curvature")
expect_conflict(base + geom_seq(aes(seq_ring = accver)) + coord_collinear(), "seq_ring")
expect_conflict(base + coord_collinear() +
  scale_seq_ring_manual(values = c(cluster_A = 1)), "seq_ring")
expect_conflict(base + geom_seq_label(seq_label_radius = 1) +
  coord_collinear(), "seq_label_radius")
expect_conflict(base + coord_collinear() + geom_seq_center_label(),
  "coord_collinear")
expect_conflict(base + coord_collinear() +
  geom_feature_label_repel(data = collinear_features, external = FALSE),
  "external = FALSE")

# Explicit gaps above the chord limit are legal linear track spacing.
wide <- base + geom_seq(seq_gap = 0.8) + coord_collinear()
wide_out <- export_ggchord_layout(wide, include = "seq")
stopifnot(all(diff(wide_out$sequence_layout$y) < 0))

# A transformed position scale still locates ticks at their original bp values.
log_plot <- ggchord(data.frame(accver = "log_sequence", length = 1000),
  validate = "none") + geom_seq() +
  scale_seq_position_continuous(transform = "log10",
    breaks = c(1, 10, 100, 1000), minor_breaks = NULL,
    labels = function(x) paste0(x, " bp")) + coord_collinear()
log_out <- export_ggchord_layout(log_plot, include = "axis")
log_ticks <- log_out$axis[!is.na(log_out$axis$label), , drop = FALSE]
log_ticks <- log_ticks[!duplicated(log_ticks$label), , drop = FALSE]
log_ticks <- log_ticks[order(log_ticks$x), , drop = FALSE]
stopifnot(identical(as.character(log_ticks$label),
  c("1 bp", "10 bp", "100 bp", "1000 bp")),
  all(diff(log_ticks$x) > 0),
  diff(tail(log_ticks$x, 2L)) > diff(head(log_ticks$x, 2L)))

# Static label vectors remain paired with breaks on shorter sequences.
static_axis <- export_ggchord_layout(base + geom_seq() +
  scale_seq_position_continuous(transform = "log10",
    breaks = c(1, 10, 100, 1000),
    labels = c("one", "ten", "hundred", "thousand")) + coord_collinear(),
  include = "axis")
stopifnot("one" %in% static_axis$axis$label,
  !any(static_axis$axis$accver == "cluster_B" &
    static_axis$axis$label %in% "thousand"))

# Fitting includes the actual region band and styled backbone vertices.
region_bounds <- export_ggchord_layout(base + geom_seq() +
  geom_seq_region(data = data.frame(accver = "cluster_A", start = 100, end = 200),
    region_offset = 2, region_side = "outside") + coord_collinear(fit = "geometry"),
  include = "feature")
stopifnot(region_bounds$metadata$ylim[2L] >= max(region_bounds$feature$y))
backbone_bounds <- export_ggchord_layout(base +
  geom_seq(seq_style = "band", seq_backbone_width = 2) +
  coord_collinear(fit = "geometry"), include = "seq")
stopifnot(backbone_bounds$metadata$ylim[1L] <= min(backbone_bounds$seq$y),
  backbone_bounds$metadata$ylim[2L] >= max(backbone_bounds$seq$y))

preview_path <- commandArgs(trailingOnly = TRUE)[1L]
if (!is.na(preview_path) && nzchar(preview_path)) {
  ggsave(preview_path, collinear_example, width = 10, height = 7, dpi = 160)
}
cat("Collinear development checks passed.\n")
