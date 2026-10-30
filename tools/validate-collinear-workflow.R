# Targeted checks for the expanded v0.14.0 workflow. Run explicitly:
# Rscript tools/validate-collinear-workflow.R
# No testthat, package checks or release device matrix is run here.
pkgload::load_all(quiet = TRUE)
source(file.path("examples", "collinear-workflow.R"), local = TRUE)

expect_conflict <- function(expr, text) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"),
    grepl(text, conditionMessage(error), fixed = TRUE))
}
seqs <- workflow_sequences
links <- workflow_links
base <- ggchord(seqs, links, validate = "none") + geom_seq()
default_plot <- base + geom_link_ribbon() + coord_collinear()
default <- export_ggchord_layout(default_plot)
stopifnot(setequal(unique(default$ribbon$source_row), c(1L, 3L)),
  identical(default$connection_report$retained, c(TRUE, FALSE, TRUE)),
  default$connection_report$reason[2L] == "non_adjacent")
drawn <- ggplot2::ggplot_build(default_plot)$data[[2L]]
stopifnot(all(drawn$linetype == "solid"))
automatic_all <- base + geom_link_ribbon(adjacent_only = FALSE) +
  coord_collinear()
automatic_lines <- unique(export_ggchord_layout(automatic_all)$ribbon[
  c("source_row", "linetype_val")])
stopifnot(identical(as.character(automatic_lines$linetype_val[
  match(seq_len(nrow(links)), automatic_lines$source_row)]),
  c("solid", "dashed", "solid")))
all_plot <- base + geom_link_ribbon(adjacent_only = FALSE,
  linetype = "solid") + coord_collinear()
all_links <- export_ggchord_layout(all_plot)
stopifnot(setequal(unique(all_links$ribbon$source_row), seq_len(nrow(links))),
  all(ggplot2::ggplot_build(all_plot)$data[[2L]]$linetype == "solid"))
selected <- export_ggchord_layout(base + geom_link_ribbon(
  pairs = data.frame(qaccver = "C", saccver = "A")) + coord_collinear())
stopifnot(identical(unique(selected$ribbon$source_row), 2L),
  sum(selected$connection_report$retained) == 1L)
expect_conflict(export_ggchord_layout(base + geom_link_ribbon(
  adjacent_only = TRUE, pairs = data.frame(qaccver = "A", saccver = "C")) +
    coord_collinear()), "non-adjacent pairs")

# Reversible anchor translations and supported display directions.
aligned <- export_ggchord_layout(workflow_aligned)
anchor_mid <- rowMeans(workflow_features[c("start", "end")])
mapped <- map_ggchord_position(workflow_aligned, seqs$accver, anchor_mid)
stopifnot(all(abs(mapped$x) < 1e-8),
  identical(aligned$sequence_layout$orientation, c(1, -1, 1)),
  all(aligned$orientation_diagnostics$status[-1L] == "supported"),
  all(aligned$label_diagnostics$status == "placed"))
expect_conflict(ggchord(seqs, links, validate = "none") +
  geom_seq(seq_orientation = -1) + coord_collinear(orientation = "auto"),
  "cannot be combined")
expect_conflict(coord_collinear(min_support = 2), "min_support")
unknown_links <- links
unknown_links$orientation_known <- FALSE
unresolved <- export_ggchord_layout(ggchord(seqs, unknown_links,
  validate = "none") + geom_seq() + geom_link_ribbon() +
    coord_collinear(orientation = "auto"))
stopifnot(all(unresolved$sequence_layout$orientation == 1),
  all(unresolved$orientation_diagnostics$status[-1L] == "no_evidence"))

# A narrow, empty channel stays usable; impossible explicit clearances fail.
narrow_seqs <- seqs[1:2, ]
narrow_links <- links[1L, ]
narrow <- ggchord(narrow_seqs, narrow_links, validate = "none") +
  geom_seq(seq_gap = .01) + geom_link_ribbon() + coord_collinear()
channel <- export_ggchord_layout(narrow)
stopifnot(all(channel$ribbon$y >= min(channel$sequence_layout$y)),
  all(channel$ribbon$y <= max(channel$sequence_layout$y)))
expect_conflict(export_ggchord_layout(ggchord(narrow_seqs, narrow_links,
  validate = "none") + geom_seq(seq_gap = .01) +
    geom_link_ribbon(link_gap = .03) + coord_collinear()), "do not fit")
expect_conflict(export_ggchord_layout(ggchord(seqs, validate = "none") +
  geom_seq(seq_order = c("A", "A", "B")) + coord_collinear()), "non-repeated")

# Tiny arrows retain genomic extents and may share a lane without overlap.
tiny <- data.frame(accver = "long", start = c(100, 110), end = c(102, 112),
  strand = "+", type = "CDS", label = c("first", "second"))
tiny_plot <- ggchord(data.frame(accver = "long", length = 1e6),
  validate = "none") + geom_seq() +
  geom_feature(data = tiny, position = position_feature_stack()) + coord_collinear()
tiny_out <- export_ggchord_layout(tiny_plot)$feature
parts <- split(tiny_out, tiny_out$source_row)
stopifnot(!any(tiny_out$display_extended), length(unique(tiny_out$lane)) == 1L,
  max(parts[[1L]]$x) < min(parts[[2L]]$x))

# A reordered/subset label layer inherits the full source geometry's lanes.
features <- data.frame(feature_id = c("a", "b", "c"), accver = "A",
  start = c(100, 200, 600), end = c(400, 500, 700), strand = "+",
  type = "CDS", label = c("alpha", "beta", "gamma"))
label_plot <- ggchord(seqs, validate = "none") + geom_seq() +
  geom_feature(data = features, position = position_feature_stack(
    spacing = .2, base_position = position_strand(.2))) +
  geom_feature_label_repel(data = features[c(2L, 1L), ]) + coord_collinear()
label_out <- export_ggchord_layout(label_plot)
geometry <- unique(label_out$feature[c("source_row", "lane", "normal_offset")])
labels <- unique(label_out$labels[c("source_row", "lane", "normal_offset")])
index <- match(c(2L, 1L), geometry$source_row)
labels <- labels[match(1:2, labels$source_row), ]
stopifnot(identical(labels$lane, geometry$lane[index]),
  isTRUE(all.equal(labels$normal_offset, geometry$normal_offset[index])))

# Multiple matching entity layers require an explicit label Position.
one <- features[1L, ]
ambiguous <- ggchord(seqs, validate = "none") + geom_seq() +
  geom_feature(data = one, position = position_strand(.2)) +
  geom_feature(data = one, position = position_strand(-.2)) + coord_collinear()
expect_conflict(export_ggchord_layout(ambiguous +
  geom_feature_label_repel(data = one)), "Multiple geometry layers")
explicit <- export_ggchord_layout(ambiguous +
  geom_feature_label_repel(data = one, position = position_strand(.3)))
stopifnot(all(explicit$labels$normal_offset == .3))

# Cropped coordinates preserve source ticks, end marks and an inverse mapping.
source_out <- export_ggchord_layout(workflow_source_axis)
stopifnot(any(source_out$seq$component == "seq_break"),
  all(source_out$sequence_layout$truncated_left),
  any(as.numeric(source_out$axis$label) >
    min(workflow_focused$seq_data$.source_origin), na.rm = TRUE))
forward <- map_ggchord_position(workflow_source_axis, c("A", "B"),
  c(400, 400), source = TRUE)
inverse <- map_ggchord_position(workflow_source_axis, forward$accver,
  forward$x, source = TRUE, inverse = TRUE)
stopifnot(isTRUE(all.equal(inverse$position, c(400, 400))))
cat("Expanded collinear workflow checks passed.\n")
