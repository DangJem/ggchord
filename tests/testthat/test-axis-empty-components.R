test_that("disabling all position breaks still renders the sequence axis", {
  plot <- ggchord(data.frame(accver = "A", length = 1000),
    validate = "none") + geom_seq() + coord_collinear() +
    scale_seq_position_continuous(breaks = NULL, minor_breaks = NULL)
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
  chord <- ggchord(data.frame(accver = "A", length = 1000)) + geom_seq() +
    scale_seq_position_continuous(breaks = NULL, minor_breaks = NULL)
  expect_equal(nrow(get_chord_layout(chord)$axis_ticks), 0L)
})

test_that("transformed axis breaks retain their genomic anchors", {
  base <- ggchord(data.frame(accver = "A", length = 1000)) + geom_seq()
  plain <- get_chord_layout(base + scale_seq_position_continuous(
    breaks = c(10, 100), minor_breaks = NULL))$axis_ticks
  logged <- get_chord_layout(base + scale_seq_position_continuous(
    breaks = c(10, 100), minor_breaks = NULL, transform = "log10"))$axis_ticks
  expect_equal(nrow(logged), 2L)
  expect_equal(as.character(logged$label), c("10", "100"))
  expect_equal(logged[c("x0", "y0")], plain[c("x0", "y0")])
})

test_that("automatic tick marks draw behind feature arrows", {
  plot <- ggchord(data.frame(accver = "format", length = 5000),
    feature_data = data.frame(accver = "format", start = 500, end = 2000,
      strand = "+", anno = "feature")) +
    geom_seq() + geom_feature() + coord_circular()
  layers <- ggplot2::ggplot_build(plot)$plot$layers
  marks <- which(vapply(layers, function(layer)
    inherits(layer$geom, "GeomChordAxis") &&
      any(layer$data$.component %in% c("major_tick", "minor_tick")), logical(1)))
  features <- which(vapply(layers, function(layer)
    identical(layer$ggchord_type, "gene_poly"), logical(1)))
  text <- which(vapply(layers, function(layer)
    inherits(layer$geom, "GeomChordAxis") &&
      any(layer$data$.component %in% "text"), logical(1)))
  expect_lt(max(marks), min(features))
  expect_gt(min(text), max(features))
})

test_that("central text is not treated as a sequence backbone obstacle", {
  # Format control: the feature forces nontrivial automatic text clearance.
  seqs <- data.frame(accver = "format", length = 5000)
  features <- data.frame(accver = "format", start = 1, end = 1000,
    strand = "+", anno = "format feature")
  plot <- ggchord(seqs, feature_data = features) +
    geom_seq(seq_style = "double") +
    geom_feature_plasmid(position = position_plasmid(offset = -.04)) +
    coord_circular() + theme_ggchord_plasmid()
  plain <- diagnose_ggchord(plot, width = 8, height = 8, device = "png")
  central <- diagnose_ggchord(plot + geom_seq_center_label(),
    width = 8, height = 8, device = "png")
  expect_true(any(plain$label_quality$axis_clearance$radial_shift_mm != 0))
  expect_equal(central$label_quality$axis_clearance,
    plain$label_quality$axis_clearance)
})
