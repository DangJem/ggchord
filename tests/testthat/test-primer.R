test_that("bundled primer catalogue preserves sequence aliases", {
  catalogue <- primer_catalog("all")
  aliases <- primer_catalog("all", aliases = TRUE)
  expect_equal(nrow(catalogue), 137L)
  expect_equal(nrow(aliases), 163L)
  duplicated_sequence <- aliases$sequence == "GGGAAACGCCTGGTATCTTT"
  expect_setequal(aliases$name[duplicated_sequence], c("hU6-F", "pBR322ori-F"))
})

test_that("primer search derives strand and circular origin", {
  primers <- c(test = "ACGTTGCA")
  target <- c(circle = "TTGCAAAAACG")
  hits <- find_primer_bindings(
    target, primers = primers, circular = TRUE,
    min_annealed_bases = 8L
  )
  expect_equal(nrow(hits), 1L)
  expect_true(hits$crosses_origin)
  expect_equal(hits$strand, "+")
  expect_equal(hits$aliases[[1L]], "test")
})

test_that("primer geometry uses the backbone centre and a head-anchored label", {
  sequence <- data.frame(accver = "circle", length = 1000)
  primers <- data.frame(
    accver = "circle", start = 100L, end = 119L,
    strand = "+", name = "test primer"
  )
  layer <- geom_primer(data = primers)
  expect_identical(layer$ggchord_params$feature_role, "primer")
  expect_true(layer$ggchord_params$is_primer)
  expect_identical(layer$ggchord_params$feature_shape, "primer_arc")
  expect_false(isTRUE(layer$ggchord_params$feature_position$ggchord_feature_stack))
  expect_equal(layer$ggchord_params$feature_position$offset, 0)
  expect_identical(layer$aes_params$feature_fill, "#A020F0")
  expect_identical(layer$aes_params$colour, "#A020F0")
  plot <- ggchord(sequence, validate = "none") +
    geom_seq(seq_style = "double") + layer +
    geom_primer_label_repel(data = primers) +
    coord_circular()
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  layout <- export_ggchord_layout(plot)
  registry <- layout$annotation_registry
  expect_true(any(registry$kind == "primer_label"))
  expect_true(all(registry$annotation_class[registry$kind == "primer_label"] ==
    "primer"))
  polygon <- layout$feature[layout$feature$.component == "polygon", ]
  expect_true(all(sqrt(polygon$x^2 + polygon$y^2) > .9875))
  expect_true(all(sqrt(polygon$x^2 + polygon$y^2) < 1.0125))
  label <- layout$labels[
    layout$labels$.component == "text" &
      layout$labels$annotation_class == "primer", , drop = FALSE
  ]
  expect_equal(unique(label$anchor_position), 119)
  segment_root <- layout$labels[
    layout$labels$.component == "segment", c("x", "y"), drop = FALSE
  ][1L, ]
  distance_to_polygon <- sqrt(
    (polygon$x - segment_root$x)^2 + (polygon$y - segment_root$y)^2
  )
  expect_equal(min(distance_to_polygon), 0, tolerance = 1e-8)
})

test_that("reference primer coordinates are one-based inclusive", {
  data(plasmid_example_pSB1C3)
  primers <- find_primer_bindings(plasmid_example_pSB1C3, set = "reference")
  expect_equal(primers$start, c(155L, 1931L))
  expect_equal(primers$end, c(174L, 1950L))
})

test_that("reverse primer labels anchor at the directional head", {
  sequence <- data.frame(accver = "circle", length = 1000)
  primer <- data.frame(
    accver = "circle", start = 700L, end = 719L,
    strand = "-", name = "reverse primer"
  )
  plot <- ggchord(sequence, validate = "none") +
    geom_seq(seq_style = "double") +
    geom_primer(data = primer) +
    geom_primer_label_repel(data = primer) +
    coord_circular()
  label <- export_ggchord_layout(plot)$labels
  label <- label[label$.component == "text" &
    label$annotation_class == "primer", , drop = FALSE]
  expect_equal(unique(label$anchor_position), 700)
})
