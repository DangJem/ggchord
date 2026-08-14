test_that("ribbon_color_by maps an arbitrary numeric column", {
  data(seq_data_example)
  data(ribbon_data_example)
  rb <- transform(ribbon_data_example,
                  bitscore = seq_len(nrow(ribbon_data_example)) * 10)

  p <- ggchord(seq_data_example, rb, validate = "none") +
    geom_seq() +
    geom_ribbon(ribbon_color_by = "bitscore", ribbon_color_name = "Bit score")
  ggplot_build(p)
  layout <- get_chord_layout()

  expect_equal(layout$ribbon_color_scheme, "value")
  expect_equal(layout$ribbon_color_name, "Bit score")
  expect_true("value" %in% names(layout$ribbon_polys))
  expect_true(all(is.finite(layout$ribbon_polys$value)))
})

test_that("ribbon alpha / outline / linetype mappings work", {
  data(seq_data_example)
  data(ribbon_data_example)
  rb <- transform(ribbon_data_example,
                  score = seq_len(nrow(ribbon_data_example)),
                  kind = rep(c("a", "b"), length.out = nrow(ribbon_data_example)))

  p <- ggchord(seq_data_example, rb, validate = "none") +
    geom_seq() +
    geom_ribbon(
      ribbon_alpha_by = "score",
      ribbon_outline_by = "kind",
      ribbon_outline_colors = c(a = "red", b = "blue"),
      ribbon_linetype_by = "kind",
      ribbon_linetypes = c(a = "solid", b = "dashed")
    )
  ggplot_build(p)
  layout <- get_chord_layout()

  expect_gt(length(unique(layout$ribbon_polys$alpha)), 1)
  expect_equal(sort(unique(layout$ribbon_polys$outline_col)), c("blue", "red"))
  expect_equal(sort(unique(layout$ribbon_polys$linetype_val)), c("dashed", "solid"))
})

test_that("ribbon_direction distinguishes same and reverse alignments", {
  data(seq_data_example)
  data(ribbon_data_example)
  rb <- ribbon_data_example
  rb$sstart[1] <- rb$send[1]
  rb$send[1] <- 1

  p <- ggchord(seq_data_example, rb, validate = "none") +
    geom_seq() +
    geom_ribbon(ribbon_direction = "alpha")
  ggplot_build(p)
  layout <- get_chord_layout()
  expect_equal(sort(unique(layout$ribbon_polys$direction)), c("reverse", "same"))
  expect_gt(length(unique(layout$ribbon_polys$alpha)), 1)
})

test_that("geom_seq_region draws bands and respects empty data", {
  data(seq_data_example)
  data(ribbon_data_example)
  regions <- data.frame(
    seq_id = c("MT108731.1", "MT118296.1"),
    start = c(100, 200), end = c(500, 900),
    label = c("A", "B"), color = c("orange", "purple"),
    stringsAsFactors = FALSE
  )

  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() + geom_ribbon() + geom_seq_region(regions = regions)
  ggplot_build(p)
  layout <- get_chord_layout()
  expect_true(nrow(layout$region_polys) > 0)
  expect_equal(sort(unique(layout$region_polys$zregionfill)), c("orange", "purple"))

  p_empty <- ggchord(seq_data_example, validate = "none") +
    geom_seq() + geom_seq_region(regions = regions[0, ])
  expect_no_error(ggplot_build(p_empty))
  expect_equal(nrow(get_chord_layout()$region_polys), 0)
})

test_that("geom_ribbon_highlight selects safely and returns an empty layer when nothing matches", {
  data(seq_data_example)
  data(ribbon_data_example)

  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() + geom_ribbon() + geom_ribbon_highlight(ribbon_ids = c(1, 3))
  ggplot_build(p)
  layout <- get_chord_layout()
  expect_equal(unique(layout$ribbon_highlight_polys$source_row), c(1, 3))

  p2 <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() + geom_ribbon() +
    geom_ribbon_highlight(predicate = function(d) d$pident > 999)
  expect_no_error(ggplot_build(p2))
  expect_equal(nrow(get_chord_layout()$ribbon_highlight_polys), 0)

  expect_error(
    ggplot_build(ggchord(seq_data_example, ribbon_data_example, validate = "none") +
      geom_seq() + geom_ribbon() + geom_ribbon_highlight(predicate = "pident > 90")),
    "predicate"
  )
})

test_that("geom_feature reuses gene geometry while geom_gene stays available", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  f <- data.frame(
    seq_id = c("MT108731.1", "MT118296.1"),
    start = c(100, 50), end = c(300, 200),
    strand = c("+", "-"), type = c("CDS", "tRNA"),
    stringsAsFactors = FALSE
  )
  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() + geom_ribbon() + geom_feature(f)
  ggplot_build(p)
  expect_true(nrow(get_chord_layout()$gene_polys) > 0)
  expect_equal(get_chord_layout()$gene_color_scheme, "manual")

  p_gene <- ggchord(seq_data_example, gene_data = gene_data_example, validate = "none") +
    geom_seq() + geom_gene()
  ggplot_build(p_gene)
  expect_true(nrow(get_chord_layout()$gene_polys) > 0)
})
