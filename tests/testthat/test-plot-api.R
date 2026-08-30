build_ggchord_smoke <- function(plot) {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  ggplot2::ggplot_build(plot)
}

test_that("ggchord creates a plot and rejects invalid sequence data", {
  data(seq_data_example)

  p <- ggchord(seq_data_example, validate = "none")
  expect_s3_class(p, "ggchord")
  expect_s3_class(coord_chord(), "Coord")

  expect_error(
    ggchord(data.frame(id = "A", length = 100)),
    "seq_data"
  )
})

test_that("the main plotting layers build together", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p <- ggchord(
    seq_data_example, ribbon_data_example, gene_data_example,
    validate = "none"
  ) +
    geom_seq() +
    geom_ribbon() +
    geom_gene() +
    geom_gene_label() +
    geom_seq_label() +
    geom_axis()

  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout()
  expect_true(length(layout$seq_arcs) > 0)
  expect_true(nrow(layout$ribbon_polys) > 0)
  expect_true(nrow(layout$gene_polys) > 0)
})

test_that("all automatic gene-label layouts build", {
  data(seq_data_example)
  data(gene_data_example)

  for (mode in c("aligned", "radial", "arc")) {
    p <- ggchord(seq_data_example, gene_data = gene_data_example,
                 validate = "none") +
      geom_seq() +
      geom_gene() +
      geom_gene_label_repel(gene_label_layout = mode)
    expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  }

  expect_error(
    geom_gene_label_repel(force = 1),
    "Removed"
  )
})

test_that("region and ribbon-highlight layers build", {
  data(seq_data_example)
  data(ribbon_data_example)

  regions <- data.frame(
    seq_id = seq_data_example$seq_id[1],
    start = 100,
    end = 500
  )
  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() +
    geom_ribbon() +
    geom_seq_region(regions = regions) +
    geom_ribbon_highlight(ribbon_ids = 1)

  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout()
  expect_true(nrow(layout$region_polys) > 0)
  expect_true(nrow(layout$ribbon_highlight_polys) > 0)
})

test_that("feature and sequence-group layers build", {
  data(seq_data_example)

  feature <- data.frame(
    seq_id = seq_data_example$seq_id[1],
    start = 100,
    end = 500,
    strand = "+",
    type = "CDS"
  )
  groups <- stats::setNames(
    rep(c("group-a", "group-b"), length.out = nrow(seq_data_example)),
    seq_data_example$seq_id
  )
  p <- ggchord(seq_data_example, validate = "none") +
    geom_seq(seq_group = groups) +
    geom_feature(feature)

  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout()
  expect_true(nrow(layout$gene_polys) > 0)
  expect_true(nrow(layout$group_labels) > 0)
})
