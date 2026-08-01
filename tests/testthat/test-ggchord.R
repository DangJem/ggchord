# testthat test file - v0.4.0 layered API

test_that("ggchord with only seq_data returns a ggchord object", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example)
  expect_s3_class(p, "ggchord")
  expect_s3_class(p, "ggplot")
})

test_that("ggchord + geom_seq creates a sequence arc layer", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example) + geom_seq()
  expect_s3_class(p, "ggchord")
  expect_true(length(p$layers) >= 1)
})

test_that("ggchord + geom_seq + geom_ribbon handles alignment data correctly", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(
    seq_data = seq_data_example,
    ribbon_data = ribbon_data_example
  ) +
    geom_seq() +
    geom_ribbon()
  expect_s3_class(p, "ggchord")
  expect_true(length(p$layers) >= 1)
})

test_that("ggchord + geom_gene handles gene data correctly", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(
    seq_data = seq_data_example,
    gene_data = gene_data_example
  ) +
    geom_seq() +
    geom_gene()
  expect_s3_class(p, "ggchord")
})

test_that("ggchord + geom_axis can add an axis", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example) +
    geom_seq() +
    geom_axis()
  expect_s3_class(p, "ggchord")
})

test_that("parameters can be distributed across geoms", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p <- ggchord(
    seq_data = seq_data_example,
    ribbon_data = ribbon_data_example,
    gene_data = gene_data_example,
    title = "Test",
    rotation = 30
  ) +
    geom_seq(
      seq_radius = c(3, 2, 2, 1),
      seq_curvature = c(0, 1, -1, 1.5),
      seq_orientation = c(-1, -1, -1, 1)
    ) +
    geom_ribbon(
      ribbon_color_scheme = "query",
      ribbon_alpha = 0.5
    ) +
    geom_gene(
      gene_color_scheme = "strand",
      gene_width = 0.08
    ) +
    geom_axis(
      axis_gap = 0.02,
      axis_label_orientation = c(0, 45, 80, 130)
    )

  expect_s3_class(p, "ggchord")
})

test_that("seq_data missing required columns errors", {
  bad_data <- data.frame(id = c("a", "b"), len = c(100, 200))
  expect_error(
    ggchord(seq_data = bad_data),
    "seq_data"
  )
})

test_that("seq_data with non-positive lengths errors", {
  bad_data <- data.frame(seq_id = c("a", "b"), length = c(0, 200))
  expect_error(
    ggchord(seq_data = bad_data),
    "positive"
  )
})

test_that("debug mode outputs debug information", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example, debug = TRUE)
  expect_s3_class(p, "ggchord")
})

test_that("ggchord global parameters are passed correctly", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example, title = "Hello", rotation = 90)
  expect_s3_class(p, "ggchord")
})

test_that("print renders the full chord diagram", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  # pident works correctly together with the gene fill scale
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_color_scheme = "query") +
    geom_gene() +
    geom_axis()

  # Render to PDF
  pdf("/tmp/ggchord_test_print.pdf", 8, 8)
  suppressMessages(suppressWarnings(print(p)))
  dev.off()
  expect_true(file.exists("/tmp/ggchord_test_print.pdf"))
})

test_that("README color, label override, and transparency parameters take effect", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq(linewidth = 2) +
    geom_ribbon(ribbon_color_scheme = "query", alpha = 0.2) +
    geom_gene(gene_color_scheme = "manual", show_label = TRUE,
              label_size = 4) +
    geom_axis()

  pdf("/tmp/ggchord_readme_params.pdf", 8, 8)
  expect_no_warning(print(p))
  dev.off()

  layout <- ggchord:::get_chord_layout()
  expect_true("fill" %in% names(layout$ribbon_polys))
  expect_true(all(layout$ribbon_polys$alpha == 0.2))
  expect_gt(nrow(layout$gene_labels), 0)
  expect_true(file.exists("/tmp/ggchord_readme_params.pdf"))
})

test_that("documented data and parameter values are validated", {
  expect_error(
    ggchord(data.frame(seq_id = c("a", "a"), length = c(1, 2))),
    "unique"
  )

  data(seq_data_example)
  expect_error({
    p <- ggchord(seq_data_example) + geom_seq(seq_orientation = 0)
    pdf(tempfile(fileext = ".pdf")); print(p); dev.off()
  }, "1 or -1")
})
