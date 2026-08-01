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

test_that("ribbon outline parameters work with sensible defaults", {
  data(seq_data_example)
  data(ribbon_data_example)
  # Defaults: black outline, width 0.05, solid line
  p <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
  l <- p$layers[[2]]
  expect_equal(l$aes_params$colour, "black")
  expect_equal(l$aes_params$linewidth, 0.05)
  expect_equal(l$aes_params$linetype, 1)

  # Custom values
  p2 <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_outline_color = "red", ribbon_outline_width = 0.8,
                ribbon_outline_linetype = "dashed")
  l2 <- p2$layers[[2]]
  expect_equal(l2$aes_params$colour, "red")
  expect_equal(l2$aes_params$linewidth, 0.8)
  expect_equal(l2$aes_params$linetype, "dashed")

  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p2)))
  dev.off()
})

test_that("plot objects are self-contained (no cross-talk between plots)", {
  data(seq_data_example)
  p1 <- ggchord(seq_data_example) + geom_seq(seq_radius = 5)
  p2 <- ggchord(seq_data_example) + geom_seq()
  expect_equal(p1$layers[[1]]$ggchord_params$seq_radius, 5)
  expect_null(p2$layers[[1]]$ggchord_params$seq_radius)
  # both build independently
  expect_s3_class(ggplot_build(p1), "ggplot_built")
  expect_s3_class(ggplot_build(p2), "ggplot_built")
})

test_that("plots survive saveRDS/readRDS and render", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq(seq_radius = 5)
  f <- tempfile(fileext = ".rds")
  saveRDS(p, f)
  p2 <- readRDS(f)
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p2)))
  dev.off()
  expect_equal(p2$layers[[1]]$ggchord_params$seq_radius, 5)
})

test_that("ggsave works directly on a ggchord plot", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
  f <- tempfile(fileext = ".png")
  expect_no_error(ggsave(f, p, width = 6, height = 6, dpi = 72))
  expect_true(file.exists(f))
})

test_that("plot objects are self-contained and repeated builds stay stable", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis()
  # The plot carries its own scales (tagged ggchord_managed) so that tools such
  # as plotly::ggplotly() that clone the plot see the same scales as a build.
  expect_gt(length(p$scales$scales), 0)
  expect_true(all(vapply(p$scales$scales,
                         function(s) !is.null(s$ggchord_managed), logical(1))))
  expect_false(is.null(p$layers[[1]]$ggchord_params))
  # repeated builds remain stable
  expect_no_error(ggplot_build(p))
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
})

test_that("ribbon data sanity checks warn on bad input", {
  data(seq_data_example)
  bad_ribbon <- data.frame(
    qaccver = c("MT108731.1", "MT108731.1"),
    saccver = c("MT118296.1", "MT118296.1"),
    length = c(100, 100), pident = c(90, 90),
    qstart = c(50000, 100), qend = c(40000, 200),   # qstart > qend
    sstart = c(1, 1), send = c(100, 100)
  )
  expect_warning(ggchord(seq_data_example, bad_ribbon), "start > end")
})

test_that("geom_seq_label places sequence labels", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) +
    geom_seq() +
    geom_seq_label(seq_label_radius = 1.25, seq_label_size = 3.5)
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
  layout <- ggchord:::get_chord_layout()
  expect_gt(nrow(layout$seq_labels_df), 0)
  expect_true(all(c("text_x", "text_y", "label") %in% names(layout$seq_labels_df)))
})

test_that("ribbon subject color scheme works", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_color_scheme = "subject")
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
})

test_that("theme customization via + works", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() +
    ggplot2::theme(legend.position = "bottom")
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
})

test_that("get_chord_layout is exported after building", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq()
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  print(p)
  dev.off()
  layout <- get_chord_layout()
  expect_true(inherits(layout, "chord_layout"))
  expect_gt(length(layout$seq_arcs), 0)
})

test_that("plots convert to plotly when plotly is installed", {
  skip_if_not_installed("plotly")
  suppressWarnings(library(plotly))
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  # sequence-only plot (the common interactive use case)
  p <- ggchord(seq_data_example) + geom_seq()
  pl <- suppressWarnings(plotly::ggplotly(p))
  expect_true(inherits(pl, "plotly"))
  expect_gt(length(pl$x$data), 0)
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
