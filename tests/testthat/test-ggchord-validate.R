# v0.7.0 tests: ggchord(validate = "warn" | "error" | "none") integration

test_that("default validate = 'warn' keeps valid input warning-free", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  expect_no_warning(
    p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example))
  expect_s3_class(p, "ggchord")
  expect_true(p$ggchord$validation$valid)
})

test_that("validate = 'warn' emits a single summary warning mentioning the problem", {
  data(seq_data_example)
  bad <- data.frame(
    qaccver = "MT108731.1", saccver = "MT118296.1",
    length = 100, pident = 90,
    qstart = 50000, qend = 40000, sstart = 1, send = 100
  )
  expect_warning(ggchord(seq_data_example, bad), "start > end")
})

test_that("validate = 'error' stops on severe problems", {
  data(seq_data_example)
  bad <- data.frame(
    qaccver = "MT108731.1", saccver = "NOT_A_SEQUENCE",
    length = 100, pident = 90,
    qstart = 1, qend = 100, sstart = 1, send = 100
  )
  expect_error(ggchord(seq_data_example, bad, validate = "error"),
               "failed validation")
})

test_that("validate = 'none' skips diagnostics without crashing", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, pident = pident + 200)
  expect_no_warning(p <- ggchord(seq_data_example, bad, validate = "none"))
  expect_null(p$ggchord$validation)
  expect_no_error(ggplot_build(p))
})

test_that("validation result is cached on the plot object", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, saccver = "UNKNOWN")
  p <- suppressWarnings(ggchord(seq_data_example, bad))
  expect_s3_class(p$ggchord$validation, "ggchord_validation")
  expect_false(p$ggchord$validation$valid)
})

test_that("structural errors are still hard errors in every mode", {
  data(seq_data_example)
  expect_error(ggchord(data.frame(id = "a", len = 1)), "seq_data")
  expect_error(ggchord(transform(seq_data_example, length = NA_real_)),
               "finite positive")
  expect_error(
    ggchord(data.frame(seq_id = c("a", "a"), length = c(1, 2))), "unique")
})

test_that("valid input produces identical layout with warn and none", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p1 <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene()
  p2 <- ggchord(seq_data_example, ribbon_data_example, gene_data_example,
                validate = "none") +
    geom_seq() + geom_ribbon() + geom_gene()
  ggplot_build(p1)
  l1 <- get_chord_layout()
  ggplot_build(p2)
  l2 <- get_chord_layout()
  expect_identical(l1$ribbon_polys, l2$ribbon_polys)
  expect_identical(l1$seq_arcs, l2$seq_arcs)
})
