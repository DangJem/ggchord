# v0.7.0 lightweight visual-regression tests.
#
# The layout-fingerprint tests are deterministic and run everywhere. PNG
# rendering and pixel-level md5 comparisons are platform/font-dependent and
# opt-in: set GGCHORD_RUN_SLOW_TESTS=1 (see helper-visual.R).

test_that("key plots render to PNG without errors and produce output", {
  skip_unless_slow_tests()
  set.seed(123)
  paths <- render_regression_set()
  for (nm in names(paths)) {
    expect_true(file.exists(paths[[nm]]), info = nm)
    expect_true(file.info(paths[[nm]])$size > 2000, info = nm)
  }
})

test_that("layout fingerprint is stable across repeated builds", {
  p <- example_ggchord_plot()
  f1 <- ggchord_layout_fingerprint(p)
  f2 <- ggchord_layout_fingerprint(p)
  expect_identical(f1, f2)
})

test_that("layout fingerprint survives saveRDS/readRDS round-trip", {
  p <- example_ggchord_plot()
  f1 <- ggchord_layout_fingerprint(p)
  tmp <- tempfile(fileext = ".rds")
  # serializing layer closures over the package namespace emits a benign
  # "package may not be available" note; the round-trip itself is what matters
  suppressWarnings(saveRDS(p, tmp))
  p2 <- suppressWarnings(readRDS(tmp))
  f2 <- ggchord_layout_fingerprint(p2)
  expect_identical(f1, f2)
})

test_that("layout fingerprint is stable across ggsave", {
  p <- example_ggchord_plot()
  f1 <- ggchord_layout_fingerprint(p)
  out <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(out, p, width = 8, height = 8)
  f2 <- ggchord_layout_fingerprint(p)
  expect_identical(f1, f2)
})

test_that("ribbon color schemes produce distinct layout fingerprints", {
  data(seq_data_example)
  data(ribbon_data_example)
  p_pid <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() + geom_ribbon()
  p_qry <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() + geom_ribbon(ribbon_color_scheme = "query")
  p_sbj <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() + geom_ribbon(ribbon_color_scheme = "subject")
  p_sng <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_color_scheme = "single", ribbon_colors = "orange")
  fps <- vapply(list(p_pid, p_qry, p_sbj, p_sng),
                ggchord_layout_fingerprint, character(1))
  expect_equal(length(unique(fps)), 4)
})

# ---------------------------------------------------------------------------
# Opt-in pixel-level regression against committed PNG md5 baselines
# ---------------------------------------------------------------------------

test_that("rendered PNGs match the committed md5 baseline (opt-in)", {
  skip_unless_slow_tests()
  baseline_file <- testthat::test_path("visual-baseline.rds")
  skip_if(!file.exists(baseline_file), "no visual-baseline.rds committed")
  baseline <- readRDS(baseline_file)
  paths <- render_regression_set()
  for (nm in names(paths)) {
    h <- tools::md5sum(paths[[nm]])
    expect_identical(unname(h), unname(baseline[[nm]]), info = nm)
  }
})
