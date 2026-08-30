test_that("label box separation stops after convergence", {
  # These overlapping boxes cannot move because both axes are fixed. The old
  # loop still ran all 500 iterations because it only tracked collision hits,
  # not whether the coordinates actually changed after boundary clamping.
  result <- ggchord:::ggchord_separate_boxes(
    x = c(0, 0), y = c(0, 0),
    bw = c(1, 1), bh = c(1, 1),
    cx_off = c(0, 0), cy_off = c(0, 0),
    max_iter = 500,
    x_lim = c(0, 0), y_lim = c(0, 0)
  )

  expect_equal(result$iterations, 1L)
  expect_equal(result$x, c(0, 0))
  expect_equal(result$y, c(0, 0))
})

test_that("vectorised label box separation removes collisions", {
  n <- 40
  pair_centres <- rep(seq(0, 9.5, by = 0.5), each = 2)
  result <- ggchord:::ggchord_separate_boxes(
    x = pair_centres + rep(c(0, 0.02), n / 2),
    y = rep(0, n),
    bw = rep(0.08, n), bh = rep(0.08, n),
    cx_off = rep(0, n), cy_off = rep(0, n),
    max_iter = 500
  )

  dx <- abs(outer(result$x, result$x, "-"))
  dy <- abs(outer(result$y, result$y, "-"))
  overlap <- upper.tri(dx) & dx < 0.08 - 1e-7 & dy < 0.08 - 1e-7
  expect_false(any(overlap))
})
