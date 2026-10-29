test_that("named tracks share domains and compose with selected axes and a bp reference", {
  data <- data.frame(accver = c("A", "B"), start = 100, end = 200, value = c(-1, 2))
  plot <- ggchord(smoke_sequences()) + geom_seq() +
    geom_track_bar(aes(track_value = value), data = data,
      position = position_track("signal")) +
    scale_track_value_continuous("signal", limits = c(-2, 3)) +
    scale_seq_position_continuous(sequences = "B") + geom_seq_scale(length = 100) +
    coord_collinear()
  out <- export_ggchord_layout(plot, include = c("track", "axis"))
  expect_equal(out$track_registry$track_id, c("signal", "signal"))
  expect_true(all(vapply(out$track_registry$limits, function(x) all(x == c(-2, 3)), logical(1))))
  expect_true(any(out$axis$component == "seq_scale"))
  expect_s3_class(build_ggchord_smoke(plot), "ggplot_built")
})
