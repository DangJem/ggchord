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
