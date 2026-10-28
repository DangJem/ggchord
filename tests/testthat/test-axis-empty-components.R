test_that("disabling all position breaks still renders the sequence axis", {
  plot <- ggchord(data.frame(accver = "A", length = 1000),
    validate = "none") + geom_seq() + coord_collinear() +
    scale_seq_position_continuous(breaks = NULL, minor_breaks = NULL)
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
})
