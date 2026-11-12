test_that("an optional layout policy preserves data and explicit guide placement", {
  seqs <- data.frame(accver = c("A", "B"), length = c(100, 200))
  base <- ggchord(seqs) + geom_seq()
  p <- base +
    guides(seq_colour = ggplot2::guide_legend(position = "bottom", nrow = 1)) +
    adapt_ggchord()
  out <- diagnose_ggchord(p, width = 8, height = 6, device = "png")
  expect_identical(out$guide_layout$status, "feasible")
  expect_equal(out$guide_layout$guides$final_position, "bottom")
  expect_equal(out$guide_layout$guides$n_items, 2L)
  expect_true(out$label_quality$summary$measured)
  expect_equal(out$label_quality$summary$n_text_over_feature, 0L)
  expect_equal(export_ggchord_layout(p, original_data = TRUE)$original_data,
    export_ggchord_layout(base, original_data = TRUE)$original_data)
})

test_that("a circular coordinate string does not report its own glyphs as collisions", {
  plot <- ggchord(data.frame(accver = "circle", length = 1000)) +
    geom_seq() + coord_circular() + theme_ggchord_plasmid() + adapt_ggchord()
  out <- diagnose_ggchord(plot, width = 6, height = 9, device = "png")
  expect_equal(out$label_quality$summary$n_text_overlaps, 0L)
})
