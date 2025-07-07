test_that("ggchord can generate a ggplot object", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example)
  expect_s3_class(p, "ggplot")  # Verify the return value is a ggplot object
})

# test_that("Gene annotation arrows are drawn correctly", {
#   data(seq_data_example, gene_data_example)
#   p <- ggchord(seq_data = seq_data_example, gene_data = gene_data_example)
#   expect_true("geom_polygon" %in% sapply(p$layers, function(x) class(x$geom)[1]))  # Verify that it contains polygon layers (gene arrows)
# })
