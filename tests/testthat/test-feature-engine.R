test_that("feature presets accept the shape alias and preserve explicit scales", {
  features <- transform(smoke_genes()[1, ], type = "CDS")
  p <- ggchord(smoke_sequences()[1, ]) + geom_seq() +
    scale_feature_fill_manual(values = c(CDS = "black")) +
    geom_feature_plasmid(data = features, shape = "block") + coord_circular()
  built <- build_ggchord_smoke(p)
  expect_true(all(built$data[[2]]$feature_fill == "black"))
})

test_that("feature factor coordinates retain genomic labels rather than level codes", {
  features <- data.frame(accver = "A",
    start = factor("100", levels = c("100", "300")),
    end = factor("300", levels = c("100", "300")), strand = "+", type = "CDS")
  p <- ggchord(smoke_sequences()[1, ]) + geom_seq() +
    geom_feature(data = features, feature_shape = "block") + coord_collinear()
  out <- export_ggchord_layout(p, include = "feature")
  expect_equal(unique(out$feature$start), 100)
  expect_equal(unique(out$feature$end), 300)
  expect_equal(as.character(features$start), "100")
  expect_equal(as.character(features$end), "300")
})
