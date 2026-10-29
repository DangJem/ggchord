test_that("feature presets accept the shape alias and preserve explicit scales", {
  features <- transform(smoke_genes()[1, ], type = "CDS")
  p <- ggchord(smoke_sequences()[1, ]) + geom_seq() +
    scale_feature_fill_manual(values = c(CDS = "black")) +
    geom_feature_plasmid(data = features, shape = "block") + coord_circular()
  built <- build_ggchord_smoke(p)
  expect_true(all(built$data[[2]]$feature_fill == "black"))
})
