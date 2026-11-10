test_that("feature presets accept fixed shape and preserve explicit scales", {
  features <- transform(smoke_genes()[1, ], type = "CDS")
  p <- ggchord(smoke_sequences()[1, ], feature_data = features) + geom_seq() +
    scale_feature_fill_manual(values = c(CDS = "black")) +
    geom_feature_plasmid(shape = "block") + coord_circular()
  built <- build_ggchord_smoke(p)
  expect_true(all(built$data[[2]]$feature_fill == "black"))
})

test_that("unified annotations select genes and preserve independent feature layers", {
  annotations <- transform(smoke_genes(), type = c("CDS", "promoter"))
  p <- ggchord(smoke_sequences(), feature_data = annotations) + geom_seq() +
    geom_gene() + geom_feature(data = function(x) x[x$type != "CDS", ], shape = "block")
  out <- export_ggchord_layout(p, include = c("gene", "feature"))
  expect_setequal(unique(out$gene$type), "CDS")
  expect_setequal(unique(out$feature$type), "promoter")
  expect_error(ggchord(smoke_sequences(), gene_data = annotations,
    feature_data = annotations), "only one")
})

test_that("feature factor coordinates retain genomic labels rather than level codes", {
  features <- data.frame(accver = "A",
    start = factor("100", levels = c("100", "300")),
    end = factor("300", levels = c("100", "300")), strand = "+", type = "CDS",
    source_row = 999L)
  p <- ggchord(smoke_sequences()[1, ]) + geom_seq() +
    geom_feature(data = features, shape = "block") + coord_collinear()
  out <- export_ggchord_layout(p, include = "feature")
  expect_equal(unique(out$feature$start), 100)
  expect_equal(unique(out$feature$end), 300)
  expect_setequal(out$feature$source_row, 1L)
  expect_equal(as.character(features$start), "100")
  expect_equal(as.character(features$end), "300")
})
