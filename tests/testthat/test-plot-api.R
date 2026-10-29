test_that("constructor and core sequence, ribbon and gene layers build", {
  p <- ggchord(smoke_sequences(), smoke_ribbon(), smoke_genes(), validate = "none") +
    geom_seq() + geom_link_ribbon() + geom_gene() + geom_seq_label()
  expect_s3_class(p, "ggchord")
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
})

test_that("the three label strategies retain representative annotations", {
  base <- ggchord(smoke_sequences(), gene_data = smoke_genes()) +
    geom_seq() + geom_gene()
  for (mode in c("radial", "auto", "arc")) {
    layout <- get_chord_layout(base +
      geom_gene_label_repel(gene_label_layout = mode))
    expect_setequal(layout$gene_labels$label, smoke_genes()$anno)
  }
  shrink <- get_chord_layout(base + geom_gene_label_repel(
    gene_label_fit = "shrink", gene_label_min_size = 1.5))
  expect_true(all(shrink$gene_labels$label_size_final >= 1.5))
})

test_that("region, feature and point-link layers retain their source rows", {
  features <- smoke_genes()
  points <- data.frame(qaccver = "A", saccver = "B", qpos = 500, spos = 600)
  p <- ggchord(smoke_sequences()) + geom_seq() +
    geom_feature(data = features, feature_shape = "block") +
    geom_seq_region(data = features[1, c("accver", "start", "end")]) +
    geom_link_line(data = points)
  out <- export_ggchord_layout(p, include = c("feature", "link"))
  expect_setequal(out$feature$source_row, 1:2)
  expect_true(any(out$feature$component == "seq_region"))
  expect_setequal(out$link$source_row, 1L)
})

test_that("independent ribbon layers retain their data and manual styles", {
  first <- smoke_ribbon()
  second <- transform(first, qstart = 500, qend = 600)
  p <- ggchord(smoke_sequences()) + geom_seq() +
    geom_link_ribbon(data = first, fill = "orange") +
    geom_link_ribbon(data = function(x) second, fill = "navy")
  built <- build_ggchord_smoke(p)
  expect_true(all(built$data[[2]]$ribbon_fill == "orange"))
  expect_true(all(built$data[[3]]$ribbon_fill == "navy"))
  out <- get_chord_layout(p)$layer_inputs
  inputs <- lapply(Filter(function(x) !is.null(x$ribbon), out), function(x) x$ribbon)
  expect_setequal(vapply(inputs, function(x) x$qstart[1], numeric(1)), c(100, 500))
})
