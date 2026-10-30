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
      geom_gene_label_repel(label_layout = mode))
    expect_setequal(layout$gene_labels$label, smoke_genes()$anno)
  }
  shrink <- get_chord_layout(base + geom_gene_label_repel(
    label_fit = "shrink", label_min_size = 1.5))
  expect_true(all(shrink$gene_labels$label_size_final >= 1.5))
})

test_that("region, feature and point-link layers retain their source rows", {
  features <- smoke_genes()
  regions <- data.frame(accver = c("missing", "A"),
    start = factor(c("1", "100"), levels = c("1", "100", "300")),
    end = factor(c("100", "300"), levels = c("1", "100", "300")))
  points <- data.frame(qaccver = "A", saccver = "B", qpos = 500, spos = 600)
  p <- ggchord(smoke_sequences()) + geom_seq() +
    geom_feature(data = features, shape = "block") +
    geom_seq_region(data = regions) +
    geom_link_line(data = points)
  out <- export_ggchord_layout(p, include = c("feature", "link"))
  expect_setequal(out$feature$source_row, 1:2)
  expect_true(any(out$feature$component == "seq_region"))
  region <- out$feature[out$feature$component == "seq_region", ]
  expect_setequal(region$source_row, 2L)
  expect_setequal(region$accver, "A")
  expect_equal(unique(region$start), 100)
  expect_equal(unique(region$end), 300)
  expect_setequal(out$link$source_row, 1L)
})

test_that("independent ribbon layers retain their data and manual styles", {
  first <- smoke_ribbon()
  second <- transform(first, qstart = 500, qend = 600)
  p <- ggchord(smoke_sequences()) + geom_seq() +
    geom_link_ribbon(aes(ribbon_alpha = pident), data = first, fill = "orange") +
    geom_link_ribbon(aes(ribbon_alpha = pident),
      data = function(x) second, fill = "navy", alpha = .7) +
    scale_ribbon_alpha_continuous()
  built <- build_ggchord_smoke(p)
  expect_true(all(built$data[[2]]$ribbon_fill == "orange"))
  expect_true(all(is.finite(built$data[[2]]$ribbon_alpha)))
  expect_true(all(built$data[[3]]$ribbon_fill == "navy"))
  expect_true(all(built$data[[3]]$ribbon_alpha == .7))
  out <- get_chord_layout(p)$layer_inputs
  inputs <- lapply(Filter(function(x) !is.null(x$ribbon), out), function(x) x$ribbon)
  expect_setequal(vapply(inputs, function(x) x$qstart[1], numeric(1)), c(100, 500))
})

test_that("source metadata remains mappable without replacing computed geometry", {
  d <- data.frame(id = "A", lo = 100, hi = 260, direction = "+",
    label = "alpha", annotation_tag = "shown", x = NA_real_)
  p <- ggchord(smoke_sequences()[1, ]) + geom_seq() +
    geom_gene(aes(accver = id, start = lo, end = hi, strand = direction,
      anno = label, gene_fill = label), data = d) +
    geom_gene_label(aes(accver = id, start = lo, end = hi,
      strand = direction, anno = annotation_tag), data = d) +
    scale_gene_fill_manual(values = c(alpha = "#113355")) + coord_collinear()
  built <- build_ggchord_smoke(p)
  expect_true(nrow(built$data[[2]]) > 0L)
  expect_true(all(built$data[[2]]$gene_fill == "#113355"))
  expect_true(all(is.finite(built$data[[2]]$x)))
  expect_true(nrow(built$data[[3]]) > 0L)
  expect_true(all(built$data[[3]]$label == "shown"))
})
