test_that("constructor and core sequence, ribbon and gene layers build", {
  genes <- smoke_genes()[c("accver", "start", "end", "strand")]
  colours <- c("+" = "#56887C", "-" = "#B58169")
  p <- ggchord(smoke_sequences(), smoke_ribbon(), genes, validate = "none") +
    geom_seq() + geom_link_ribbon() + geom_gene() + geom_seq_label() +
    scale_gene_fill_manual(values = colours)
  expect_s3_class(p, "ggchord")
  built <- build_ggchord_smoke(p)
  expect_s3_class(built, "ggplot_built")
  band <- ggchord(smoke_sequences()[1, ]) + geom_seq(seq_style = "band") + coord_circular()
  expect_s3_class(ggplot2::ggplotGrob(band), "gtable")
  # F53: optional annotation text must not disable the strand fill mapping.
  expect_setequal(unique(built$data[[3]]$gene_fill), unname(colours))
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
  first <- transform(smoke_ribbon(), pident = 95)
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

test_that("Dotplot retains reverse HSP endpoints and reversible axis positions", {
  hits <- smoke_ribbon()
  hits$sstart <- 600
  hits$send <- 200
  p <- ggchord(smoke_sequences(), hits, validate = "none") +
    geom_alignment_segment() + coord_dotplot()
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  out <- export_ggchord_layout(p)$alignment
  expect_true(all(out$y > out$yend))
  back <- map_ggchord_position(p, out$saccver, out$y,
    inverse = TRUE, axis = "y")
  expect_equal(back$position, hits$sstart)
})

test_that("Radial spokes retain complete intervals and reversible native positions", {
  p <- ggchord(smoke_sequences(), smoke_ribbon(), validate = "none") +
    geom_seq() + geom_link_ribbon() + coord_radial_synteny()
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  expect_equal(nrow(export_ggchord_layout(p)$radial_endpoints), 1L)
  mapped <- map_ggchord_position(p, "A", 100)
  back <- map_ggchord_position(p, "A", mapped[c("x", "y")], inverse = TRUE)
  expect_equal(back$position, 100)
})

test_that("Fusion retains source intervals across independently measured panels", {
  center <- ggchord(smoke_sequences()[1, ]) + geom_seq() + coord_circular() + ggplot2::scale_x_reverse()
  detail <- ggchord(smoke_sequences()) + geom_seq() + coord_collinear()
  source <- smoke_ribbon()
  fusion <- ggchord_fusion(center) + fusion_panel(detail, "detail") +
    fusion_ribbon(source, qpanel = "center", spanel = "detail", fill = "steelblue") +
    ggplot2::labs(title = "Source intervals")
  out <- export_ggchord_fusion(fusion)
  expect_true(all(out$connections$retained))
  expect_identical(out$metadata$plot_decoration$labels$title, "Source intervals")
  expect_identical(out$relationship_data[[1]], source)
  query <- out$fronts[out$fronts$endpoint == "q", ]
  expect_lt(stats::cor(query$local_x, query$x), 0)
  locus <- data.frame(accver = c("A", "B"), start = 100, end = 200,
    instance_note = "ordinary metadata", detail_instance_note = "also metadata")
  focused <- focus_ggchord_data(smoke_sequences(), loci = locus)
  local <- ggchord(focused$seq_data) + geom_seq() + coord_collinear()
  nested <- ggchord_fusion(local) + fusion_panel(local, "leaf")
  locus$detail_panel <- c("first/leaf", "second")
  magnified <- export_ggchord_fusion(magnify_ggchord(detail, locus,
    list(first = nested, second = local), side = c(first = "top", second = "bottom")))
  expect_identical(magnified$focus_regions$detail_panel, c("first/leaf", "second"))
  expect_equal(magnified$focus_regions$detail_origin, c(100, 100))
  expect_identical(magnified$focus_regions$region_row, 1:2)
  expect_identical(magnified$focus_region_data, locus)
  expect_equal(nrow(magnified$leaders), 4L)
  expect_equal(nrow(magnified$connections), 0L)
})
