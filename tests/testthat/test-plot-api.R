build_ggchord_smoke <- function(plot) {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  ggplot2::ggplot_build(plot)
}

test_that("constructor, coordinate and public layers build", {
  data(seq_data_example); data(ribbon_data_example); data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example,
               validate = "none") +
    geom_seq() + geom_link_ribbon() + geom_gene() +
    geom_gene_label_repel() + geom_seq_label()
  expect_s3_class(p, "ggchord")
  expect_s3_class(coord_chord(), "CoordChord")
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  expect_gt(nrow(get_chord_layout(p)$ribbon_polys), 0L)
  expect_error(coord_chord(fit = "manual"), "requires xlim and ylim")
  expect_error(get_chord_layout(), "plot must be supplied")

  feature <- data.frame(accver = seq_data_example$accver[1], start = 1,
    end = 100, strand = "+", type = "CDS")
  region <- feature[c("accver", "start", "end")]
  layers <- list(geom_seq(), geom_link_ribbon(), geom_link_line(), geom_gene(),
    geom_feature(data = feature), geom_seq_label(), geom_gene_label(),
    geom_gene_label_repel(), geom_seq_region(data = region),
    stat_ribbon_bundle(), stat_ribbon_density())
  expect_true(all(vapply(layers, inherits, logical(1), "LayerInstance")))
})

test_that("standard data, mapping and independent layers are preserved", {
  data(seq_data_example); data(ribbon_data_example)
  subset_a <- ribbon_data_example[1:2, ]
  subset_b <- ribbon_data_example[3:4, ]
  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() +
    geom_link_ribbon(data = subset_a, aes(ribbon_alpha = pident)) +
    geom_link_ribbon(data = function(x) subset_b, fill = "orange")
  built <- build_ggchord_smoke(p)
  layout <- get_chord_layout(p)
  ribbon_layers <- Filter(function(x) "ribbon" %in% names(x),
                          layout$layer_geometry)
  expect_equal(unname(vapply(ribbon_layers, function(x) nrow(x$ribbon), integer(1))),
               c(400L, 400L))
  expect_equal(unname(vapply(layout$layer_inputs[names(ribbon_layers)],
    function(x) nrow(x$ribbon), integer(1))), c(2L, 2L))
  expect_false(identical(built$data[[2]]$ribbon_fill,
                         built$data[[3]]$ribbon_fill))
})

test_that("role scales, themes and colour aliases remain composable", {
  data(seq_data_example); data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq(color = "black") + geom_link_ribbon(color = "navy") +
    scale_ribbon_fill_stepsn(colours = c("white", "steelblue")) +
    theme_ggchord_minimal(legend.ribbon.position = "bottom")
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  expect_error(geom_seq(colour = "red", color = "blue"), "only one")
  expect_error(geom_link_ribbon(colour = "red", color = "blue"), "only one")
  expect_identical(aes, ggplot2::aes)
  expect_identical(ggsave, ggplot2::ggsave)
})

test_that("label modes and overlap policies build", {
  data(seq_data_example); data(gene_data_example)
  base <- ggchord(seq_data_example, gene_data = gene_data_example,
                  validate = "none") + geom_seq() + geom_gene()
  for (mode in c("radial", "auto", "arc")) {
    expect_s3_class(build_ggchord_smoke(
      base + geom_gene_label_repel(gene_label_layout = mode)
    ), "ggplot_built")
  }
  for (overlap in c("fade", "clip", "show")) {
    expect_s3_class(build_ggchord_smoke(
      base + geom_gene_label_repel(gene_label_segment_overlap = overlap)
    ), "ggplot_built")
  }
})

test_that("feature shapes, regions and layout export are public behavior", {
  data(seq_data_example)
  feature <- data.frame(accver = seq_data_example$accver[1],
    start = c(100, 300, 500, 700), end = c(200, 400, 600, 800),
    strand = "+", type = c("CDS", "tRNA", "repeat", "promoter"))
  region <- feature[1, c("accver", "start", "end")]
  p <- ggchord(seq_data_example) + geom_seq() +
    geom_feature(aes(feature_shape = type), data = feature) +
    scale_feature_shape_manual(values = c(CDS = "arrow", tRNA = "block",
      "repeat" = "chevron", promoter = "lollipop")) +
    geom_seq_region(data = region)
  layout <- get_chord_layout(p)
  expect_setequal(unique(layout$gene_polys$feature_shape),
                  c("arrow", "block", "chevron", "lollipop"))
  exported <- export_ggchord_layout(p, include = c("seq", "feature"))
  expect_s3_class(exported, "ggchord_layout_export")
  expect_gt(nrow(exported$feature), 0L)
})

test_that("canonical ribbon statistics do not use a compatibility entry", {
  expect_false("geom_ribbon" %in% getNamespaceExports("ggchord"))
  expect_no_warning(stat_ribbon_bundle())
})
