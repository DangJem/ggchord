build_ggchord_smoke <- function(plot) {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  ggplot2::ggplot_build(plot)
}

test_that("constructor and coordinate follow the ggplot2 contract", {
  data(seq_data_example)
  p <- ggchord(seq_data_example, validate = "none") + geom_seq()
  expect_s3_class(p, "ggchord")
  expect_s3_class(coord_chord(), "CoordChord")
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  expect_error(coord_chord(fit = "manual"), "requires xlim and ylim")
  expect_error(get_chord_layout(), "plot must be supplied")
})

test_that("public geoms and stats return one standard layer", {
  data(seq_data_example)
  feature <- data.frame(
    seq_id = seq_data_example$seq_id[1], start = 1, end = 100,
    strand = "+", type = "CDS"
  )
  region <- feature[c("seq_id", "start", "end")]
  layers <- list(
    geom_seq(), geom_ribbon(), geom_gene(), geom_feature(data = feature),
    geom_axis(), geom_seq_label(), geom_gene_label(),
    geom_gene_label_repel(), geom_seq_region(data = region),
    geom_ribbon_highlight(), stat_ribbon_bundle(), stat_ribbon_density()
  )
  expect_true(all(vapply(layers, inherits, logical(1), "LayerInstance")))
  for (layer_fun in list(
    geom_seq, geom_ribbon, geom_gene, geom_feature, geom_axis,
    geom_seq_label, geom_gene_label, geom_gene_label_repel,
    geom_seq_region, geom_ribbon_highlight
  )) {
    expect_true(all(c(
      "mapping", "data", "position", "show.legend", "inherit.aes", "..."
    ) %in% names(formals(layer_fun))))
  }
})

test_that("main plotting layers build together as single layers", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(
    seq_data_example, ribbon_data_example, gene_data_example,
    validate = "none"
  ) +
    geom_seq() + geom_ribbon() + geom_gene() +
    geom_gene_label_repel() + geom_seq_label() + geom_axis()
  built <- build_ggchord_smoke(p)
  layout <- get_chord_layout(p)
  expect_length(built$data, 6L)
  expect_true(length(layout$seq_arcs) > 0L)
  expect_true(nrow(layout$ribbon_polys) > 0L)
  expect_true(nrow(layout$gene_polys) > 0L)
  expect_setequal(
    unique(built$data[[4]]$.component), c("segment", "text")
  )
  expect_gt(sum(built$data[[4]]$.component == "segment"), 0L)
  expect_setequal(
    unique(built$data[[6]]$.component), c("line", "tick", "text")
  )
})

test_that("removed parameters fail with direct migrations", {
  data(seq_data_example)
  expect_error(ggchord(seq_data_example, title = "x"), "labs")
  expect_error(geom_seq(seq_colors = "red"), "scale_seq_colour_manual")
  expect_error(geom_ribbon(ribbon_alpha = 0.5), "alpha")
  expect_error(geom_gene(gene_colors = "red"), "scale_gene_fill_manual")
  expect_error(geom_feature(type = "kind"), "feature_type")
  expect_error(geom_axis(show_axis = FALSE), "omitting geom_axis")
  expect_error(geom_seq_region(regions = data.frame()), "data")
  expect_error(geom_gene_label_repel(force = 1), "Removed")
})

test_that("complete themes own ggchord-specific elements", {
  themes <- list(
    theme_ggchord(), theme_ggchord_minimal(),
    theme_ggchord_publication(), theme_ggchord_dark()
  )
  expect_true(all(vapply(themes, inherits, logical(1), "theme")))
  custom <- theme_ggchord_publication(
    gene_label = element_text(colour = "purple", size = 8),
    axis_ticks = element_blank()
  )
  expect_equal(
    ggplot2::calc_element("ggchord.gene.label", custom)@colour, "purple"
  )
  expect_s3_class(
    ggplot2::calc_element("ggchord.axis.ticks", custom), "element_blank"
  )
  expect_false("theme_ggchord_elements" %in% getNamespaceExports("ggchord"))
})

test_that("curated ggplot2 helpers are exact reexports", {
  useful <- c(
    "aes", "after_stat", "after_scale", "annotate", "labs", "ggtitle",
    "guides", "theme", "element_blank", "element_line", "element_rect",
    "element_text", "margin", "rel", "ggsave", "last_plot", "waiver",
    "expansion"
  )
  exports <- getNamespaceExports("ggchord")
  expect_true(all(useful %in% exports))
  expect_false(any(c(
    "geom_point", "facet_wrap", "coord_cartesian", "theme_minimal"
  ) %in% exports))
  expect_identical(getExportedValue("ggchord", "theme"), ggplot2::theme)
  expect_identical(getExportedValue("ggchord", "ggtitle"), ggplot2::ggtitle)
})

test_that("data functions and role mappings are resolved once", {
  seq <- data.frame(seq_id = c("A", "B"), length = c(1000, 1000))
  feature <- data.frame(
    chromosome = "A", from = 100, to = 300, direction = "+",
    kind = "CDS", text = "gene A"
  )
  p <- ggchord(seq, validate = "none") +
    geom_seq(data = function(x) x) +
    geom_feature(
      aes(
        seq_id = chromosome, start = from, end = to, strand = direction,
        feature_type = kind, feature_label = text
      ),
      data = function(x) feature
    )
  expect_equal(unique(get_chord_layout(p)$gene_polys$anno), "CDS")
  expect_error(
    build_ggchord_smoke(
      ggchord(seq, validate = "none") +
        geom_seq(aes(seq_ring = after_stat(seq_id)))
    ),
    "cannot use after_stat"
  )
})

test_that("user scales win and missing scales are inferred", {
  seq <- data.frame(seq_id = c("A", "B"), length = c(1000, 1000))
  ribbon <- data.frame(
    qaccver = "A", saccver = "B", length = 101, pident = 90,
    qstart = 100, qend = 200, sstart = 300, send = 400,
    category = "shared"
  )
  inferred <- ggchord(seq, ribbon, validate = "none") +
    geom_seq() + geom_ribbon(aes(ribbon_fill = category))
  prepared <- ggchord:::prepare_ggchord_plot(inferred)
  expect_s3_class(prepared$scales$get_scales("ribbon_fill"), "ScaleDiscrete")

  manual <- inferred + scale_ribbon_fill_manual(values = c(shared = "black"))
  prepared_manual <- ggchord:::prepare_ggchord_plot(manual)
  expect_identical(unname(
    prepared_manual$scales$get_scales("ribbon_fill")$palette(1)
  ), "black")
})

test_that("automatic labels and advanced layout tools build", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  for (mode in c("aligned", "radial", "arc")) {
    p <- ggchord(
      seq_data_example, ribbon_data_example, gene_data_example,
      validate = "none"
    ) + geom_seq() + geom_gene() +
      geom_gene_label_repel(gene_label_layout = mode)
    expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  }
  stat_plot <- ggchord(
    seq_data_example, ribbon_data_example, validate = "none"
  ) + geom_seq() +
    stat_ribbon_bundle(aes(ribbon_alpha = after_stat(density)), bins = 10)
  stat_build <- build_ggchord_smoke(stat_plot)
  expect_true(all(is.finite(stat_build$data[[2]]$ribbon_alpha)))

  rings <- transform(seq_data_example[1:2, ], ring = c("outer", "inner"))
  ring_plot <- ggchord(rings, validate = "none") +
    geom_seq(aes(seq_ring = ring)) +
    scale_seq_ring_manual(values = c(outer = 2, inner = 1))
  expect_equal(unname(get_chord_layout(ring_plot)$seq_ring_radius), c(2, 1))
})

test_that("feature shapes, regions and highlights remain composable", {
  data(seq_data_example)
  data(ribbon_data_example)
  feature <- data.frame(
    seq_id = seq_data_example$seq_id[1],
    start = c(100, 700, 1300, 1900), end = c(500, 1100, 1700, 2300),
    strand = c("+", "-", "+", "-"),
    type = c("CDS", "tRNA", "repeat", "promoter")
  )
  region <- data.frame(
    seq_id = seq_data_example$seq_id[1], start = 2500, end = 3000
  )
  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() + geom_ribbon() +
    geom_feature(aes(feature_shape = type), data = feature) +
    scale_feature_shape_manual(values = c(
      CDS = "arrow", tRNA = "block", "repeat" = "chevron",
      promoter = "lollipop"
    )) +
    geom_seq_region(data = region) +
    geom_ribbon_highlight(ribbon_ids = 1)
  layout <- get_chord_layout(p)
  expect_setequal(
    unique(layout$gene_polys$feature_shape),
    c("arrow", "block", "chevron", "lollipop")
  )
  expect_true(nrow(layout$region_polys) > 0L)
  expect_true(nrow(layout$ribbon_highlight_polys) > 0L)
})

test_that("layout export and static viewer use explicit plots", {
  seq <- data.frame(seq_id = c("A", "B"), length = c(1000, 1000))
  ribbon <- data.frame(
    qaccver = "A", saccver = "B", length = 101, pident = 90,
    qstart = 100, qend = 200, sstart = 300, send = 400
  )
  p <- ggchord(seq, ribbon, validate = "none") + geom_seq() + geom_ribbon()
  exported <- export_ggchord_layout(p, include = c("seq", "ribbon"))
  expect_s3_class(exported, "ggchord_layout_export")
  expect_true(all(c("layer_id", "source_row") %in% names(exported$ribbon)))

  preview <- view_ggchord(
    p, width = 2, height = 1, units = "in", dpi = 50,
    device = "png", viewer = "none"
  )
  expect_true(file.exists(preview))
  expect_true(file.exists(sub("\\.png$", ".html", preview)))
})

test_that("plot-owned layouts are deterministic and isolated", {
  p1 <- ggchord(
    data.frame(seq_id = "A", length = 100), validate = "none"
  ) + geom_seq()
  p2 <- ggchord(
    data.frame(seq_id = "B", length = 200), validate = "none"
  ) + geom_seq()
  first <- get_chord_layout(p1)
  repeated <- get_chord_layout(p1)
  expect_identical(first$seq_arcs, repeated$seq_arcs)
  expect_equal(first$seqs, "A")
  expect_equal(get_chord_layout(p2)$seqs, "B")
})
