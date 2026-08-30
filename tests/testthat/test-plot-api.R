build_ggchord_smoke <- function(plot) {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  ggplot2::ggplot_build(plot)
}

test_that("ggchord creates a plot and rejects invalid sequence data", {
  data(seq_data_example)

  p <- ggchord(seq_data_example, validate = "none")
  expect_s3_class(p, "ggchord")
  expect_s3_class(coord_chord(), "Coord")

  expect_error(
    ggchord(data.frame(id = "A", length = 100)),
    "seq_data"
  )
})

test_that("coord_chord controls rotation, fit and explicit limits", {
  seq <- data.frame(seq_id = c("A", "B"), length = c(1000, 1200))
  p <- ggchord(seq, validate = "none") +
    geom_seq() +
    coord_chord(
      rotation = 90, xlim = c(-5, 6), ylim = c(-7, 8),
      fit = "labels", expand = FALSE
    )

  layout <- get_chord_layout(p)
  expect_equal(layout$rotation, 90)
  prepared <- ggchord:::prepare_ggchord_plot(p)
  expect_equal(prepared$coordinates$limits$x, c(-5, 6))
  expect_equal(prepared$coordinates$limits$y, c(-7, 8))
  expect_false(prepared$coordinates$expand)

  expect_error(coord_chord(fit = "manual"), "requires xlim and ylim")

  p_cartesian <- suppressMessages(
    ggchord(seq, validate = "none") + geom_seq() +
      ggplot2::coord_cartesian(xlim = c(-2, 2))
  )
  prepared_cartesian <- ggchord:::prepare_ggchord_plot(p_cartesian)
  expect_s3_class(prepared_cartesian$coordinates, "CoordCartesian")
  expect_false(isTRUE(prepared_cartesian$coordinates$ggchord_coord))
  expect_equal(prepared_cartesian$coordinates$limits$x, c(-2, 2))
})

test_that("ggchord themes and guides use registered role elements", {
  for (fun in list(
    theme_ggchord, theme_ggchord_minimal,
    theme_ggchord_dark, theme_ggchord_publication
  )) {
    expect_s3_class(fun(), "theme")
  }
  expect_s3_class(
    ggplot2::calc_element("ggchord.axis.line", theme_ggchord()),
    "element_line"
  )
  expect_s3_class(guide_ggchord_legend(), "GuideLegend")
  expect_s3_class(guide_ggchord_colourbar(), "GuideColourbar")
})

test_that("the main plotting layers build together", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p <- ggchord(
    seq_data_example, ribbon_data_example, gene_data_example,
    validate = "none"
  ) +
    geom_seq() +
    geom_ribbon() +
    geom_gene() +
    geom_axis() +
    geom_gene_label() +
    geom_seq_label() +
    geom_axis()

  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout(p)
  expect_true(length(layout$seq_arcs) > 0)
  expect_true(nrow(layout$ribbon_polys) > 0)
  expect_true(nrow(layout$gene_polys) > 0)
})

test_that("all automatic gene-label layouts build", {
  data(seq_data_example)
  data(gene_data_example)

  for (mode in c("aligned", "radial", "arc")) {
    p <- ggchord(seq_data_example, gene_data = gene_data_example,
                 validate = "none") +
      geom_seq() +
      geom_gene() +
      geom_gene_label_repel(gene_label_layout = mode)
    expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  }

  expect_error(
    geom_gene_label_repel(force = 1),
    "Removed"
  )
})

test_that("region and ribbon-highlight layers build", {
  data(seq_data_example)
  data(ribbon_data_example)

  regions <- data.frame(
    seq_id = seq_data_example$seq_id[1],
    start = 100,
    end = 500
  )
  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() +
    geom_ribbon() +
    geom_seq_region(regions = regions) +
    geom_ribbon_highlight(ribbon_ids = 1)

  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout(p)
  expect_true(nrow(layout$region_polys) > 0)
  expect_true(nrow(layout$ribbon_highlight_polys) > 0)
})

test_that("feature and sequence-group layers build", {
  data(seq_data_example)

  feature <- data.frame(
    seq_id = seq_data_example$seq_id[1],
    start = 100,
    end = 500,
    strand = "+",
    type = "CDS"
  )
  groups <- stats::setNames(
    rep(c("group-a", "group-b"), length.out = nrow(seq_data_example)),
    seq_data_example$seq_id
  )
  p <- ggchord(seq_data_example, validate = "none") +
    geom_seq(seq_group = groups) +
    geom_feature(feature)

  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout(p)
  expect_true(nrow(layout$gene_polys) > 0)
  expect_true(nrow(layout$group_labels) > 0)
})

test_that("feature category and region outline survive geometry generation", {
  data(seq_data_example)
  feature <- data.frame(
    seq_id = seq_data_example$seq_id[1], start = 100, end = 500,
    strand = "+", type = "CDS", category = "coding", label = "display"
  )
  region <- data.frame(
    seq_id = seq_data_example$seq_id[1], start = 600, end = 900
  )
  p <- ggchord(seq_data_example, validate = "none") +
    geom_seq(seq_curvature = 0.4) +
    geom_feature(feature, category = "category") +
    geom_seq_region(regions = region, region_color = "#123456",
                    region_side = "auto")
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout(p)
  expect_true(all(layout$gene_polys$anno == "coding"))
  expect_true(all(layout$region_polys$colour == "#123456"))

  expect_error(geom_ribbon_highlight(ribbon_ids = 0), "positive")
})

test_that("same-type layers keep independent data and mapped columns", {
  seq <- data.frame(seq_id = c("A", "B"), length = c(1000, 1000))
  first <- data.frame(
    chromosome = "A", from = 100, to = 200, direction = "+",
    category = "first"
  )
  second <- data.frame(
    seq_id = "B", start = 300, end = 400, strand = "-", anno = "second"
  )
  r1 <- data.frame(seq_id = "A", start = 450, end = 500, category = "r1")
  r2 <- data.frame(seq_id = "B", start = 550, end = 600, category = "r2")

  p <- ggchord(seq, validate = "none") +
    geom_seq() +
    geom_gene(
      data = first,
      mapping = aes(seq_id = chromosome, start = from, end = to,
                    strand = direction, anno = category,
                    gene_fill = category)
    ) +
    geom_gene(data = second, mapping = aes(gene_fill = anno)) +
    geom_seq_region(data = r1, mapping = aes(region_fill = category)) +
    geom_seq_region(data = r2, mapping = aes(region_fill = category))

  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout(p, build = FALSE)
  gene_ids <- vapply(
    p$layers[vapply(p$layers, function(x) identical(x$ggchord_type, "gene_poly"),
                    logical(1))],
    function(x) x$ggchord_layer_id, character(1)
  )
  region_ids <- vapply(
    p$layers[vapply(p$layers, function(x) identical(x$ggchord_type, "seq_region"),
                    logical(1))],
    function(x) x$ggchord_layer_id, character(1)
  )
  expect_equal(length(unique(gene_ids)), 2)
  expect_equal(length(unique(region_ids)), 2)
  expect_equal(unique(layout$layer_geometry[[gene_ids[1]]]$gene_poly$anno), "first")
  expect_equal(unique(layout$layer_geometry[[gene_ids[2]]]$gene_poly$anno), "second")
  expect_equal(unique(layout$layer_geometry[[region_ids[1]]]$seq_region$category), "r1")
  expect_equal(unique(layout$layer_geometry[[region_ids[2]]]$seq_region$category), "r2")
})

test_that("plot-specific layout retrieval does not use another plot's cache", {
  p1 <- ggchord(data.frame(seq_id = "A", length = 100), validate = "none") +
    geom_seq()
  p2 <- ggchord(data.frame(seq_id = "B", length = 200), validate = "none") +
    geom_seq()
  expect_equal(get_chord_layout(p1)$seqs, "A")
  expect_equal(get_chord_layout(p2)$seqs, "B")
  expect_equal(get_chord_layout(p1, build = FALSE)$seqs, "A")
})

test_that("sequence, ribbon, axis and label data mappings are honoured", {
  seq_base <- data.frame(seq_id = c("A", "B"), length = c(1000, 1000))
  seq_mapped <- data.frame(chromosome = c("A", "B"), bases = c(1000, 1000))
  ribbons <- data.frame(
    query = "A", subject = "B", span = 100, identity = 95,
    q_from = 1, q_to = 100, s_from = 200, s_to = 101,
    score = 0.6
  )
  seq_subset <- data.frame(seq_id = "A", length = 1000)
  p <- ggchord(seq_base, validate = "none") +
    geom_seq(
      data = seq_mapped,
      mapping = aes(seq_id = chromosome, length = bases)
    ) +
    geom_ribbon(
      data = ribbons,
      mapping = aes(qaccver = query, saccver = subject, length = span,
                    pident = identity, qstart = q_from, qend = q_to,
                    sstart = s_from, send = s_to, ribbon_alpha = score)
    ) +
    geom_seq_label(data = seq_subset) +
    geom_axis(data = seq_subset)

  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  layout <- get_chord_layout(p, build = FALSE)
  ribbon_id <- p$layers[[which(vapply(
    p$layers, function(x) identical(x$ggchord_type, "ribbon"), logical(1)
  ))]]$ggchord_layer_id
  label_id <- p$layers[[which(vapply(
    p$layers, function(x) identical(x$ggchord_type, "seq_label"), logical(1)
  ))]]$ggchord_layer_id
  expect_equal(unique(layout$layer_geometry[[ribbon_id]]$ribbon$score), 0.6)
  expect_equal(unique(layout$layer_geometry[[label_id]]$seq_label$seq_id), "A")
})

test_that("role-specific scales coexist without replacing one another", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  region <- data.frame(
    seq_id = seq_data_example$seq_id[1], start = 100, end = 500,
    category = "focus"
  )
  seq_values <- stats::setNames(
    rep_len(c("#0072B2", "#D55E00"), nrow(seq_data_example)),
    seq_data_example$seq_id
  )
  p <- ggchord(
    seq_data_example, ribbon_data_example, gene_data_example,
    validate = "none"
  ) +
    geom_seq() +
    geom_ribbon() +
    geom_gene() +
    geom_seq_region(
      data = region, mapping = aes(region_fill = category), show_legend = TRUE
    ) +
    scale_seq_colour_manual(values = seq_values) +
    scale_ribbon_fill_gradientn(colours = c("#F7FBFF", "#08306B")) +
    scale_gene_fill_manual(values = c("+" = "#D55E00", "-" = "#0072B2")) +
    scale_region_fill_manual(values = c(focus = "#E69F00")) +
    scale_seq_position_continuous(
      breaks = c(0, 100), labels = c("start", "100 bp")
    )

  expect_true(all(vapply(
    c("seq_colour", "ribbon_fill", "gene_fill", "region_fill", "seq_position"),
    p$scales$has_scale, logical(1)
  )))
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  axis_labels <- unique(stats::na.omit(get_chord_layout(p, FALSE)$axis_ticks$label))
  expect_true(all(c("start", "100 bp") %in% axis_labels))
  expect_s3_class(scale_seq_position_continuous(), "ScaleContinuous")
})

test_that("legacy scale arguments warn and conflict with role scales", {
  seq <- data.frame(seq_id = c("A", "B"), length = c(100, 100))
  old <- ggchord(seq, validate = "none") +
    geom_seq(seq_colors = c(A = "red", B = "blue"))
  expect_warning(build_ggchord_smoke(old), "seq_colors.*deprecated")

  conflict <- old + scale_seq_colour_manual(
    values = c(A = "black", B = "grey50")
  )
  expect_error(build_ggchord_smoke(conflict), "conflicts.*seq_colour")
})
