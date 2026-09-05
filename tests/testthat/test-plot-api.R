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
    geom_seq(), geom_link_ribbon(), geom_gene(), geom_feature(data = feature),
    geom_seq_label(), geom_gene_label(),
    geom_gene_label_repel(), geom_seq_region(data = region),
    geom_ribbon_highlight(), stat_ribbon_bundle(), stat_ribbon_density()
  )
  expect_true(all(vapply(layers, inherits, logical(1), "LayerInstance")))
  for (layer_fun in list(
    geom_seq, geom_link_ribbon, geom_gene, geom_feature,
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
    geom_seq() + geom_link_ribbon() + geom_gene() +
    geom_gene_label_repel() + geom_seq_label()
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
    unique(built$data[[6]]$.component),
    c("line", "major_tick", "minor_tick", "text")
  )
})

test_that("removed parameters fail with direct migrations", {
  data(seq_data_example)
  expect_error(ggchord(seq_data_example, title = "x"), "labs")
  expect_error(geom_seq(seq_colors = "red"), "scale_seq_colour_manual")
  expect_error(geom_link_ribbon(ribbon_alpha = 0.5), "alpha")
  expect_error(geom_gene(gene_colors = "red"), "scale_gene_fill_manual")
  expect_error(geom_feature(type = "kind"), "feature_type")
  expect_false("geom_axis" %in% getNamespaceExports("ggchord"))
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
    gene.label = element_text(colour = "purple", size = 8),
    axis.ticks = element_blank()
  )
  expect_equal(
    ggplot2::calc_element("ggchord.gene.label", custom)@colour, "purple"
  )
  expect_s3_class(
    ggplot2::calc_element("ggchord.axis.ticks", custom), "element_blank"
  )
  expect_false("theme_ggchord_elements" %in% getNamespaceExports("ggchord"))
  expect_false(any(c("base_size", "base_family") %in%
                     names(formals(theme_ggchord))))
  expect_error(theme_ggchord(axis_ticks = element_blank()), "axis.ticks")
  inherited <- theme_ggchord(text = element_text(size = 10))
  expect_equal(ggplot2::calc_element("ggchord.axis.text", inherited)@size, 7.2)
})

test_that("automatic sequence axes are theme controlled", {
  data(seq_data_example)
  p <- ggchord(seq_data_example, validate = "none") + geom_seq()
  layout <- get_chord_layout(p)
  expect_gt(nrow(layout$axis_lines), 0L)
  expect_true(all(c(TRUE, FALSE) %in% layout$axis_ticks$is_major))

  hidden <- p + theme_ggchord(axis = element_blank())
  expect_equal(nrow(get_chord_layout(hidden)$axis_ticks), 0L)
  expect_error(
    theme_ggchord(axis.gap = grid::unit(c(1, 2), "mm")),
    "one grid::unit"
  )

  scaled <- p + scale_seq_position_continuous(
    breaks = c(0, 100), minor_breaks = 50,
    labels = c("start", "100")
  )
  ticks <- get_chord_layout(scaled)$axis_ticks
  expect_true(all(na.omit(unique(ticks$label)) %in% c("start", "100")))
})

test_that("role legends inherit independently and explicit guides win", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() + geom_link_ribbon() +
    theme_ggchord(
      legend.seq.position = "right",
      legend.ribbon.position = "left",
      legend.gene = element_blank(),
      legend.feature.position = "top",
      legend.region.position = "bottom",
      legend.ribbon.key.height = grid::unit(42, "mm"),
      legend.ribbon.ticks.length = grid::unit(1, "mm")
    )
  expect_equal(ggchord:::ggchord_role_guide_spec(p, "seq")$position, "right")
  ribbon <- ggchord:::ggchord_role_guide_spec(p, "ribbon", TRUE)
  expect_equal(ribbon$position, "left")
  expect_equal(as.numeric(ribbon$theme$legend.key.height), 42)
  expect_true(ggchord:::ggchord_role_guide_spec(p, "gene")$hidden)
  expect_equal(
    ggchord:::ggchord_role_guide_spec(p, "feature")$position, "top"
  )
  expect_equal(
    ggchord:::ggchord_role_guide_spec(p, "region")$position, "bottom"
  )
  ordinary <- p + theme(legend.key.width = grid::unit(9, "mm"))
  expect_equal(
    as.numeric(ggchord:::ggchord_role_guide_spec(ordinary, "seq")$
      theme$legend.key.width), 9
  )

  explicit <- p + scale_ribbon_fill_gradientn(
    colours = c("navy", "gold"), guide = "none"
  )
  prepared <- ggchord:::prepare_ggchord_plot(explicit)
  expect_identical(prepared$scales$get_scales("ribbon_fill")$guide, "none")
})

test_that("default role guides occupy opposite plot edges", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, validate = "none") +
    geom_seq() + geom_link_ribbon()
  expect_equal(
    ggchord:::ggchord_role_guide_spec(p, "ribbon", TRUE)$position,
    "left"
  )
  expect_equal(ggchord:::ggchord_role_guide_spec(p, "seq")$position, "right")
  common <- p + theme_ggchord(legend.position = "bottom")
  expect_equal(
    ggchord:::ggchord_role_guide_spec(common, "ribbon", TRUE)$position,
    "bottom"
  )
  expect_equal(ggchord:::ggchord_role_guide_spec(common, "seq")$position, "bottom")
})

test_that("colour and color aliases are symmetric and unambiguous", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, validate = "none") +
    geom_seq(aes(seq_color = seq_id))
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
  expect_identical(scale_seq_color_manual, scale_seq_colour_manual)
  expect_identical(scale_ribbon_color_manual, scale_ribbon_colour_manual)
  expect_identical(guide_ggchord_colorbar, guide_ggchord_colourbar)
  ribbon_layer <- geom_link_ribbon(aes(ribbon_color = pident))
  expect_true("ribbon_colour" %in% names(ribbon_layer$ggchord_input_mapping))
  expect_s3_class(
    scale_ribbon_fill_gradientn(colors = c("navy", "gold")),
    "ScaleContinuous"
  )
  expect_error(suppressWarnings(
    geom_seq(aes(seq_color = seq_id, seq_colour = seq_id))
  ), "normalisation")
  expect_error(geom_link_ribbon(colour = "red", color = "blue"), "only one")
  expect_error(
    scale_ribbon_fill_gradientn(
      colours = c("navy", "gold"), colors = c("black", "white")
    ),
    "only one"
  )
  region_color <- data.frame(
    seq_id = seq_data_example$seq_id[1], start = 1, end = 100,
    color = "red"
  )
  region_colour <- region_color
  names(region_colour)[names(region_colour) == "color"] <- "colour"
  for (region in list(region_color, region_colour)) {
    rp <- ggchord(seq_data_example, validate = "none") +
      geom_seq() + geom_seq_region(data = region)
    expect_equal(unique(get_chord_layout(rp)$region_polys$zregionfill), "red")
  }
  both <- cbind(region_color, colour = "blue")
  expect_error(get_chord_layout(
    ggchord(seq_data_example, validate = "none") +
      geom_seq_region(data = both)
  ), "only one")
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
    geom_seq() + geom_link_ribbon(aes(ribbon_fill = category))
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
  for (mode in c("auto", "radial", "arc")) {
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

test_that("automatic labels use compact sequence rails and adaptive fitting", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(
    seq_data_example, ribbon_data_example, gene_data_example,
    validate = "none"
  ) +
    geom_seq(
      seq_radius = c(3.3, 2.5, 1.8, 1.25),
      seq_orientation = c(1, -1, 1, -1)
    ) +
    geom_gene() +
    geom_gene_label_repel(gene_label_layout = "auto", gene_label_fit = "none")
  layout <- get_chord_layout(p)
  labels <- layout$gene_labels
  frame <- ggchord:::ggchord_label_curve_frame(labels, layout$seq_arcs)
  direction <- ifelse(abs(frame$outward_x) >= abs(frame$outward_y),
    ifelse(frame$outward_x < 0, "left", "right"),
    ifelse(frame$outward_y < 0, "bottom", "top"))
  expect_setequal(unique(direction), c("left", "right", "top", "bottom"))
  vertical <- split(labels[direction %in% c("left", "right"), ],
                    direction[direction %in% c("left", "right")])
  expect_true(all(vapply(vertical, function(x) {
    length(unique(round(x$text_x, 8))) == 1L
  }, logical(1))))
  expect_false(ggchord:::ggchord_label_box_conflicts(
    labels,
    units_per_inch = layout$text_units_per_inch,
    box_padding = 0
  ))
  count_crossings <- function(leaders) {
    crossings <- 0L
    if (nrow(leaders) < 2L) return(crossings)
    for (i in seq_len(nrow(leaders) - 1L)) {
      for (j in (i + 1L):nrow(leaders)) {
        if (leaders$group[i] == leaders$group[j]) next
        crossings <- crossings + ggchord:::ggchord_segments_cross(
          leaders$x0[i], leaders$y0[i], leaders$x1[i], leaders$y1[i],
          leaders$x0[j], leaders$y0[j], leaders$x1[j], leaders$y1[j]
        )
      }
    }
    crossings
  }
  expect_equal(count_crossings(layout$gene_label_segments), 0L)

  # Rotating curved sequences can assign labels from different sequences to
  # one side column. Their shared endpoint order must remain crossing-free.
  rotated <- get_chord_layout(
    ggchord(
      seq_data_example, ribbon_data_example, gene_data_example,
      validate = "none"
    ) +
      geom_seq(
        seq_radius = c(3.3, 2.5, 1.8, 1.25),
        seq_orientation = c(1, -1, 1, -1),
        seq_curvature = c(0.8, 1.2, 0.7, 1.1),
        seq_gap = c(0.03, 0.06, 0.04, 0.08)
      ) +
      geom_gene() +
      geom_gene_label_repel(gene_label_layout = "auto", gene_label_fit = "none") +
      coord_chord(rotation = 35)
  )
  expect_equal(count_crossings(rotated$gene_label_segments), 0L)

  fitted <- function(method) {
    get_chord_layout(
      ggchord(
        seq_data_example, ribbon_data_example, gene_data_example,
        validate = "none"
      ) +
        geom_seq(
          seq_radius = c(3.3, 2.5, 1.8, 1.25),
          seq_orientation = -1
        ) +
        geom_gene_label_repel(gene_label_fit = method)
    )$gene_labels$text
  }
  expect_true(any(grepl("\n", fitted("wrap"), fixed = TRUE)))
  expect_true(any(grepl("…$", fitted("ellipsis"))))
  expect_error(
    geom_gene_label_repel(gene_label_max_lines = 1.5),
    "positive integer"
  )
})

test_that("covered gene leaders can fade, clip or remain visible", {
  labels <- data.frame(
    text = c("target", "blocking label"),
    text_x = c(2, 1), text_y = c(0, 0), text_angle = 0,
    size = 2.5, hjust = 0.5, vjust = 0.5
  )
  segment <- data.frame(
    x0 = 0, y0 = 0, x1 = 2, y1 = 0, group = 1L
  )
  split_leader <- function(mode, alpha = 0.18) {
    ggchord:::ggchord_clip_segments_to_labels(
      segment, labels, units_per_inch = 1,
      overlap = mode, overlap_alpha = alpha
    )
  }
  faded <- split_leader("fade", 0.3)
  clipped <- split_leader("clip")
  shown <- split_leader("show")
  expect_true(any(faded$occluded))
  expect_equal(unique(faded$alpha[faded$occluded]), 0.3)
  expect_true(all(clipped$alpha == 1) && !any(clipped$occluded))
  expect_equal(nrow(shown), 1L)
  expect_equal(shown[c("x0", "y0", "x1", "y1")], segment[1:4])

  layer <- geom_gene_label_repel(
    gene_label_segment_overlap = "clip",
    gene_label_segment_overlap_alpha = 0.4
  )
  expect_equal(layer$ggchord_params$gene_label_segment_overlap, "clip")
  expect_equal(layer$ggchord_params$gene_label_segment_overlap_alpha, 0.4)
  expect_error(
    geom_gene_label_repel(gene_label_segment_overlap = "route"),
    "should be one of"
  )
  expect_error(
    geom_gene_label_repel(gene_label_segment_overlap_alpha = 2),
    "finite number in \\[0, 1\\]"
  )
})

test_that("manual auto-side labels extend toward their actual side", {
  seq <- data.frame(seq_id = "A", length = 1000)
  genes <- data.frame(
    seq_id = "A", start = c(150, 650), end = c(300, 800),
    strand = c("+", "-"), anno = c("outside", "inside")
  )
  p <- ggchord(seq, gene_data = genes, validate = "none") +
    geom_seq() +
    geom_gene_label(
      gene_label_side = "auto",
      gene_label_radial_offset = 0.25
    )
  layout <- get_chord_layout(p)
  labels <- layout$gene_labels
  frame <- ggchord:::ggchord_label_curve_frame(labels, layout$seq_arcs)
  expect_true(any(frame$signed_distance < 0))
  dx <- labels$text_x - frame$curve_x
  dy <- labels$text_y - frame$curve_y
  horizontal_side <- abs(dx) >= 0.75 * abs(dy)
  if (any(horizontal_side)) {
    expect_equal(
      labels$hjust[horizontal_side],
      ifelse(dx[horizontal_side] >= 0, 0, 1)
    )
  }
  if (any(!horizontal_side)) {
    expect_equal(
      labels$vjust[!horizontal_side],
      ifelse(dy[!horizontal_side] >= 0, 0, 1)
    )
  }
})

test_that("single-strand gene guides train with matching arrow keys", {
  seq <- data.frame(seq_id = "A", length = 1000)
  genes <- data.frame(
    seq_id = "A", start = c(100, 500), end = c(250, 700),
    strand = "+", anno = c("one", "two")
  )
  p <- ggchord(seq, gene_data = genes, validate = "none") +
    geom_seq() + geom_gene()
  expect_s3_class(build_ggchord_smoke(p), "ggplot_built")
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
    geom_seq() + geom_link_ribbon() +
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
  p <- ggchord(seq, ribbon, validate = "none") + geom_seq() + geom_link_ribbon()
  exported <- export_ggchord_layout(p, include = c("seq", "ribbon", "axis"))
  expect_s3_class(exported, "ggchord_layout_export")
  expect_true(all(c("layer_id", "source_row") %in% names(exported$ribbon)))
  expect_true(nrow(exported$axis) > 0L)

  preview <- view_ggchord(
    p, width = 2, height = 2, units = "in", dpi = 50,
    device = "png", viewer = "none"
  )
  expect_true(file.exists(preview))
  expect_true(file.exists(sub("\\.png$", ".html", preview)))
  expect_equal(formals(view_ggchord)$width, 11)
  expect_null(formals(view_ggchord)$height)
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


test_that("ribbon compatibility name preserves geometry and explicit styles", {
  data(seq_data_example)
  data(ribbon_data_example)
  base <- ggchord(seq_data_example, ribbon_data_example) + geom_seq()
  expect_warning(
    legacy <- geom_ribbon(fill = "red", color = "blue", ribbon_gap = 0.08),
    "deprecated; use geom_link_ribbon"
  )
  canonical <- geom_link_ribbon(fill = "red", color = "blue", ribbon_gap = 0.08)
  expect_equal(legacy$ggchord_params, canonical$ggchord_params)
  expect_equal(
    get_chord_layout(base + legacy)$ribbon_polys,
    get_chord_layout(base + canonical)$ribbon_polys
  )
  expect_no_warning(stat_ribbon_bundle())
  expect_no_warning(stat_ribbon_density())
})
