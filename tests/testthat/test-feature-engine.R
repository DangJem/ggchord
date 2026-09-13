test_that("gene and arrow feature use one geometry engine", {
  seq <- data.frame(accver = "circle", length = 1000)
  data <- data.frame(
    accver = "circle", start = c(950, 250, 450), end = c(50, 350, 451),
    strand = c("+", "-", "+"), anno = c("cross", "minus", "short")
  )
  common <- function(layer, category) {
    export_ggchord_layout(
      ggchord(seq, validate = "none") + geom_seq() + layer + coord_circular(),
      include = category
    )[[category]]
  }
  gene <- common(geom_gene(
    data = data, position = "plasmid", arrow_head_length = .04,
    arrow_head_width = 1.25, short_feature = "auto"
  ), "gene")
  feature <- common(geom_feature(
    data = data, position = "plasmid", feature_shape = "arrow",
    arrow_head_length = .04, arrow_head_width = 1.25,
    short_feature = "auto"
  ), "feature")
  for (column in c(
    "x", "y", "group", "source_row", "position_name", "base_offset",
    "lane", "lane_offset", "normal_offset"
  )) expect_equal(gene[[column]], feature[[column]])
  expect_gt(length(unique(gene$group[gene$source_row == 1L])), 1L)
  expect_true(all(is.finite(gene$x) & is.finite(gene$y)))
  expect_s3_class(geom_feature(data = data, shape = "arrow"), "LayerInstance")
  expect_error(
    geom_feature(data = data, shape = "arrow", feature_shape = "block"),
    "only one"
  )
})

test_that("short arrows preserve body shoulders and biological midpoint", {
  seq <- data.frame(accver = "circle", length = 1000)
  data <- data.frame(
    accver = "circle", start = 100, end = 101, strand = "+", anno = "short"
  )
  draw <- function(fallback) export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature(data = data, short_feature = fallback) + coord_circular(),
    include = "feature"
  )$feature
  auto <- draw("auto")
  wedge <- draw("wedge")
  block <- draw("block")
  expect_gt(nrow(auto), 3L)
  expect_true(all(is.finite(auto$x) & is.finite(auto$y)))
  display <- ggchord_arrow_display_interval(0, 2 * pi / 1000, .9, .07, .04)
  expect_gt(diff(display), 2 * pi / 1000)
  expect_equal(mean(display), pi / 1000)
  auto_polygon <- auto[auto$.component == "polygon", , drop = FALSE]
  auto_radius <- sqrt(auto_polygon$x^2 + auto_polygon$y^2)
  expect_gte(length(unique(round(auto_radius, 4))), 3L)

  medium <- data
  medium$end <- 115
  medium_auto <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature(data = medium, short_feature = "auto") +
      coord_circular(), include = "feature"
  )$feature
  medium_auto <- medium_auto[medium_auto$.component == "polygon", ]
  expect_gt(nrow(medium_auto), 3L)
  # A shortened shouldered arrow retains a body and multiple radial levels;
  # it must not collapse to the three-vertex wedge fallback.
  radius <- sqrt(medium_auto$x^2 + medium_auto$y^2)
  expect_gte(length(unique(round(radius, 4))), 3L)
  expect_true(nrow(wedge) >= 3L)
  expect_equal(nrow(block), 120L)
})

test_that("all feature shape factories build and preserve source identity", {
  seq <- data.frame(accver = "circle", length = 1000)
  data <- data.frame(
    accver = "circle", start = seq(40, 740, length.out = 8),
    end = seq(110, 810, length.out = 8), strand = rep(c("+", "-"), 4),
    type = c("CDS", "regulatory", "promoter", "primer", "marker",
      "protein_bind", "RBS", "terminator")
  )
  shapes <- c(CDS = "arrow", regulatory = "compact_arrow",
    promoter = "promoter_arrow", primer = "primer_arrow", marker = "marker",
    protein_bind = "block", RBS = "chevron", terminator = "lollipop")
  p <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature(
      aes(feature_shape = type, feature_fill = type), data = data,
      position = "plasmid"
    ) + scale_feature_shape_manual(values = shapes) +
    scale_feature_fill_plasmid() +
    coord_circular()
  out <- export_ggchord_layout(p, include = "feature")$feature
  expect_equal(sort(unique(out$source_row)), seq_len(nrow(data)))
  expect_setequal(unique(out$feature_shape), unname(shapes))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("directional variants share topology with type width ratios", {
  seq <- data.frame(accver = "circle", length = 1000)
  feature <- data.frame(
    accver = "circle", start = c(50, 300, 550, 800),
    end = c(130, 380, 630, 880), strand = "+",
    anno = letters[1:4],
    feature_shape = c(
      "arrow", "compact_arrow", "promoter_arrow", "primer_arrow"
    )
  )
  out <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature(data = feature, position = "plasmid", feature_width = .07) +
      coord_circular(),
    include = "feature"
  )$feature
  out <- out[out$.component == "polygon", , drop = FALSE]
  radial_extent <- lapply(split(out, out$source_row), function(x) {
    radius <- sqrt(x$x^2 + x$y^2)
    c(vertices = nrow(x), width = diff(range(radius)))
  })
  extent <- do.call(rbind, radial_extent)
  expect_equal(length(unique(extent[, "vertices"])), 1L)
  expect_equal(
    unname(extent[, "width"] / extent[1L, "width"]),
    c(1, .78, .72, .65), tolerance = 1e-5
  )
})

test_that("feature label fitting uses final font metrics and exports modes", {
  seq <- data.frame(accver = "circle", length = 1000)
  feature <- data.frame(
    accver = "circle", start = 100, end = 220, strand = "+",
    anno = "measured label", type = "CDS"
  )
  draw <- function(size) export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = feature) +
      geom_feature_label_repel(
        data = feature, size = size, family = "mono",
        fontface = "bold", lineheight = 1.35
      ) + coord_circular(), include = "labels"
  )$labels
  small <- draw(2)
  large <- draw(8)
  small <- small[small$.component == "text", , drop = FALSE]
  large <- large[large$.component == "text", , drop = FALSE]
  expect_equal(small$feature_label_mode, "inside")
  expect_true(large$feature_label_mode %in% c("adjacent", "external"))
  expect_equal(large$family, "mono")
  expect_equal(large$fontface, "bold")
  expect_equal(large$lineheight, 1.35)
})

test_that("inside feature text contrast follows fixed and mapped fill", {
  seq <- data.frame(accver = "circle", length = 1000)
  feature <- data.frame(
    accver = "circle", start = c(100, 400), end = c(300, 600),
    strand = "+", anno = c("dark", "light"), fill_role = c("a", "b")
  )
  tracks <- position_feature_stack(base_position = position_plasmid())
  fixed <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(
        data = feature[1, ], feature_fill = "#202020", position = tracks
      ) +
      geom_feature_label_repel(data = feature[1, ], position = tracks) +
      coord_circular(),
    include = "labels"
  )$labels
  expect_equal(unique(fixed$feature_label_colour), "#FFFFFF")

  mapped <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(
        aes(feature_fill = fill_role), data = feature, position = tracks
      ) +
      scale_feature_fill_manual(values = c(a = "#202020", b = "#F4F4F4")) +
      geom_feature_label_repel(data = feature, position = tracks) +
      coord_circular(),
    include = "labels"
  )$labels
  mapped <- mapped[mapped$.component == "text", , drop = FALSE]
  expect_equal(mapped$feature_label_colour, c("#FFFFFF", "#202020"))
})

test_that("dense plasmid feature labels avoid labels and feature glyphs", {
  data(plasmid_example_pBluescript_II_SK_plus)
  features <- find_common_features(plasmid_example_pBluescript_II_SK_plus)
  tracks <- position_feature_stack(
    spacing = .10, base_position = position_plasmid()
  )
  plot <- ggchord(plasmid_example_pBluescript_II_SK_plus) + geom_seq() +
    geom_feature_plasmid(data = features, position = tracks) +
    geom_feature_label_repel(data = features, position = tracks) +
    coord_circular(rotation = 90)
  first <- export_ggchord_layout(
    plot, include = c("feature", "labels")
  )
  second <- export_ggchord_layout(
    plot, include = c("feature", "labels")
  )
  expect_identical(first$labels, second$labels)
  labels <- first$labels[first$labels$.component == "text", , drop = FALSE]
  expect_true(all(labels$feature_label_mode %in%
    c("inside", "adjacent", "external")))
  expect_true(any(labels$feature_label_mode == "inside"))
  expect_true(any(labels$feature_label_mode != "inside"))
  expect_false(any(labels$feature_label_mode == "external"))
  expect_true(any(labels$label_track > 1L))
  feature_lanes <- unique(first$feature[
    first$feature$.component == "polygon", c("anno", "lane")
  ])
  lane_of <- stats::setNames(feature_lanes$lane, feature_lanes$anno)
  expect_equal(unname(lane_of[c(
    "lacZα", "lac operator", "lac promoter", "ori", "AmpR",
    "AmpR promoter"
  )]), rep(0L, 6L))
  expect_equal(unname(lane_of[c(
    "f1 ori", "M13 fwd", "T7 promoter", "MCS", "T3 promoter", "M13 rev"
  )]), rep(1L, 6L))
  expect_equal(unname(lane_of[c("KS primer", "SK primer")]), c(2L, 2L))
  polygon <- first$feature[first$feature$.component == "polygon", ]
  polygon$radius <- sqrt(polygon$x^2 + polygon$y^2)
  band_bounds <- lapply(split(polygon, polygon$lane), function(x) {
    range(x$radius)
  })
  expect_gt(band_bounds[["0"]][1L] - band_bounds[["1"]][2L], .10)
  expect_gt(band_bounds[["1"]][1L] - band_bounds[["2"]][2L], .015)
  expect_lt(band_bounds[["1"]][1L] - band_bounds[["2"]][2L], .06)
  adjacent <- labels[labels$feature_label_mode == "adjacent", , drop = FALSE]
  label_radius <- sqrt(adjacent$x^2 + adjacent$y^2)
  # One label-track id means one physical circle, rather than a collection of
  # per-label offsets that merely happen to be called tracks.
  radius_spread <- tapply(label_radius, adjacent$label_track,
    function(value) diff(range(value)))
  expect_true(all(radius_spread < 1e-8))
  tangent_angle <- atan2(adjacent$y, adjacent$x) * 180 / pi + 90
  tangent_error <- abs((adjacent$angle - tangent_angle + 90) %% 180 - 90)
  expect_true(all(tangent_error < 1e-6))
  expect_true(any(first$labels$.component == "segment"))
  leader_features <- unique(first$labels$anno[
    first$labels$.component == "segment"
  ])
  expect_true(all(c(
    "M13 fwd", "T7 promoter", "KS primer", "SK primer",
    "T3 promoter", "M13 rev", "lac operator"
  ) %in% leader_features))
  expect_true(all(labels$.draw_as_arc))
  glyphs <- first$labels[
    first$labels$.component == "arc_text", , drop = FALSE
  ]
  expect_equal(nrow(glyphs), sum(nchar(labels$label, type = "chars")))
  glyph_radius <- sqrt(glyphs$x^2 + glyphs$y^2)
  radius_by_label <- tapply(glyph_radius, glyphs$.arc_parent_group,
    function(value) diff(range(value)))
  expect_true(all(radius_by_label < 1e-8))
  glyph_tangent <- atan2(glyphs$y, glyphs$x) * 180 / pi + 90
  glyph_tangent_error <- abs((glyphs$angle - glyph_tangent + 90) %% 180 - 90)
  expect_true(all(glyph_tangent_error < 1e-6))
  glyph_boxes <- ggchord_text_boxes(
    glyphs, x_col = "x", y_col = "y", text_col = "label",
    angle_col = "angle", size_col = "size", units_per_inch = .25,
    box_padding = .005
  )
  cross_label_overlaps <- 0L
  if (nrow(glyph_boxes) > 1L) for (i in seq_len(nrow(glyph_boxes) - 1L)) {
    other <- seq.int(i + 1L, nrow(glyph_boxes))
    other <- other[glyphs$.arc_parent_group[other] !=
      glyphs$.arc_parent_group[i]]
    if (length(other)) cross_label_overlaps <- cross_label_overlaps + sum(
      ggchord_oriented_box_overlaps(
        glyph_boxes[i, , drop = FALSE], glyph_boxes[other, , drop = FALSE]
      )
    )
  }
  expect_equal(cross_label_overlaps, 0L)
  long_glyphs <- glyphs[glyphs$label != " " & glyphs$source_row ==
    labels$source_row[labels$label == "AmpR promoter"], , drop = FALSE]
  expect_gt(diff(range(long_glyphs$angle)), 1)
  internal <- get_chord_layout(plot)
  internal_labels <- internal$gene_labels
  polygons <- split(
    internal$gene_polys[
      internal$gene_polys$.component == "polygon", , drop = FALSE
    ],
    internal$gene_polys$group[
      internal$gene_polys$.component == "polygon"
    ]
  )
  for (i in seq_len(nrow(internal_labels))) {
    internal_boxes <- ggchord_arc_text_layout(
      internal_labels[i, , drop = FALSE],
      units_per_inch = internal$text_units_per_inch, box_padding = .04
    )$boxes
    expect_false(ggchord_boxes_hit_features(
      internal_boxes, polygons,
      source_row = internal_labels$source_row[i],
      allow_own = internal_labels$feature_label_mode[i] == "inside"
    ))
  }
  expect_true(all(labels$angle <= 90 | labels$angle >= 270))
  expect_true("feature_class" %in% names(features))
  expect_false("preferred_lane" %in% names(features))
  expected_shape <- ifelse(features$strand == ".", "block", "arrow")
  expected_shape[features$type == "promoter" & features$strand != "."] <-
    "promoter_arrow"
  expected_shape[features$type == "primer_bind" & features$strand != "."] <-
    "primer_arrow"
  expected_shape[features$type == "protein_bind" & features$strand != "."] <-
    "compact_arrow"
  expect_equal(features$feature_shape, expected_shape)
  boundary_features <- table(first$feature$anno[
    first$feature$.component == "boundary"
  ])
  expect_true(all(c("AmpR", "lac promoter") %in% names(boundary_features)))
  expect_equal(unname(boundary_features[c("AmpR", "lac promoter")]), c(8L, 16L))
})

test_that("plasmid feature preset merges segment joins and cleavage marks", {
  seq <- data.frame(accver = "circle", length = 1000)
  segments <- data.frame(
    segment_index = 1:2, segment_type = "standard",
    start = c(100, 181), end = c(180, 300), color = "#CCFFCC"
  )
  feature <- data.frame(
    accver = "circle", start = 100, end = 300, strand = "+",
    anno = "joined", feature_color = "#CCFFCC",
    feature_shape = "arrow", segments = I(list(segments))
  )
  feature$cleavage_arrows <- I(list(180))
  p <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_plasmid(data = feature) + coord_circular()
  out <- export_ggchord_layout(p, include = "feature")$feature
  expect_equal(unique(out$.component), c("polygon", "boundary"))
  expect_equal(length(unique(out$group[out$.component == "polygon"])), 1L)
  boundary <- out[out$.component == "boundary", ]
  expect_equal(length(unique(boundary$group)), 4L)
  expect_equal(nrow(boundary), 8L)
  expect_equal(unique(boundary$boundary_linetype), "dotted")
  expect_equal(unique(boundary$boundary_draw_linetype), "solid")
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")

  feature$cleavage_arrows <- I(list(numeric()))
  joined_only <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = feature) + coord_circular(),
    include = "feature"
  )$feature
  joined_boundary <- joined_only[joined_only$.component == "boundary", ]
  expect_equal(length(unique(joined_boundary$group)), 4L)
  expect_equal(unique(joined_boundary$boundary_linetype), "dotted")

  styled_segments <- segments
  styled_segments$line_style <- c("solid", "dashed")
  feature$segments <- I(list(styled_segments))
  styled <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = feature) + coord_circular(),
    include = "feature"
  )$feature
  expect_equal(unique(stats::na.omit(styled$boundary_linetype)), "dashed")

  no_join <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_plasmid(data = feature, segment_boundaries = FALSE) +
    coord_circular()
  expect_s3_class(ggplot2::ggplotGrob(no_join), "gtable")
  expect_s3_class(geom_feature(data = feature, arrow_head_style = "flush"),
    "LayerInstance")
  expect_s3_class(geom_feature(data = feature, arrow_head_style = "triangle"),
    "LayerInstance")
})

test_that("feature labels reuse fixed and repel label layout", {
  seq <- data.frame(accver = "circle", length = 1000)
  features <- data.frame(
    accver = "circle", start = seq(100, 170, by = 10),
    end = seq(106, 176, by = 10), strand = rep(c("+", "-"), 4),
    type = rep(c("CDS", "promoter"), 4),
    anno = paste("dense feature", seq_len(8))
  )
  fixed_plot <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_label(
      aes(label = anno), data = features, position = "plasmid",
      label_orientation = "tangent", label_side = "outside",
      label_overlap = "allow"
    ) + coord_circular()
  fixed <- export_ggchord_layout(fixed_plot, include = "labels")$labels
  expect_true(any(grepl("dense feature", fixed$text)))
  expect_false(any(fixed$text %in% c("CDS", "promoter")))
  expect_equal(unique(fixed$normal_offset), -.1)

  repel_plot <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_label_repel(
      aes(label = anno), data = features, position = "plasmid",
      label_layout = "callout"
    ) + coord_circular()
  first <- export_ggchord_layout(repel_plot, include = "labels")$labels
  second <- export_ggchord_layout(repel_plot, include = "labels")$labels
  expect_identical(first, second)
  expect_true(any(first$label_layout == "callout", na.rm = TRUE))
  expect_s3_class(ggplot2::ggplotGrob(repel_plot), "gtable")

  feature_plot <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_label_repel(
      aes(label = anno), data = features, position = "plasmid"
    ) + coord_circular()
  staged <- export_ggchord_layout(
    feature_plot, include = "labels"
  )$labels
  expect_true(any(staged$.component == "text"))
  expect_true(all(staged$label_layout[staged$.component == "text"] ==
    "feature"))
  expect_s3_class(ggplot2::ggplotGrob(feature_plot), "gtable")

  displaced <- features[rep(1L, 2L), , drop = FALSE]
  displaced$start <- 100
  displaced$end <- 110
  displaced$anno <- c("long feature alpha", "long feature beta")
  displaced_plot <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_label_repel(data = displaced, position = "plasmid") +
    coord_circular()
  moved <- export_ggchord_layout(
    displaced_plot, include = "labels"
  )$labels
  moved_text <- moved[moved$.component == "text", , drop = FALSE]
  expect_true(all(moved_text$feature_label_mode %in%
    c("adjacent", "external")))
  if (any(moved_text$feature_label_mode == "external")) {
    expect_true(any(moved$.component == "segment"))
  }
})

test_that("plasmid tracks use priority, span and compact interval reuse", {
  seq <- data.frame(accver = "circle", length = 1000)
  features <- data.frame(
    accver = "circle", start = c(100, 150), end = c(300, 220),
    strand = c("+", "+"), anno = c("outer", "nested")
  )
  out <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = features) + coord_circular(),
    include = "feature"
  )$feature
  lanes <- unique(out[, c("anno", "lane")])
  expect_equal(lanes$lane[lanes$anno == "outer"], 0L)
  expect_equal(lanes$lane[lanes$anno == "nested"], 1L)

  typed <- transform(features, type = c("promoter", "CDS"))
  typed_out <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = typed) + coord_circular(),
    include = "feature"
  )$feature
  expect_equal(unique(typed_out[, c("anno", "lane")])$lane, lanes$lane)

  preferred <- transform(features, display_priority = c(0, 10))
  preferred_out <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = preferred) + coord_circular(),
    include = "feature"
  )$feature
  preferred_lanes <- unique(preferred_out[, c("anno", "lane")])
  expect_equal(preferred_lanes$lane[preferred_lanes$anno == "nested"], 0L)
  expect_equal(preferred_lanes$lane[preferred_lanes$anno == "outer"], 1L)

  chain <- data.frame(
    accver = "circle", start = c(100, 220, 340), end = c(240, 360, 470),
    strand = "+", anno = c("A", "B", "C")
  )
  chain_out <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = chain) + coord_circular(),
    include = "feature"
  )$feature
  chain_lanes <- unique(chain_out[, c("anno", "lane")])
  expect_equal(chain_lanes$lane[chain_lanes$anno %in% c("A", "C")], c(0L, 0L))
  expect_equal(chain_lanes$lane[chain_lanes$anno == "B"], 1L)
})

test_that("feature directions and segment-specific styles stay semantic", {
  seq <- data.frame(accver = "circle", length = 1000)
  segments <- data.frame(
    segment_index = 1:3, segment_type = "standard",
    start = c(100, 151, 201), end = c(150, 200, 300),
    color = c("#FF0000", "#FF0000", "#00FF00"),
    line_style = c("solid", "dotted", "dashed")
  )
  feature <- data.frame(
    accver = "circle", start = 100, end = 300,
    direction = "bidirectional", anno = "split",
    feature_color = "#0000FF"
  )
  feature$segments <- I(list(segments))
  plot <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature(data = feature) + coord_circular()
  out <- export_ggchord_layout(plot, include = "feature")$feature
  polygons <- out[out$.component == "polygon", , drop = FALSE]
  expect_true(all(polygons$biological_strand == "+/-"))
  expect_equal(unique(polygons$feature_fill_explicit),
    c("#FF0000", "#00FF00"))
  expect_equal(unique(stats::na.omit(out$boundary_linetype)), "dotted")
  boundary <- out[out$.component == "boundary", , drop = FALSE]
  expect_equal(length(unique(boundary$group)), 4L)
  expect_true(all(boundary$boundary_draw_linetype == "solid"))
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
})

test_that("point and multi-segment features retain distinct layout semantics", {
  seq <- data.frame(accver = "circle", length = 1000)
  point <- data.frame(
    accver = "circle", start = 80, end = 80,
    directionality = "nondirectional", anno = "point"
  )
  point_out <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = point) + coord_circular(),
    include = "feature"
  )$feature
  expect_equal(unique(point_out$.component), "point")
  expect_equal(nrow(point_out), 2L)
  expect_equal(unique(point_out$feature_shape), "point")

  segments <- data.frame(
    segment_index = 1:3, segment_type = c("standard", "gap", "standard"),
    start = c(100, 181, 221), end = c(180, 220, 320),
    color = c("#FF0000", NA, "#00AA00"),
    segment_name = c("left", NA, "right")
  )
  multi <- data.frame(
    accver = "circle", start = 100, end = 320,
    directionality = "forward", anno = "one feature",
    segments = I(list(segments))
  )
  p <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_plasmid(data = multi) +
    geom_feature_label_repel(data = multi) + coord_circular()
  layout <- export_ggchord_layout(p, include = c("feature", "labels"))
  polygons <- layout$feature[layout$feature$.component == "polygon", ]
  expect_equal(unique(polygons$source_row), 1L)
  expect_equal(unique(polygons$lane), 0L)
  expect_setequal(unique(polygons$feature_fill_explicit),
    c("#FF0000", "#00AA00"))
  expect_equal(sum(layout$labels$.component == "text"), 1L)
})

test_that("feature labels use external fallback or hide by policy", {
  seq <- data.frame(accver = "circle", length = 1000)
  feature <- data.frame(
    accver = "circle", start = rep(100, 3), end = rep(105, 3),
    strand = "+", anno = paste("very long feature label", 1:3)
  )
  draw <- function(external) export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = feature) +
      geom_feature_label_repel(data = feature, external = external,
        max_overlaps = 0) + coord_circular(),
    include = "labels"
  )$labels
  external <- draw(TRUE)
  external <- external[external$.component == "text", , drop = FALSE]
  internal_only <- draw(FALSE)
  internal_only <- internal_only[
    internal_only$.component == "text", , drop = FALSE
  ]
  expect_true(any(external$feature_label_mode == "external"))
  expect_lte(nrow(internal_only), nrow(external))
  expect_false(any(internal_only$feature_label_mode == "external"))
})

test_that("plasmid preset scales defer to manual scales in either order", {
  seq <- data.frame(accver = "circle", length = 100)
  base <- ggchord(seq, validate = "none")
  manual_fill <- scale_feature_fill_manual(values = c(CDS = "black"))
  manual_shape <- scale_feature_shape_manual(values = c(CDS = "block"))
  first <- base + scale_feature_fill_plasmid() + manual_fill +
    scale_feature_shape_plasmid() + manual_shape
  second <- base + manual_fill + scale_feature_fill_plasmid() +
    manual_shape + scale_feature_shape_plasmid()
  map_value <- function(plot, aesthetic, value) {
    scale <- plot$scales$get_scales(aesthetic)$clone()
    scale$train(value)
    unname(scale$map(value))
  }
  expect_equal(map_value(first, "feature_fill", "CDS"), "black")
  expect_equal(map_value(second, "feature_fill", "CDS"), "black")
  expect_equal(map_value(first, "feature_shape", "CDS"), "block")
  expect_equal(map_value(second, "feature_shape", "CDS"), "block")
  expect_s3_class(theme_ggchord_plasmid(), "theme")
})

test_that("circular feature plot renders to standard vector and raster devices", {
  seq <- data.frame(accver = "circle", label = "Example", length = 1000)
  data <- data.frame(
    accver = "circle", start = c(100, 400), end = c(250, 550),
    strand = c("+", "-"), type = c("CDS", "promoter")
  )
  p <- ggchord(seq, validate = "none") +
    geom_seq(seq_style = "double") +
    geom_feature(aes(feature_fill = type), data = data, position = "plasmid") +
    geom_seq_center_label() + scale_feature_fill_plasmid() +
    coord_circular(rotation = 90) + theme_ggchord_plasmid()
  png_file <- tempfile(fileext = ".png")
  pdf_file <- tempfile(fileext = ".pdf")
  ggplot2::ggsave(png_file, p, width = 4, height = 4, dpi = 72)
  ggplot2::ggsave(pdf_file, p, width = 4, height = 4)
  expect_gt(file.info(png_file)$size, 0)
  expect_gt(file.info(pdf_file)$size, 0)
  if (requireNamespace("svglite", quietly = TRUE)) {
    svg_file <- tempfile(fileext = ".svg")
    ggplot2::ggsave(svg_file, p, width = 4, height = 4)
    expect_gt(file.info(svg_file)$size, 0)
  }
})
