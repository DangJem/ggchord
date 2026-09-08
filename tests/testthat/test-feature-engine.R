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

test_that("short arrows fall back without swallowing the interval", {
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
  expect_equal(nrow(auto), 120L)
  expect_true(all(is.finite(auto$x) & is.finite(auto$y)))

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
    accver = "circle", start = c(50, 250, 450, 650),
    end = c(180, 380, 580, 780), strand = c("+", "-", "+", "-"),
    type = c("CDS", "RBS", "terminator", "protein_bind")
  )
  p <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature(
      aes(feature_shape = type, feature_fill = type), data = data,
      position = "plasmid"
    ) + scale_feature_shape_plasmid() + scale_feature_fill_plasmid() +
    coord_circular()
  out <- export_ggchord_layout(p, include = "feature")$feature
  expect_equal(sort(unique(out$source_row)), seq_len(nrow(data)))
  expect_setequal(unique(out$feature_shape),
    c("arrow", "chevron", "lollipop", "block"))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})

test_that("plasmid feature preset keeps one arrowhead and draws segment joins", {
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
  p <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature_plasmid(data = feature) + coord_circular()
  out <- export_ggchord_layout(p, include = "feature")$feature
  expect_equal(unique(out$.component), c("polygon", "boundary"))
  expect_equal(length(unique(out$group[out$.component == "polygon"])), 1L)
  expect_equal(nrow(out[out$.component == "boundary", ]), 2L)
  expect_s3_class(ggplot2::ggplotGrob(p), "gtable")

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
  expect_true(any(moved$.component == "segment"))
})

test_that("plasmid feature preset stacks overlapping intervals by default", {
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
  expect_setequal(unique(out$lane), 0:1)
  expect_true(all(out$lane_offset %in% c(0, -.10)))

  grouped <- transform(features, feature_group = "related")
  grouped_out <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = grouped) + coord_circular(),
    include = "feature"
  )$feature
  expect_equal(unique(grouped_out$lane), 0L)

  preferred <- transform(features, preferred_lane = c(0L, 2L))
  preferred_out <- export_ggchord_layout(
    ggchord(seq, validate = "none") + geom_seq() +
      geom_feature_plasmid(data = preferred) + coord_circular(),
    include = "feature"
  )$feature
  expect_setequal(unique(preferred_out$lane), c(0L, 2L))
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
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
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
