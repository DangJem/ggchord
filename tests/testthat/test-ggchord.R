# testthat test file - v0.4.0 layered API
#
# Slow rendering / integration tests are opt-in: set GGCHORD_RUN_SLOW_TESTS=1.

test_that("ggchord builds with minimal and layered inputs", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p0 <- ggchord(seq_data = seq_data_example)
  expect_s3_class(p0, "ggchord")
  expect_s3_class(p0, "ggplot")

  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis()
  expect_s3_class(p, "ggchord")
  expect_gte(length(p$layers), 4)
})

test_that("parameters can be distributed across geoms", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p <- ggchord(
    seq_data = seq_data_example,
    ribbon_data = ribbon_data_example,
    gene_data = gene_data_example,
    title = "Test",
    rotation = 30
  ) +
    geom_seq(
      seq_radius = c(3, 2, 2, 1),
      seq_curvature = c(0, 1, -1, 1.5),
      seq_orientation = c(-1, -1, -1, 1)
    ) +
    geom_ribbon(
      ribbon_color_scheme = "query",
      ribbon_alpha = 0.5
    ) +
    geom_gene(
      gene_color_scheme = "strand",
      gene_width = 0.08
    ) +
    geom_axis(
      axis_gap = 0.02,
      axis_label_orientation = c(0, 45, 80, 130)
    )

  expect_s3_class(p, "ggchord")
})

test_that("seq_data missing required columns errors", {
  bad_data <- data.frame(id = c("a", "b"), len = c(100, 200))
  expect_error(
    ggchord(seq_data = bad_data),
    "seq_data"
  )
})

test_that("seq_data with non-positive lengths errors", {
  bad_data <- data.frame(seq_id = c("a", "b"), length = c(0, 200))
  expect_error(
    ggchord(seq_data = bad_data),
    "positive"
  )
})

test_that("missing and non-finite input values are rejected before layout", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  expect_error(ggchord(transform(seq_data_example, length = NA_real_)),
               "finite positive")
  expect_error(ggchord(transform(seq_data_example, length = Inf)),
               "finite positive")
  expect_error(
    ggchord(seq_data_example, transform(ribbon_data_example, qstart = NA_real_)),
    "numeric columns"
  )
  expect_error(
    ggchord(seq_data_example, gene_data = transform(gene_data_example, start = Inf)),
    "finite numbers"
  )
})

test_that("invalid numeric layout parameters report their parameter names", {
  data(seq_data_example)
  expect_error(
    ggplot_build(ggchord(seq_data_example) + geom_seq(seq_radius = NA_real_)),
    "seq_radius"
  )
  expect_error(
    ggplot_build(ggchord(seq_data_example) + geom_seq() +
                   geom_axis(axis_tick_minor_number = -1)),
    "axis_tick_minor_number"
  )
  expect_error(
    ggplot_build(ggchord(seq_data_example) + geom_seq() +
                   geom_seq_label(seq_label_size = Inf)),
    "seq_label_size"
  )
})

test_that("print renders the full chord diagram", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  # pident works correctly together with the gene fill scale
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_color_scheme = "query") +
    geom_gene() +
    geom_axis()

  # Render to PDF (use a session temp file so the tests pass on any platform
  # and do not leave artifacts behind for R CMD check)
  out <- tempfile(fileext = ".pdf")
  pdf(out, 8, 8)
  suppressMessages(suppressWarnings(print(p)))
  dev.off()
  expect_true(file.exists(out))
})

test_that("README color, label override, and transparency parameters take effect", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq(linewidth = 2) +
    geom_ribbon(ribbon_color_scheme = "query", alpha = 0.2) +
    geom_gene(gene_color_scheme = "manual") +
    geom_gene_label(gene_label_size = 4) +
    geom_axis()

  ggplot_build(p)
  layout <- ggchord:::get_chord_layout()
  expect_true("fill" %in% names(layout$ribbon_polys))
  expect_true(all(layout$ribbon_polys$alpha == 0.2))
  expect_gt(nrow(layout$gene_labels), 0)
})

test_that("ribbon outline parameters work with sensible defaults", {
  data(seq_data_example)
  data(ribbon_data_example)
  # Defaults: black outline, width 0.05, solid line
  p <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
  l <- p$layers[[2]]
  expect_equal(l$aes_params$colour, "black")
  expect_equal(l$aes_params$linewidth, 0.05)
  expect_equal(l$aes_params$linetype, 1)

  # Custom values
  p2 <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_outline_color = "red", ribbon_outline_width = 0.8,
                ribbon_outline_linetype = "dashed")
  l2 <- p2$layers[[2]]
  expect_equal(l2$aes_params$colour, "red")
  expect_equal(l2$aes_params$linewidth, 0.8)
  expect_equal(l2$aes_params$linetype, "dashed")

  expect_no_error(ggplot_build(p2))
})

test_that("plot objects are self-contained (no cross-talk between plots)", {
  data(seq_data_example)
  p1 <- ggchord(seq_data_example) + geom_seq(seq_radius = 5)
  p2 <- ggchord(seq_data_example) + geom_seq()
  expect_equal(p1$layers[[1]]$ggchord_params$seq_radius, 5)
  expect_null(p2$layers[[1]]$ggchord_params$seq_radius)
  # both build independently
  expect_s3_class(ggplot_build(p1), "ggplot_built")
  expect_s3_class(ggplot_build(p2), "ggplot_built")
})

test_that("plots survive saveRDS/readRDS and render", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq(seq_radius = 5)
  f <- tempfile(fileext = ".rds")
  saveRDS(p, f)
  p2 <- readRDS(f)
  expect_no_error(ggplot_build(p2))
  expect_equal(p2$layers[[1]]$ggchord_params$seq_radius, 5)
})



test_that("plot objects are self-contained and repeated builds stay stable", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis()
  # Layout/scales are computed lazily; prepare explicitly to verify the plot
  # carries its own scales (tagged ggchord_managed) so that tools such as
  # plotly::ggplotly() that clone the plot see the same scales as a build.
  p <- ggchord:::prepare_ggchord_plot(p)
  expect_gt(length(p$scales$scales), 0)
  expect_true(all(vapply(p$scales$scales,
                         function(s) !is.null(s$ggchord_managed), logical(1))))
  expect_false(is.null(p$layers[[1]]$ggchord_params))
  # repeated builds remain stable
  expect_no_error(ggplot_build(p))
})

test_that("geom_seq_label places sequence labels", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) +
    geom_seq() +
    geom_seq_label(seq_label_radius = 1.25, seq_label_size = 3.5)
  ggplot_build(p)
  layout <- ggchord:::get_chord_layout()
  expect_gt(nrow(layout$seq_labels_df), 0)
  expect_true(all(c("text_x", "text_y", "label") %in% names(layout$seq_labels_df)))
})

test_that("seq_label_radius places labels outside/inside the arcs", {
  data(seq_data_example)
  label_multiplier <- function(radius) {
    p <- ggchord(seq_data_example) + geom_seq() +
      geom_seq_label(seq_label_radius = radius)
    ggplot_build(p)
    l <- get_chord_layout()
    arc_r <- vapply(l$seq_arcs, function(a) {
      median(sqrt(a$x^2 + a$y^2))
    }, numeric(1))
    lbl_r <- sqrt(l$seq_labels_df$text_x^2 + l$seq_labels_df$text_y^2)
    lbl_r / arc_r
  }
  # > 1 = outside, 1 = on the arc, < 1 = inside (as documented)
  expect_true(all(label_multiplier(1.2) > 1))
  expect_true(all(abs(label_multiplier(1) - 1) < 1e-6))
  expect_true(all(label_multiplier(0.8) < 1))
})

test_that("sequence labels stay readable under global rotation", {
  data(seq_data_example)
  data(gene_data_example)
  for (rot in c(0, 45, 90, 135)) {
    p <- ggchord(seq_data_example, rotation = rot) + geom_seq() + geom_seq_label()
    ggplot_build(p)
    a <- get_chord_layout()$seq_labels_df$text_angle
    expect_equal(sum(a > 90 & a < 270), 0,
                 info = paste("seq labels upside down at rotation", rot))
  }
  # the same readability fix applies to fixed gene labels
  p <- ggchord(seq_data_example, gene_data = gene_data_example, rotation = 90) +
    geom_seq() + geom_gene() + geom_gene_label()
  ggplot_build(p)
  a <- get_chord_layout()$gene_labels$text_angle
  expect_equal(sum(a > 90 & a < 270), 0)
})

test_that("seq_label_orientation horizontal draws horizontal labels", {
  data(seq_data_example)
  p <- ggchord(seq_data_example, rotation = 45) + geom_seq() +
    geom_seq_label(seq_label_orientation = "horizontal")
  ggplot_build(p)
  l <- get_chord_layout()
  expect_true(all(l$seq_labels_df$text_angle == 0))
  # text extends away from the chord center: hjust 0 on the right, 1 on the left
  expect_equal(l$seq_labels_df$hjust, ifelse(l$seq_labels_df$text_x >= 0, 0, 1))
  # invalid orientation is rejected
  expect_error(geom_seq_label(seq_label_orientation = "vertical"), "should be")
})

test_that("seq_label_hjust and seq_label_vjust are applied", {
  data(seq_data_example)
  # rotation 0 avoids the readability flips that toggle hjust
  p <- ggchord(seq_data_example, rotation = 0) + geom_seq() +
    geom_seq_label(seq_label_hjust = c(0.1, 0.9, 0.1, 0.9), seq_label_vjust = 1)
  ggplot_build(p)
  l <- get_chord_layout()
  # user-supplied hjust values are used as-is; when the readability flip turns
  # a label by 180 degrees the justification toggles (1 - h) so the text box
  # stays anchored at the arc midpoint
  expect_equal(unname(l$seq_labels_df$hjust), c(0.9, 0.9, 0.1, 0.1),
               tolerance = 1e-6)
  expect_equal(unique(l$seq_labels_df$vjust), 1)
  # a centered justification is flip-invariant
  p2 <- ggchord(seq_data_example, rotation = 0) + geom_seq() +
    geom_seq_label(seq_label_hjust = 0.5)
  ggplot_build(p2)
  expect_equal(unique(get_chord_layout()$seq_labels_df$hjust), 0.5)
})

test_that("seq_label check_overlap renders", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() +
    geom_seq_label(check_overlap = TRUE)
  expect_no_error(ggplot_build(p))
})

test_that("ribbon subject color scheme works", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_color_scheme = "subject")
  expect_no_error(ggplot_build(p))
})

test_that("theme customization via + works", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() +
    ggplot2::theme(legend.position = "bottom")
  expect_no_error(ggplot_build(p))
})



test_that("legend keys are transparent regardless of panel.background", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis() +
    theme(panel.background = element_rect(fill = "grey95"))
  b <- ggplot_build(p)
  # legend keys have no fill of their own, so they blend into the background
  expect_true(is.na(b$plot$theme$legend.key$fill))
  # the default theme has no grid lines
  expect_true(is.null(b$plot$theme$panel.grid) ||
                inherits(b$plot$theme$panel.grid, "element_blank"))
})

test_that("Identity colourbar stays visible with a horizontal bottom legend", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  b <- ggplot_build(p)
  # find the colourbar guide and check it uses a fixed (non-null) key height
  heights <- character(0)
  for (gd in b$plot$guides$guides) {
    if (identical(class(gd)[1], "GuideColourbar")) {
      heights <- c(heights, as.character(gd$params$theme$legend.key.height))
    }
  }
  expect_gt(length(heights), 0)
  expect_false(any(grepl("null", heights)))
})





test_that("gene parameters accept all flexible input formats", {
  data(seq_data_example)
  seqs <- seq_data_example$seq_id
  # 1) single value / 2) strand vector / list forms all produce per-seq lists
  f <- ggchord:::process_gene_param
  expect_identical(f(20, seqs, "p", 0)[[1]], c("+" = 20, "-" = 20))
  expect_identical(f(c("+" = -15, "-" = -45), seqs, "p", 0)[[1]],
                   c("+" = -15, "-" = -45))
  expect_identical(f(list("1" = c("+" = -15, "-" = -45),
                          "2" = c("+" = 30, "-" = -30),
                          "3" = c("+" = 15, "-" = -15),
                          "4" = c("+" = 0, "-" = 0)), seqs, "p", 0),
                   setNames(list(c("+" = -15, "-" = -45),
                                 c("+" = 30, "-" = -30),
                                 c("+" = 15, "-" = -15),
                                 c("+" = 0, "-" = 0)), seqs))
  expect_identical(f(list(20), seqs, "p", 0)[[1]], c("+" = 20, "-" = 20))
  expect_identical(f(list(c("+" = -15, "-" = -45)), seqs, "p", 0)[[1]],
                   c("+" = -15, "-" = -45))
  # a named list by sequence ID
  res <- f(list("MT108731.1" = c("+" = 1, "-" = 2)), seqs, "p", 0)
  expect_identical(res[["MT108731.1"]], c("+" = 1, "-" = 2))
  expect_identical(res[["MT118296.1"]], c("+" = 0, "-" = 0))
})

test_that("sequence parameters accept list formats", {
  data(seq_data_example)
  seqs <- seq_data_example$seq_id
  f <- ggchord:::process_sequence_param
  expect_identical(f(list(5), seqs, "p"), setNames(rep(5, length(seqs)), seqs))
  expect_identical(f(list("1" = 3, "4" = 1), seqs, "p", 2),
                   setNames(c(3, 2, 2, 1), seqs))
  expect_identical(f(list(3, 2, 2, 1), seqs, "p"),
                   setNames(c(3, 2, 2, 1), seqs))
})



test_that("axis_gap defaults to 0.05", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() + geom_axis()
  # the layout stores per-sequence axis gaps; all should be 0.05
  b <- ggplot_build(p)
  axis_layer <- b$data[[length(b$data) - 1]]  # axis segments
  expect_true(nrow(axis_layer) > 0)
})





test_that("gene_label_size is applied to the gene text layer", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() + geom_gene_label()
  p <- ggchord:::prepare_ggchord_plot(p)
  layout <- p$ggchord$ref$layout
  # the layout carries the (default 2.5) label size
  expect_true(all(layout$gene_labels$size == 2.5))
  # the gene text layer maps the size aesthetic
  text_layer <- NULL
  for (lyr in p$layers) {
    if (!is.null(lyr$ggchord_type) && identical(lyr$ggchord_type, "gene_text")) text_layer <- lyr
  }
  expect_false(is.null(text_layer))
  expect_true("size" %in% names(text_layer$mapping))
  # a custom size flows through
  p2 <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() + geom_gene_label(gene_label_size = 4)
  p2 <- ggchord:::prepare_ggchord_plot(p2)
  expect_true(all(p2$ggchord$ref$layout$gene_labels$size == 4))
})

test_that("geom_gene_label_repel repels gene labels with leader lines", {
  data(seq_data_example)
  data(gene_data_example)
  # dense gene data to force overlapping labels (keep positions in range)
  gd <- gene_data_example
  dup <- transform(gd, start = pmax(1, start - 200),
                   end = pmax(200, end - 200),
                   anno = paste0(anno, " (copy)"))
  gd <- rbind(gd, dup)
  p0 <- suppressWarnings(ggchord(seq_data_example, gene_data = gd) +
                           geom_seq() + geom_gene() + geom_gene_label())
  p1 <- suppressWarnings(ggchord(seq_data_example, gene_data = gd) +
                           geom_seq() + geom_gene() +
                           geom_gene_label_repel())
  p0 <- ggchord:::prepare_ggchord_plot(p0)
  p1 <- ggchord:::prepare_ggchord_plot(p1)
  # repel should move at least one label and create leader-line segments
  expect_true(is.list(p1$ggchord$ref$layout))
  moved <- any(p0$ggchord$ref$layout$gene_labels$text_x != p1$ggchord$ref$layout$gene_labels$text_x |
                 p0$ggchord$ref$layout$gene_labels$text_y != p1$ggchord$ref$layout$gene_labels$text_y)
  expect_true(moved)
  expect_gt(nrow(p1$ggchord$ref$layout$gene_label_segments), 0)
  expect_no_error(ggplot_build(p1))
})

test_that("gene label repulsion separates rotated text boxes with varied leaders", {
  # A tightly packed arc of vertical labels exposed the former width-only,
  # radial collision test: boxes could still overlap and the layout tended to
  # preserve a common leader-line offset.
  anchors <- seq(0, 0.14, length.out = 8)
  gl <- data.frame(
    text = paste("dense label", seq_len(8)),
    text_x = anchors, text_y = rep(0, 8),
    anchor_x = anchors, anchor_y = rep(0, 8),
    text_angle = rep(90, 8), size = rep(2.5, 8)
  )
  res <- ggchord:::ggchord_repel_labels(
    gl, seed = 42, min_segment_length = 0.05,
    repel_points = data.frame(x = numeric(0), y = numeric(0))
  )
  out <- res$labels
  boxes <- ggchord:::ggchord_text_boxes(
    out, units_per_inch = 0.35, box_padding = 0.25
  )
  bw <- boxes$bw
  bh <- boxes$bh
  for (i in seq_len(nrow(out) - 1)) {
    for (j in (i + 1):nrow(out)) {
      expect_false(
        abs(out$text_x[i] - out$text_x[j]) < (bw[i] + bw[j]) / 2 &&
          abs(out$text_y[i] - out$text_y[j]) < (bh[i] + bh[j]) / 2
      )
    }
  }
  seg_len <- with(res$segments, sqrt((x1 - x0)^2 + (y1 - y0)^2))
  expect_gt(length(unique(round(seg_len, 3))), 2)
})

test_that("geom_seq_label seq_labels maps unnamed vectors positionally", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() +
    geom_seq_label(seq_labels = c("S1", "S2", "S3", "S4"))
  ggplot_build(p)
  l <- get_chord_layout()
  expect_equal(l$seq_labels_df$label, c("S1", "S2", "S3", "S4"))
  # named vectors are matched by sequence ID
  p2 <- ggchord(seq_data_example) + geom_seq() +
    geom_seq_label(seq_labels = c("MT118296.1" = "B2", "MT108731.1" = "A1"))
  ggplot_build(p2)
  l2 <- get_chord_layout()
  expect_equal(unname(l2$seq_labels_df$label[l2$seq_labels_df$seq_id == "MT118296.1"]), "B2")
  expect_equal(unname(l2$seq_labels_df$label[l2$seq_labels_df$seq_id == "MT108731.1"]), "A1")
})

test_that("geom_gene_label_repel layouts are deterministic", {
  data(seq_data_example)
  data(gene_data_example)
  build_mode <- function(mode) {
    p <- ggchord(seq_data_example, gene_data = gene_data_example) +
      geom_seq() + geom_gene() +
      geom_gene_label_repel(gene_label_layout = mode)
    ggplot_build(p)
    get_chord_layout()$gene_labels[, c("text_x", "text_y", "text_angle")]
  }
  for (mode in c("aligned", "radial", "arc")) {
    expect_equal(build_mode(mode), build_mode(mode), info = mode)
  }
  expect_error(geom_gene_label_repel(gene_label_layout = "unknown"),
               "arg.*aligned")
})

test_that("geom_gene_label_repel exposes the v0.9.0 focused interface", {
  expect_equal(
    names(formals(geom_gene_label_repel)),
    c("mapping", "data", "gene_label_layout", "gene_label_size",
      "gene_label_wrap", "gene_label_side", "max_overlaps",
      "gene_label_segment_linetype", "show_legend", "...")
  )
})

test_that("removed repel arguments fail with migration guidance", {
  position_args <- c(
    "gene_label_rotation", "gene_label_radial_offset",
    "gene_label_circum_offset", "gene_label_circum_limit"
  )
  solver_args <- c(
    "box_padding", "point_padding", "min_segment_length", "force", "seed",
    "gene_label_orientation", "gene_label_segment"
  )
  for (arg in position_args) {
    call <- setNames(list(1), arg)
    expect_error(do.call(geom_gene_label_repel, call), "geom_gene_label\\(\\)")
  }
  for (arg in solver_args) {
    call <- setNames(list(1), arg)
    expect_error(do.call(geom_gene_label_repel, call), "gene_label_layout")
  }
})

test_that("gene_label_wrap wraps long annotations", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() + geom_gene_label(gene_label_wrap = 10)
  p <- ggchord:::prepare_ggchord_plot(p)
  gl <- p$ggchord$ref$layout$gene_labels
  expect_true(any(grepl("\n", gl$text)))
  expect_no_error(ggplot_build(p))
})

test_that("gene_label_layout modes have their documented geometry", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  build_mode <- function(mode) {
    p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
      geom_seq() + geom_ribbon() + geom_gene() +
      geom_gene_label_repel(gene_label_layout = mode) + geom_axis()
    ggchord:::prepare_ggchord_plot(p)$ggchord$ref$layout
  }
  aligned <- build_mode("aligned")
  radial <- build_mode("radial")
  arc <- build_mode("arc")
  expect_true(all(aligned$gene_labels$text_angle == 0))
  expect_true(all(radial$gene_labels$text_angle == 0))
  expect_true(any(abs(arc$gene_labels$text_angle) > 1e-8))
  expect_true(all(arc$gene_labels$text_angle <= 90 |
                    arc$gene_labels$text_angle >= 270))
  expect_true(all(radial$gene_labels$label_track >= 1))
  expect_true(all(arc$gene_labels$label_track >= 1))
  expect_true(any(arc$gene_labels$label_track == 1))
  first_track <- which(arc$gene_labels$label_track == 1)
  expect_false(any(arc$gene_label_segments$group %in% first_track))
  seg <- aligned$gene_label_segments
  gl <- aligned$gene_labels
  expect_gt(nrow(seg), 0)
  # Clipping can split either leg around another label, so an elbow need not
  # remain exactly two rows. Its final visible piece must still approach the
  # label on the cardinal rail: horizontally for left/right rails and
  # vertically for top/bottom rails.
  final_rows <- vapply(seq_len(nrow(gl)), function(i) {
    candidates <- which(
      seg$group == i &
        abs(seg$x1 - gl$text_x[i]) < 1e-10 &
        abs(seg$y1 - gl$text_y[i]) < 1e-10
    )
    if (length(candidates)) candidates[length(candidates)] else NA_integer_
  }, integer(1))
  final_rows <- final_rows[!is.na(final_rows)]
  expect_gt(length(final_rows), 0)
  final <- seg[final_rows, , drop = FALSE]
  dx_stub <- final$x1 - final$x0
  dy_stub <- final$y1 - final$y0
  expect_true(all(abs(dx_stub) < 1e-10 | abs(dy_stub) < 1e-10))
  expect_true(all(sqrt(dx_stub^2 + dy_stub^2) < 0.2))
})

test_that("finite max_overlaps hides only unresolved conflicts", {
  labels <- data.frame(
    text = c("same", "same"), text_x = c(0, 0), text_y = c(0, 0),
    text_angle = c(0, 0), size = c(2.5, 2.5),
    hjust = c(0.5, 0.5), vjust = c(0.5, 0.5)
  )
  kept <- ggchord:::ggchord_hide_conflicted_labels(
    labels, max_overlaps = Inf, units_per_inch = 0.5
  )
  hidden <- ggchord:::ggchord_hide_conflicted_labels(
    labels, max_overlaps = 0, units_per_inch = 0.5
  )
  expect_equal(kept$text, labels$text)
  expect_true(all(is.na(hidden$text)))
})

test_that("horizontal repelled labels sit on the far side of the leader line", {
  data(seq_data_example)
  data(gene_data_example)
  # elbow mode
  p <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() +
    geom_gene_label_repel(gene_label_wrap = 15)
  ggplot_build(p)
  l <- get_chord_layout()
  gl <- l$gene_labels
  seg <- l$gene_label_segments
  if (nrow(seg) > 0 && nrow(gl) > 0) {
    text_boxes <- ggchord:::ggchord_text_boxes(
      gl,
      units_per_inch = max(1, diff(range(c(seg$x0, seg$x1)))) / 6
    )
    for (i in seq_len(nrow(gl))) {
      g <- which(seg$group == i)
      if (length(g) == 0) next
      stub <- seg[g[length(g)], ]
      # the elbow stub must never cross the interior of the text box
      if (abs(stub$x1 - stub$x0) >= abs(stub$y1 - stub$y0)) {
        expect_false(
          max(stub$x0, stub$x1) > text_boxes$xmin[i] + 1e-3 &&
            min(stub$x0, stub$x1) < text_boxes$xmax[i] - 1e-3,
          info = paste("label", i, "(", gl$text[i], ") stub crosses its text")
        )
      } else {
        expect_false(
          max(stub$y0, stub$y1) > text_boxes$ymin[i] + 1e-3 &&
            min(stub$y0, stub$y1) < text_boxes$ymax[i] - 1e-3,
          info = paste("label", i, "(", gl$text[i], ") stub crosses its text")
        )
      }
    }
    # the elbow must not double back: the oblique segment and the horizontal
    # stub may both be slanted/vertical, but never point in opposite
    # horizontal directions
    for (g in unique(seg$group)) {
      rows <- seg[seg$group == g, ]
      if (nrow(rows) < 2) next
      horizontal <- abs(rows$x1[2] - rows$x0[2]) >=
        abs(rows$y1[2] - rows$y0[2])
      long_component <- if (horizontal) {
        sign(rows$x1[1] - rows$x0[1])
      } else {
        sign(rows$y1[1] - rows$y0[1])
      }
      stub_component <- if (horizontal) {
        sign(rows$x1[2] - rows$x0[2])
      } else {
        sign(rows$y1[2] - rows$y0[2])
      }
      if (long_component != 0 && stub_component != 0 &&
          long_component != stub_component) {
        fail(paste("elbow for label", g, "doubles back on itself"))
      }
    }
  }
  # Default aligned mode justifies text away from the gene anchor.
  p2 <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() +
    geom_gene_label_repel()
  ggplot_build(p2)
  l2 <- get_chord_layout()
  seg2 <- l2$gene_label_segments
  gl2 <- l2$gene_labels
  if (nrow(seg2) > 0) {
    # In elbow mode the first leg may be vertical. The final leg is the one
    # that approaches the text and therefore determines its justification.
    final_leg <- !duplicated(seg2$group, fromLast = TRUE)
    final <- seg2[final_leg, ]
    horizontal <- abs(final$x1 - final$x0) >= abs(final$y1 - final$y0)
    if (any(horizontal)) {
      moved_right <- (final$x1[horizontal] - final$x0[horizontal]) >= 0
      expect_equal(
        gl2$hjust[final$group[horizontal]], ifelse(moved_right, 0, 1)
      )
    }
    if (any(!horizontal)) {
      moved_up <- (final$y1[!horizontal] - final$y0[!horizontal]) >= 0
      expect_equal(
        gl2$vjust[final$group[!horizontal]], ifelse(moved_up, 0, 1)
      )
    }
  }
})

test_that("gene_label_side moves labels to the requested arc side", {
  data(seq_data_example)
  data(gene_data_example)
  for (side in c("outside", "inside")) {
    p <- ggchord(seq_data_example, gene_data = gene_data_example) +
      geom_seq() + geom_gene() +
      geom_gene_label_repel(
        gene_label_side = side
      )
    ggplot_build(p)
    l <- get_chord_layout()
    gl <- l$gene_labels
    seg <- l$gene_label_segments
    expect_true("side_flipped" %in% names(gl))
    expect_true(any(gl$side_flipped))
    # Every label sits on the requested side of its actual sequence curve.
    # A global radius comparison is not valid for straight or strongly curved
    # paths, nor when individual sequence radii differ.
    label_delta <- ggchord:::ggchord_label_curve_frame(
      gl, l$seq_arcs
    )$signed_distance
    if (side == "outside") {
      expect_true(all(label_delta > -0.02))
    } else {
      expect_true(all(label_delta < 0.02))
    }
    # leader lines: flipped labels get dashed, the others stay solid
    expect_true(nrow(seg) > 0)
    flipped <- gl$side_flipped[match(seg$group, seq_len(nrow(gl)))]
    expect_equal(unique(seg$linetype[flipped]), "dashed")
    expect_equal(unique(seg$linetype[!flipped]), "solid")
  }
})

test_that("gene_label_segment_linetype overrides the auto dash behaviour", {
  data(seq_data_example)
  data(gene_data_example)
  build_lt <- function(...) {
    p <- ggchord(seq_data_example, gene_data = gene_data_example) +
      geom_seq() + geom_gene() +
      geom_gene_label_repel(gene_label_side = "outside", ...)
    ggplot_build(p)
    get_chord_layout()$gene_label_segments$linetype
  }
  # a forced linetype applies to every leader line, flipped or not
  expect_equal(unique(build_lt(gene_label_segment_linetype = "dotted")), "dotted")
  expect_equal(unique(build_lt(gene_label_segment_linetype = "solid")), "solid")
  expect_equal(unique(build_lt(gene_label_segment_linetype = 2)), 2)
  # invalid linetypes are rejected at layer creation
  expect_error(
    geom_gene_label_repel(gene_label_segment_linetype = "zigzag"),
    "linetype"
  )
  # default "auto" stays solid when nothing was moved to the other side
  p0 <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() +
    geom_gene_label_repel(gene_label_side = "auto")
  ggplot_build(p0)
  expect_equal(unique(get_chord_layout()$gene_label_segments$linetype), "solid")
})

test_that("elbow leader lines adapt their segment lengths per label", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() +
    geom_gene_label_repel(
      gene_label_side = "outside"
    )
  ggplot_build(p)
  seg <- get_chord_layout()$gene_label_segments
  expect_true(nrow(seg) >= 4)
  groups <- unique(seg$group)
  stub_len <- vapply(groups, function(g) {
    rows <- seg[seg$group == g, ]
    sqrt((rows$x1[2] - rows$x0[2])^2 +
           (rows$y1[2] - rows$y0[2])^2)
  }, numeric(1))
  # span between the gene anchor and the label along the final approach axis
  span <- vapply(groups, function(g) {
    rows <- seg[seg$group == g, ]
    horizontal <- abs(rows$x1[2] - rows$x0[2]) >=
      abs(rows$y1[2] - rows$y0[2])
    if (horizontal) {
      abs(rows$x1[2] - rows$x0[1])
    } else {
      abs(rows$y1[2] - rows$y0[1])
    }
  }, numeric(1))
  # stubs scale with each label's position instead of being one fixed length
  expect_gt(length(unique(round(stub_len, 4))), 1)
  # Stubs stay >= 0.02 unless the label sits almost vertically above/below
  # its gene, or the stub is collapsed to keep two leaders from crossing.
  expect_true(all(
    stub_len < 1e-7 | stub_len >= pmin(0.02, span) - 1e-6
  ))
  # the bend never lands beyond the gene anchor (no doubled-back elbows)
  for (g in groups) {
    rows <- seg[seg$group == g, ]
    horizontal <- abs(rows$x1[2] - rows$x0[2]) >=
      abs(rows$y1[2] - rows$y0[2])
    long_component <- if (horizontal) {
      sign(rows$x1[1] - rows$x0[1])
    } else {
      sign(rows$y1[1] - rows$y0[1])
    }
    stub_component <- if (horizontal) {
      sign(rows$x1[2] - rows$x0[2])
    } else {
      sign(rows$y1[2] - rows$y0[2])
    }
    expect_true(
      long_component == 0 || stub_component == 0 ||
        long_component == stub_component,
      info = paste("elbow for label group", g, "doubles back")
    )
  }
})

test_that("axis labels are aligned outward and can be hidden on overlap", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() + geom_axis()
  p <- ggchord:::prepare_ggchord_plot(p)
  layout <- p$ggchord$ref$layout
  at <- layout$axis_ticks
  # labels carry outward hjust/vjust (rotated labels may be centered on one axis)
  expect_true(all(at$label_hjust[!is.na(at$label)] %in% c(0, 0.5, 1)))
  expect_true(all(at$label_vjust[!is.na(at$label)] %in% c(0, 0.5, 1)))
  # hide-overlaps option renders without error
  p2 <- ggchord(seq_data_example) + geom_seq() +
    geom_axis(axis_label_hide_overlaps = TRUE)
  expect_no_error(ggplot_build(p2))
})

test_that("axis_label_orientation rotates axis labels", {
  data(seq_data_example)
  built_angles <- function(p) {
    b <- ggplot_build(p)
    unlist(lapply(b$data, function(d) {
      if (!is.null(d) && "label" %in% names(d) && any(!is.na(d$label))) d$angle else NULL
    }))
  }

  # default: parallel to the axis (text direction follows the arc)
  p0 <- ggchord(seq_data_example) + geom_seq() + geom_axis()
  angles0 <- built_angles(p0)
  expect_true(all(angles0 != 0))  # parallel labels are not horizontal

  # numeric vector: one angle per sequence
  p1 <- ggchord(seq_data_example) + geom_seq() +
    geom_axis(axis_label_orientation = c(0, 45, 80, 130))
  angles1 <- built_angles(p1)
  expect_setequal(angles1, c(0, 45, 80, 130))

  # named vector mixes numeric angles and "horizontal"
  p2 <- ggchord(seq_data_example) + geom_seq() +
    geom_axis(axis_label_orientation = c("MT108731.1" = 90, "MT118296.1" = "horizontal"))
  angles2 <- built_angles(p2)
  expect_setequal(angles2, c(0, 90))

  # "horizontal" keeps the text horizontal (angle 0)
  p3 <- ggchord(seq_data_example) + geom_seq() +
    geom_axis(axis_label_orientation = "horizontal")
  expect_equal(unique(built_angles(p3)), 0)

  # "parallel" and "perpendicular" align with / normal to the axis direction
  l4 <- get_chord_layout()
  align_dev <- function(mode, expected) {
    p4 <- ggchord(seq_data_example) + geom_seq() +
      geom_axis(axis_label_orientation = mode)
    ggplot_build(p4)
    tl <- get_chord_layout()$axis_ticks
    tl <- tl[!is.na(tl$label), ]
    al <- get_chord_layout()$axis_lines
    d2 <- outer(al$x, tl$label_x, function(a, b) (a - b)^2) +
      outer(al$y, tl$label_y, function(a, b) (a - b)^2)
    k <- apply(d2, 2, which.min)
    k2 <- ifelse(k < nrow(al), k + 1, k - 1)
    seg_ang <- atan2(al$y[k2] - al$y[k], al$x[k2] - al$x[k]) * 180 / pi
    dev <- abs(tl$label_angle %% 180 - seg_ang %% 180)
    dev <- ifelse(dev > 90, 180 - dev, dev)
    max(abs(dev - expected))
  }
  expect_lt(align_dev("parallel", 0), 10)
  expect_lt(align_dev("perpendicular", 90), 10)

  # rotated labels carry outward justification in the layout
  rot <- l4$axis_ticks[!is.na(l4$axis_ticks$label), ]
  expect_true(all(rot$label_hjust %in% c(0, 0.5, 1)))
  expect_true(all(rot$label_vjust %in% c(0, 0.5, 1)))
})

test_that("legend_key_width and legend_key_height control the Identity colourbar key", {
  data(seq_data_example)
  data(ribbon_data_example)
  p0 <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
  p1 <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() +
    geom_ribbon(legend_key_height = 5, legend_key_width = 1)
  b0 <- ggplot_build(p0)
  b1 <- ggplot_build(p1)
  find_key <- function(b) {
    for (s in b$plot$scales$non_position_scales()$scales) {
      if (identical(class(s$guide)[1], "GuideColourbar")) {
        return(list(
          height = as.character(s$guide$params$theme$legend.key.height),
          width = as.character(s$guide$params$theme$legend.key.width)
        ))
      }
    }
    list(height = "none", width = "none")
  }
  expect_equal(find_key(b0)$height, "1null")
  expect_equal(find_key(b1)$height, "5cm")
  expect_equal(find_key(b1)$width, "1cm")
})

test_that("documented data and parameter values are validated", {
  expect_error(
    ggchord(data.frame(seq_id = c("a", "a"), length = c(1, 2))),
    "unique"
  )

  data(seq_data_example)
  expect_error({
    p <- ggchord(seq_data_example) + geom_seq(seq_orientation = 0)
    ggplot_build(p)
  }, "1 or -1")
})
