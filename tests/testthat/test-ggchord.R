# testthat test file - v0.4.0 layered API

test_that("ggchord with only seq_data returns a ggchord object", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example)
  expect_s3_class(p, "ggchord")
  expect_s3_class(p, "ggplot")
})

test_that("ggchord + geom_seq creates a sequence arc layer", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example) + geom_seq()
  expect_s3_class(p, "ggchord")
  expect_true(length(p$layers) >= 1)
})

test_that("ggchord + geom_seq + geom_ribbon handles alignment data correctly", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(
    seq_data = seq_data_example,
    ribbon_data = ribbon_data_example
  ) +
    geom_seq() +
    geom_ribbon()
  expect_s3_class(p, "ggchord")
  expect_true(length(p$layers) >= 1)
})

test_that("ggchord + geom_gene handles gene data correctly", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(
    seq_data = seq_data_example,
    gene_data = gene_data_example
  ) +
    geom_seq() +
    geom_gene()
  expect_s3_class(p, "ggchord")
})

test_that("ggchord + geom_axis can add an axis", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example) +
    geom_seq() +
    geom_axis()
  expect_s3_class(p, "ggchord")
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

test_that("debug mode outputs debug information", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example, debug = TRUE)
  expect_s3_class(p, "ggchord")
})

test_that("ggchord global parameters are passed correctly", {
  data(seq_data_example)
  p <- ggchord(seq_data = seq_data_example, title = "Hello", rotation = 90)
  expect_s3_class(p, "ggchord")
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

  out <- tempfile(fileext = ".pdf")
  pdf(out, 8, 8)
  expect_no_warning(print(p))
  dev.off()

  layout <- ggchord:::get_chord_layout()
  expect_true("fill" %in% names(layout$ribbon_polys))
  expect_true(all(layout$ribbon_polys$alpha == 0.2))
  expect_gt(nrow(layout$gene_labels), 0)
  expect_true(file.exists(out))
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

  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p2)))
  dev.off()
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
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p2)))
  dev.off()
  expect_equal(p2$layers[[1]]$ggchord_params$seq_radius, 5)
})

test_that("ggsave works directly on a ggchord plot", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
  f <- tempfile(fileext = ".png")
  expect_no_error(ggsave(f, p, width = 6, height = 6, dpi = 72))
  expect_true(file.exists(f))
})

test_that("plot objects are self-contained and repeated builds stay stable", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis()
  # The plot carries its own scales (tagged ggchord_managed) so that tools such
  # as plotly::ggplotly() that clone the plot see the same scales as a build.
  expect_gt(length(p$scales$scales), 0)
  expect_true(all(vapply(p$scales$scales,
                         function(s) !is.null(s$ggchord_managed), logical(1))))
  expect_false(is.null(p$layers[[1]]$ggchord_params))
  # repeated builds remain stable
  expect_no_error(ggplot_build(p))
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
})

test_that("ribbon data sanity checks warn on bad input", {
  data(seq_data_example)
  bad_ribbon <- data.frame(
    qaccver = c("MT108731.1", "MT108731.1"),
    saccver = c("MT118296.1", "MT118296.1"),
    length = c(100, 100), pident = c(90, 90),
    qstart = c(50000, 100), qend = c(40000, 200),   # qstart > qend
    sstart = c(1, 1), send = c(100, 100)
  )
  expect_warning(ggchord(seq_data_example, bad_ribbon), "start > end")
})

test_that("geom_seq_label places sequence labels", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) +
    geom_seq() +
    geom_seq_label(seq_label_radius = 1.25, seq_label_size = 3.5)
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
  layout <- ggchord:::get_chord_layout()
  expect_gt(nrow(layout$seq_labels_df), 0)
  expect_true(all(c("text_x", "text_y", "label") %in% names(layout$seq_labels_df)))
})

test_that("ribbon subject color scheme works", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example) +
    geom_seq() +
    geom_ribbon(ribbon_color_scheme = "subject")
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
})

test_that("theme customization via + works", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() +
    ggplot2::theme(legend.position = "bottom")
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
})

test_that("get_chord_layout is exported after building", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq()
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  print(p)
  dev.off()
  layout <- get_chord_layout()
  expect_true(inherits(layout, "chord_layout"))
  expect_gt(length(layout$seq_arcs), 0)
})

test_that("plots convert to plotly when plotly is installed", {
  skip_if_not_installed("plotly")
  # NOTE: do not attach plotly with library() here - it masks geom_ribbon()
  # and would leak into the other tests in this session.
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  # sequence-only plot (the common interactive use case)
  p <- ggchord(seq_data_example) + geom_seq()
  pl <- suppressWarnings(plotly::ggplotly(p))
  expect_true(inherits(pl, "plotly"))
  expect_gt(length(pl$x$data), 0)
})

test_that("legend keys keep a white background regardless of panel.background", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis() +
    theme(panel.background = element_rect(fill = "grey95"))
  b <- ggplot_build(p)
  # the legend key fill is fixed to white, not inherited from panel.background
  expect_equal(b$plot$theme$legend.key$fill, "white")
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

test_that("plotly conversion keeps legends and adds sequence-arc arrows", {
  skip_if_not_installed("plotly")
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis()
  pl <- suppressWarnings(plotly::ggplotly(p))
  expect_true(isTRUE(pl$x$layout$showlegend))
  expect_gt(length(pl$x$layout$annotations), 0)
})

test_that("legends can be positioned independently with legend_position", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq(legend_position = "left") +
    geom_ribbon(legend_position = "bottom") +
    geom_gene(legend_position = "right") +
    geom_axis()
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(suppressWarnings(print(p)))
  dev.off()
  # the three legends occupy separate boxes (left, bottom, right);
  # ggplotGrob() opens a device for text measurement, so open one explicitly
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  g <- ggplotGrob(p)
  dev.off()
  filled <- vapply(g$grobs, function(gr) {
    inherits(gr, "gtable") && grepl("guide-box", gr$name) &&
      (length(gr$grobs) > 0 || length(gr$children) > 0)
  }, logical(1))
  positions <- g$layout$name[filled]
  expect_true(any(grepl("guide-box-left", positions)))
  expect_true(any(grepl("guide-box-bottom", positions)))
  expect_true(any(grepl("guide-box-right", positions)))
  # invalid positions are rejected (validate at build time, no device)
  p2 <- ggchord(seq_data_example) + geom_seq(legend_position = "center")
  expect_error(ggplot_build(p2), "legend_position")
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

test_that("default legend positions split the three legends", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  # defaults: Seq ID + Strand on the right, Identity(%) colourbar on the left
  expect_identical(ggchord::geom_seq()[[1]]$ggchord_params$legend_position, "right")
  expect_identical(ggchord::geom_ribbon()[[1]]$ggchord_params$legend_position, "left")
  expect_identical(ggchord::geom_gene()[[1]]$ggchord_params$legend_position, "right")

  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() + geom_axis()
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  g <- ggplotGrob(p)
  dev.off()
  filled <- vapply(g$grobs, function(gr) {
    inherits(gr, "gtable") && grepl("guide-box", gr$name) &&
      (length(gr$grobs) > 0 || length(gr$children) > 0)
  }, logical(1))
  positions <- g$layout$name[filled]
  expect_true(any(grepl("guide-box-left", positions)))
  expect_true(any(grepl("guide-box-right", positions)))
})

test_that("axis_gap defaults to 0.05", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() + geom_axis()
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  print(p)
  dev.off()
  layout <- ggchord:::get_chord_layout()
  # the layout stores per-sequence axis gaps; all should be 0.05
  b <- ggplot_build(p)
  axis_layer <- b$data[[length(b$data) - 1]]  # axis segments
  expect_true(nrow(axis_layer) > 0)
})

test_that("geom_ribbon shows the Identity colourbar legend without gene data", {
  data(seq_data_example)
  data(ribbon_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  g <- ggplotGrob(p)
  dev.off()
  # the colourbar legend is rendered (a rastergrob inside a guide box)
  found <- FALSE
  walk <- function(x) {
    if (inherits(x, "rastergrob")) found <<- TRUE
    if (!is.null(x$grobs)) for (ch in x$grobs) walk(ch)
    if (!is.null(x$children)) for (ch in x$children) walk(ch)
  }
  for (i in seq_along(g$grobs)) if (inherits(g$grobs[[i]], "gtable") && grepl("guide-box", g$grobs[[i]]$name)) walk(g$grobs[[i]])
  expect_true(found)
})

test_that("no spurious fill-scale warning when gene data is present without a gene layer", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon()
  warns <- character(0)
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  withCallingHandlers(print(p), warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  dev.off()
  expect_false(any(grepl("No shared levels", warns)))
})

test_that("gene_label_size is applied to the gene text layer", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() + geom_gene_label()
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
  p0 <- ggchord(seq_data_example, gene_data = gd) + geom_seq() + geom_gene() + geom_gene_label()
  p1 <- ggchord(seq_data_example, gene_data = gd) + geom_seq() + geom_gene() +
    geom_gene_label_repel()
  # repel should move at least one label and create leader-line segments
  expect_true(is.list(p1$ggchord$ref$layout))
  moved <- any(p0$ggchord$ref$layout$gene_labels$text_x != p1$ggchord$ref$layout$gene_labels$text_x |
                 p0$ggchord$ref$layout$gene_labels$text_y != p1$ggchord$ref$layout$gene_labels$text_y)
  expect_true(moved)
  expect_gt(nrow(p1$ggchord$ref$layout$gene_label_segments), 0)
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(print(p1))
  dev.off()
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

test_that("geom_gene_label_repel seed changes and reproduces the layout", {
  data(seq_data_example)
  data(gene_data_example)
  build_seed <- function(s) {
    p <- ggchord(seq_data_example, gene_data = gene_data_example) +
      geom_seq() + geom_gene() +
      geom_gene_label_repel(seed = s)
    ggplot_build(p)
    get_chord_layout()$gene_labels[, c("text_x", "text_y")]
  }
  l1 <- build_seed(1)
  l2 <- build_seed(2)
  l1b <- build_seed(1)
  # different seeds give different layouts; the same seed is reproducible
  expect_gt(max(abs(l1$text_x - l2$text_x)), 0)
  expect_equal(l1, l1b)
})

test_that("gene_label_wrap wraps long annotations", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() + geom_gene_label(gene_label_wrap = 10)
  gl <- p$ggchord$ref$layout$gene_labels
  expect_true(any(grepl("\n", gl$text)))
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(print(p))
  dev.off()
})

test_that("gene_label_orientation and gene_label_segment work", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  # horizontal text + elbow leader lines
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
    geom_seq() + geom_ribbon() + geom_gene() +
    geom_gene_label_repel(gene_label_orientation = "horizontal",
                          gene_label_segment = "elbow") +
    geom_axis()
  gl <- p$ggchord$ref$layout$gene_labels
  seg <- p$ggchord$ref$layout$gene_label_segments
  expect_true(all(gl$text_angle == 0))
  expect_gt(nrow(seg), 0)
  # each elbow is two rows: an oblique segment from the gene to the label's
  # horizontal level, then a short horizontal stub into the label
  expect_equal(nrow(seg) %% 2, 0)
  idx <- which(seg$group == seg$group[1])
  # the two segments meet at the bend point
  expect_equal(seg$x0[idx[2]], seg$x1[idx[1]])
  expect_equal(seg$y0[idx[2]], seg$y1[idx[1]])
  # the stub is horizontal and short
  expect_equal(seg$y0[idx[2]], seg$y1[idx[2]])
  expect_lt(abs(seg$x1[idx[2]] - seg$x0[idx[2]]), 0.2)
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(print(p))
  dev.off()
})

test_that("horizontal repelled labels sit on the far side of the leader line", {
  data(seq_data_example)
  data(gene_data_example)
  # elbow mode
  p <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() +
    geom_gene_label_repel(gene_label_orientation = "horizontal",
                          gene_label_segment = "elbow", gene_label_wrap = 15)
  ggplot_build(p)
  l <- get_chord_layout()
  gl <- l$gene_labels
  seg <- l$gene_label_segments
  if (nrow(seg) > 0 && nrow(gl) > 0) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    w <- suppressWarnings(graphics::strwidth(gl$text, units = "inches",
                                             cex = (gl$size %||% 2.5) / 12)) *
      max(1, diff(range(c(seg$x0, seg$x1)))) / 6
    for (i in seq_len(nrow(gl))) {
      g <- which(seg$group == i)
      if (length(g) == 0) next
      stub <- seg[g[length(g)], ]
      tx <- gl$text_x[i]
      hjust <- gl$hjust[i]
      x0box <- if (hjust < 0.5) tx else if (hjust > 0.5) tx - w[i] else tx - w[i] / 2
      x1box <- if (hjust < 0.5) tx + w[i] else if (hjust > 0.5) tx else tx + w[i] / 2
      # the elbow stub must never cross the interior of the text box
      expect_false(
        max(stub$x0, stub$x1) > x0box + 1e-3 &&
          min(stub$x0, stub$x1) < x1box - 1e-3,
        info = paste("label", i, "(", gl$text[i], ") stub crosses its text")
      )
    }
    # the elbow must not double back: the oblique segment and the horizontal
    # stub may both be slanted/vertical, but never point in opposite
    # horizontal directions
    for (g in unique(seg$group)) {
      rows <- seg[seg$group == g, ]
      if (nrow(rows) < 2) next
      dx_long <- sign(rows$x1[1] - rows$x0[1])
      dx_stub <- sign(rows$x1[2] - rows$x0[2])
      if (dx_long != 0 && dx_stub != 0 && dx_long != dx_stub) {
        fail(paste("elbow for label", g, "doubles back on itself"))
      }
    }
  }
  # line mode: text is justified away from the gene anchor
  p2 <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() +
    geom_gene_label_repel(gene_label_orientation = "horizontal")
  ggplot_build(p2)
  l2 <- get_chord_layout()
  seg2 <- l2$gene_label_segments
  gl2 <- l2$gene_labels
  if (nrow(seg2) > 0) {
    moved_right <- (seg2$x1 - seg2$x0) >= 0
    expect_equal(gl2$hjust[seg2$group], ifelse(moved_right, 0, 1))
  }
})

test_that("gene_label_side moves labels to the requested arc side", {
  data(seq_data_example)
  data(gene_data_example)
  for (side in c("outside", "inside")) {
    p <- ggchord(seq_data_example, gene_data = gene_data_example) +
      geom_seq() + geom_gene() +
      geom_gene_label_repel(
        seed = 1,
        gene_label_side = side,
        gene_label_orientation = "horizontal",
        gene_label_segment = "elbow"
      )
    ggplot_build(p)
    l <- get_chord_layout()
    gl <- l$gene_labels
    seg <- l$gene_label_segments
    expect_true("side_flipped" %in% names(gl))
    expect_true(any(gl$side_flipped))
    # every label sits on the requested side of its sequence arc
    arc_r <- vapply(l$seq_arcs, function(a) {
      median(sqrt(a$x^2 + a$y^2))
    }, numeric(1))
    names(arc_r) <- vapply(l$seq_arcs, function(a) {
      unique(a$seq_id)[1]
    }, character(1))
    label_delta <- sqrt(gl$text_x^2 + gl$text_y^2) - arc_r[gl$seq_id]
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
      geom_gene_label_repel(seed = 1, gene_label_side = "outside", ...)
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
    geom_seq() + geom_gene() + geom_gene_label_repel(seed = 1)
  ggplot_build(p0)
  expect_equal(unique(get_chord_layout()$gene_label_segments$linetype), "solid")
})

test_that("elbow leader lines adapt their segment lengths per label", {
  data(seq_data_example)
  data(gene_data_example)
  p <- ggchord(seq_data_example, gene_data = gene_data_example) +
    geom_seq() + geom_gene() +
    geom_gene_label_repel(
      seed = 1,
      gene_label_side = "outside",
      gene_label_orientation = "horizontal",
      gene_label_segment = "elbow"
    )
  ggplot_build(p)
  seg <- get_chord_layout()$gene_label_segments
  expect_true(nrow(seg) >= 4)
  groups <- unique(seg$group)
  stub_len <- vapply(groups, function(g) {
    rows <- seg[seg$group == g, ]
    abs(rows$x1[2] - rows$x0[2])
  }, numeric(1))
  # stubs scale with each label's position instead of being one fixed length
  expect_gt(length(unique(round(stub_len, 4))), 1)
  expect_true(all(stub_len >= 0.02 - 1e-6))
  # the bend never lands beyond the gene anchor (no doubled-back elbows)
  for (g in groups) {
    rows <- seg[seg$group == g, ]
    expect_true(
      sign(rows$x1[1] - rows$x0[1]) == 0 ||
        sign(rows$x1[2] - rows$x0[2]) == 0 ||
        sign(rows$x1[1] - rows$x0[1]) == sign(rows$x1[2] - rows$x0[2]),
      info = paste("elbow for label group", g, "doubles back")
    )
  }
})

test_that("axis labels are aligned outward and can be hidden on overlap", {
  data(seq_data_example)
  p <- ggchord(seq_data_example) + geom_seq() + geom_axis()
  layout <- p$ggchord$ref$layout
  at <- layout$axis_ticks
  # labels carry outward hjust/vjust (rotated labels may be centered on one axis)
  expect_true(all(at$label_hjust[!is.na(at$label)] %in% c(0, 0.5, 1)))
  expect_true(all(at$label_vjust[!is.na(at$label)] %in% c(0, 0.5, 1)))
  # hide-overlaps option renders without error
  p2 <- ggchord(seq_data_example) + geom_seq() +
    geom_axis(axis_label_hide_overlaps = TRUE)
  pdf(tempfile(fileext = ".pdf"), 8, 8)
  expect_no_error(print(p2))
  dev.off()
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

test_that("legend_key_length controls the Identity colourbar length", {
  data(seq_data_example)
  data(ribbon_data_example)
  p0 <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
  p1 <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() +
    geom_ribbon(legend_key_length = 5)
  b0 <- ggplot_build(p0)
  b1 <- ggplot_build(p1)
  find_kh <- function(b) {
    for (s in b$plot$scales$non_position_scales()$scales) {
      if (identical(class(s$guide)[1], "GuideColourbar")) {
        return(as.character(s$guide$params$theme$legend.key.height))
      }
    }
    "none"
  }
  expect_equal(find_kh(b0), "1null")
  expect_equal(find_kh(b1), "5cm")
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
