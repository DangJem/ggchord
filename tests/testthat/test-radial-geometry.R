count_radial_crossings <- function(paths) {
  if (nrow(paths) < 2L) return(0L)
  sum(vapply(seq_len(nrow(paths) - 1L), function(i) {
    j <- seq.int(i + 1L, nrow(paths))
    j <- j[paths$group[j] != paths$group[i]]
    sum(ggchord:::ggchord_segments_cross(paths$x0[i], paths$y0[i],
      paths$x1[i], paths$y1[i], paths$x0[j], paths$y0[j], paths$x1[j], paths$y1[j]))
  }, integer(1)))
}

test_that("radial is default and removed aligned is rejected", {
  expect_identical(formals(geom_gene_label_repel)$gene_label_layout, "radial")
  expect_error(geom_gene_label_repel(gene_label_layout = "aligned"), "auto")
  seq <- data.frame(seq_id = "circle", length = 10000)
  genes <- data.frame(seq_id = "circle", start = seq(300, 9000, length.out = 12),
    end = seq(300, 9000, length.out = 12) + 100, strand = "+",
    anno = paste0("Site", seq_len(12)))
  p <- ggchord(seq, gene_data = genes, validate = "none") + geom_seq() + geom_gene()
  implicit <- get_chord_layout(p + geom_gene_label_repel())
  explicit <- get_chord_layout(p + geom_gene_label_repel(gene_label_layout = "radial"))
  expect_equal(implicit$gene_labels, explicit$gene_labels)
  expect_equal(nrow(implicit$gene_labels), nrow(genes))
  expect_true(all(implicit$gene_labels$text_angle == 0))
  expect_false(ggchord:::ggchord_label_box_conflicts(implicit$gene_labels,
    units_per_inch = implicit$text_units_per_inch, box_padding = 0))
  expect_equal(count_radial_crossings(implicit$gene_label_segments), 0L)
})

test_that("signed curvature has fixed endpoints and continuous symmetric bow", {
  curve <- function(c) ggchord:::generate_curvature_path(-0.7, 0.7, 2, c, 101)
  zero <- curve(0)
  for (c in c(0.001, 0.5, 1, 2, 12)) {
    pos <- curve(c)
    neg <- curve(-c)
    expect_equal(pos[c(1, 101), ], zero[c(1, 101), ], tolerance = 1e-12)
    expect_equal(neg[c(1, 101), ], zero[c(1, 101), ], tolerance = 1e-12)
    expect_equal(pos$x + neg$x, 2 * zero$x, tolerance = 1e-12)
    expect_equal(pos$y, neg$y, tolerance = 1e-12)
  }
  expect_lt(max(abs(as.matrix(curve(1e-8) - zero))), 1e-6)
  expect_lt(max(abs(as.matrix(curve(1 + 1e-8) - curve(1)))), 1e-6)
  expect_gt(abs(curve(12)$x[51] - curve(10)$x[51]), 0.1)
  expect_error(curve(Inf), "finite")
  flat <- ggchord:::generate_curvature_path(0, 1.8 * pi, 2, 0, 101)
  expect_equal(flat$x, seq(flat$x[1], flat$x[101], length.out = 101))
  expect_equal(flat$y, seq(flat$y[1], flat$y[101], length.out = 101))
})

test_that("negative and large bows keep feature, axis and ribbon coordinates finite", {
  seq <- data.frame(seq_id = c("A", "B"), length = 1000)
  genes <- data.frame(seq_id = c("A", "B"), start = c(200, 500),
                      end = c(350, 650), strand = c("+", "-"), anno = c("one", "two"))
  ribbons <- data.frame(qaccver = "A", saccver = "B", qstart = 200,
                        qend = 350, sstart = 500, send = 650, pident = 95, length = 151)
  for (curvature in list(c(-0.5, 0.5), c(-2, 2), c(-5, 5))) {
    p <- ggchord(seq, ribbons, genes, validate = "none") +
      geom_seq(seq_curvature = curvature) + geom_gene() + geom_ribbon()
    layout <- get_chord_layout(p)
    for (data in c(layout$seq_arcs, list(layout$gene_polys, layout$ribbon_polys))) {
      expect_true(all(is.finite(c(data$x, data$y))))
    }
    expect_true(all(is.finite(c(layout$axis_ticks$x0, layout$axis_ticks$y0))))
  }
})

test_that("normal leaders keep their side when sequence orientation reverses", {
  seq <- data.frame(seq_id = c("A", "B"), length = 1000)
  genes <- data.frame(seq_id = rep(c("A", "B"), each = 2),
    start = c(200, 650, 200, 650), end = c(280, 730, 280, 730),
    strand = c("+", "-", "+", "-"), anno = c("a", "b", "c", "d"))
  for (side in c("outside", "inside")) {
    p <- ggchord(seq, gene_data = genes, validate = "none") +
      geom_seq(seq_orientation = c(1, -1), seq_curvature = c(0.5, -0.5)) +
      geom_gene() + geom_gene_label_repel(gene_label_side = side)
    layout <- get_chord_layout(p)
    frame <- ggchord:::ggchord_label_curve_frame(layout$gene_labels, layout$seq_arcs)
    expected <- if (side == "outside") 1 else -1
    expect_true(all(expected * frame$signed_distance > 0))
    expect_equal(count_radial_crossings(layout$gene_label_segments), 0L)
  }
})

test_that("content-fitted previews keep explicit dimensions optional", {
  p <- ggchord(data.frame(seq_id = "circle", length = 1000), validate = "none") +
    geom_seq() + ggplot2::theme(legend.position = "none")
  height <- ggchord:::ggchord_preview_height(p, 5)
  expect_true(is.finite(height) && height > 1 && height < 8)
  expect_null(formals(view_ggchord)$height)
})

test_that("preview preserves output warnings, dimensions and the caller device", {
  grDevices::pdf(NULL, width = 6, height = 4)
  on.exit(grDevices::dev.off())
  caller_device <- grDevices::dev.cur()
  p <- ggplot2::ggplot(data.frame(x = c(1, NA), y = c(1, 2)),
    ggplot2::aes(x, y)) + ggplot2::geom_point()
  expect_warning(file <- view_ggchord(p, width = 5.08, height = 2.54,
    units = "cm", dpi = 100, viewer = "none"), "Removed 1 row")
  expect_identical(grDevices::dev.cur(), caller_device)
  # PNG IHDR stores the actual width and height as big-endian integers.
  connection <- file(file, "rb")
  on.exit(close(connection), add = TRUE)
  readBin(connection, "raw", n = 16)
  expect_identical(readBin(connection, "integer", n = 2, size = 4,
    endian = "big"), c(200L, 100L))
})

test_that("content fitting reports a legend that cannot fit the requested width", {
  p <- ggplot2::ggplot(data.frame(x = 1:3, group = c("a", "b", "c")),
    ggplot2::aes(x, x, colour = group)) + ggplot2::geom_point() +
    ggplot2::scale_colour_discrete(labels = rep("A long legend description", 3)) +
    ggplot2::theme(legend.position = "top")
  expect_warning(view_ggchord(p, width = 2, viewer = "none"),
    "legends exceed the requested width")
})

test_that("radial balances opposite sectors without losing collision guarantees", {
  grDevices::pdf(NULL, width = 11, height = 7)
  on.exit(grDevices::dev.off())
  p <- ggchord(seq_data_example, ribbon_data_example, gene_data_example,
    validate = "none") + geom_seq() + geom_gene() + geom_gene_label_repel()
  layout <- get_chord_layout(p)
  segments <- layout$gene_label_segments
  lengths <- tapply(sqrt((segments$x1 - segments$x0)^2 +
                         (segments$y1 - segments$y0)^2), segments$group, sum)
  means <- tapply(lengths[as.character(seq_len(nrow(layout$gene_labels)))],
                   layout$gene_labels$seq_id, mean)
  expect_lt(means[["MT108731.1"]] / means[["OQ646790.1"]], 1.6)
  expect_equal(count_radial_crossings(segments), 0L)
  expect_false(ggchord:::ggchord_label_box_conflicts(layout$gene_labels,
    units_per_inch = layout$text_units_per_inch, box_padding = 0))
  fitted <- ggchord:::ggchord_preview_layout(p, 11)
  expect_s3_class(fitted$plot, "gtable")
  expect_true(is.finite(fitted$height) && fitted$height > 0)
})
