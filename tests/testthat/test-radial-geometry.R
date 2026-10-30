count_radial_crossings <- function(paths) {
  if (nrow(paths) < 2L) return(0L)
  sum(vapply(seq_len(nrow(paths) - 1L), function(i) {
    j <- seq.int(i + 1L, nrow(paths))
    j <- j[paths$group[j] != paths$group[i]]
    sum(ggchord:::ggchord_segments_cross(paths$x0[i], paths$y0[i],
      paths$x1[i], paths$y1[i], paths$x0[j], paths$y0[j],
      paths$x1[j], paths$y1[j]))
  }, integer(1)))
}

test_that("nested measurements preserve the caller graphics device", {
  old_dir <- setwd(tempdir())
  on.exit(setwd(old_dir), add = TRUE)
  previous <- grDevices::dev.cur()
  original <- grDevices::dev.list()
  on.exit({
    for (device in setdiff(grDevices::dev.list(), original))
      grDevices::dev.off(device)
    if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
  }, add = TRUE)
  grDevices::pdf(NULL, width = 6, height = 4)
  caller <- grDevices::dev.cur()
  p <- ggchord(smoke_sequences(), gene_data = smoke_genes(),
               validate = "none") +
    geom_seq() + geom_gene() + geom_gene_label_repel()
  view_ggchord(p, viewer = "none")
  expect_identical(grDevices::dev.cur(), caller)
  expect_gt(nrow(get_chord_layout(p, build = FALSE)$gene_labels), 0L)
})


test_that("Unicode annotation measurement preserves content and caller device", {
  previous <- grDevices::dev.cur()
  annotation <- data.frame(text = "lacZ\u03b1", text_x = 0, text_y = 0, size = 2.5)
  boxes <- ggchord:::ggchord_text_boxes(annotation)
  expect_true(all(is.finite(boxes$w) & boxes$w > 0))
  expect_identical(grDevices::dev.cur(), previous)
})

test_that("radial labels retain inputs without overlapping or crossing", {
  seq <- data.frame(accver = "circle", length = 10000)
  genes <- data.frame(accver = "circle",
    start = seq(300, 9000, length.out = 6),
    end = seq(300, 9000, length.out = 6) + 100,
    strand = "+", anno = paste0("Site", seq_len(6)))
  p <- ggchord(seq, gene_data = genes, validate = "none") +
    geom_seq() + geom_gene() + geom_gene_label_repel()
  layout <- get_chord_layout(p)
  expect_equal(nrow(layout$gene_labels), nrow(genes))
  expect_false(ggchord:::ggchord_label_box_conflicts(layout$gene_labels,
    units_per_inch = layout$text_units_per_inch, box_padding = 0))
  expect_equal(count_radial_crossings(layout$gene_label_segments), 0L)
})

test_that("reversed orientation keeps leaders on their requested side", {
  seq <- data.frame(accver = c("A", "B"), length = 1000)
  genes <- data.frame(accver = rep(c("A", "B"), each = 2),
    start = c(200, 650, 200, 650), end = c(280, 730, 280, 730),
    strand = c("+", "-", "+", "-"), anno = letters[1:4])
  p <- ggchord(seq, gene_data = genes, validate = "none") +
    geom_seq(seq_orientation = c(1, -1)) + geom_gene() +
    geom_gene_label_repel(label_side = "outside")
  layout <- get_chord_layout(p)
  frame <- ggchord:::ggchord_label_curve_frame(layout$gene_labels,
                                                layout$seq_arcs)
  expect_true(all(frame$signed_distance > 0))
  expect_equal(count_radial_crossings(layout$gene_label_segments), 0L)
})
