seqs <- data.frame(accver = c("A", "B"), length = 1000)
links <- data.frame(qaccver = "A", saccver = "B", qstart = 100,
  qend = 200, sstart = 300, send = 400, qstrand = "+", sstrand = "+")
genes <- data.frame(accver = c("A", "B"), start = c(100, 400),
  end = c(200, 500), strand = c("+", "-"), label = c("alpha", "beta"))

test_that("core layers build and keep manual gene colours", {
  colours <- c("+" = "#56887C", "-" = "#B58169")
  p <- ggchord() + geom_seq(data = seqs) + geom_link_ribbon(data = links) +
    geom_gene(ggplot2::aes(fill = strand), data = genes) +
    ggplot2::scale_fill_manual(values = colours)
  built <- ggplot2::ggplot_build(p)
  fills <- unlist(lapply(built$data, function(layer) layer$fill))
  expect_true(all(colours %in% fills))
  expect_true(ggplot2::is_waiver(p$data))
})

test_that("fixed and repelled labels retain representative annotations", {
  p <- ggchord(seqs) + geom_gene(data = genes, id = "genes") + coord_collinear()
  for (layer in list(geom_gene_label(anchor_layer = "genes"),
                     geom_gene_label(anchor_layer = "genes", label_layout = layout_label_repel()))) {
    built <- ggplot2::ggplot_build(p + layer)
    labels <- unlist(lapply(built$data, function(layer) as.character(layer$label)))
    expect_true(all(genes$label %in% labels))
  }
})

test_that("feature, region and point-link layers retain their source tables", {
  points <- data.frame(qaccver = "A", saccver = "B", qpos = 100, spos = 400)
  p <- ggchord(seqs) + geom_feature(data = genes, id = "genes") +
    geom_seq_region(data = genes, id = "region") +
    geom_link_line(data = points, id = "points")
  out <- export_ggchord_layout(p, stages = "source")
  expect_setequal(out$stage_status$layer_id, c("genes", "region", "points"))
  expect_true(all(out$stage_status$status == "available"))
})

test_that("explicit ribbon styles follow all sequence relationship projections", {
  for (coord in list(coord_chord(), coord_collinear(), coord_radial_synteny(), coord_sankey())) {
    p <- ggchord(seqs) + geom_link_ribbon(data = links, colour = "black",
      linetype = "dotted") + coord
    built <- ggplot2::ggplot_build(p)
    expect_true("dotted" %in% built$data[[1]]$linetype)
  }
})

test_that("empty primer results compose with circular layers", {
  hits <- find_primer_bindings(c(A = "AAAAAAAAAAAAAAAAAAAA"),
    primers = c(test = "ACGTTGCA"), min_annealed_bases = 8)
  p <- ggchord(data.frame(accver = "A", length = 20)) + geom_seq() +
    geom_primer(data = hits) + geom_primer_label(data = hits, label_layout = layout_label_repel()) +
    coord_circular()
  out <- export_ggchord_layout(p, stages = "source")
  expect_true(all(out$stage_status$status == "available"))
  expect_equal(sum(vapply(out$source, function(layer) nrow(layer$data) == 0L, logical(1))), 2L)
})

test_that("Unicode measurement preserves content and the active PDF device", {
  grDevices::pdf(NULL)
  caller <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(caller), add = TRUE)
  annotation <- data.frame(text = "lacZ\u03b1", text_x = 0, text_y = 0, size = 2.5)
  boxes <- ggchord:::ggchord_text_boxes(annotation)
  expect_true(all(is.finite(boxes$w) & boxes$w > 0))
  expect_identical(annotation$text, "lacZ\u03b1")
  expect_identical(grDevices::dev.cur(), caller)
})
