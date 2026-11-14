# Tiny format fixtures; no external programs or optional biological databases.
seqs <- data.frame(accver = c("A", "B"), length = 1000)
links <- data.frame(qaccver = "A", saccver = "B",
  qstart = 100, qend = 200, sstart = 300, send = 400)
genes <- data.frame(accver = c("A", "B"), start = c(100, 400),
  end = c(200, 500), strand = c("+", "-"), anno = c("alpha", "beta"))

test_that("core layers build and keep manual gene colours", {
  colours <- c("+" = "#56887C", "-" = "#B58169")
  p <- ggchord(seqs, links, genes[c("accver", "start", "end", "strand")]) +
    geom_seq() + geom_link_ribbon() + geom_gene() + geom_seq_label() +
    scale_gene_fill_manual(values = colours)
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
  # Automatic axes may insert layers before genes; select the actual fill output.
  fills <- unlist(lapply(built$data, function(layer) layer$gene_fill))
  expect_setequal(fills, unname(colours))
})

test_that("three label layouts retain representative annotations", {
  p <- ggchord(seqs, gene_data = genes) + geom_seq() + geom_gene()
  for (mode in c("radial", "auto", "arc")) {
    labels <- get_chord_layout(p + geom_gene_label_repel(label_layout = mode))$gene_labels
    expect_setequal(labels$label, genes$anno)
  }
})

test_that("feature, region and point-link layers retain source rows", {
  p <- ggchord(seqs) + geom_seq() + geom_feature(data = genes) +
    geom_seq_region(data = genes) +
    geom_link_line(data = data.frame(qaccver = "A", saccver = "B", qpos = 100, spos = 400))
  out <- export_ggchord_layout(p, include = c("feature", "link"))
  expect_setequal(out$feature$source_row, 1:2)
  expect_setequal(out$feature$source_row[out$feature$component == "seq_region"], 1:2)
  expect_equal(unique(out$link$source_row), 1L)
})

test_that("only Collinear uses automatic cross-row ribbon styles", {
  s <- data.frame(accver = LETTERS[1:4], length = 1000)
  h <- data.frame(qaccver = "A", saccver = c("B", "C", "D"),
    qstart = 100, qend = 200, sstart = 300, send = 400)
  coords <- list(coord_chord(), coord_collinear(), coord_radial_synteny(), coord_sankey())
  for (i in seq_along(coords)) {
    p <- ggchord(s, h) + geom_seq() +
      geom_link_ribbon(adjacent_only = FALSE) +
      geom_link_ribbon(adjacent_only = FALSE, linetype = "dotted") + coords[[i]]
    built <- ggplot2::ggplot_build(p)
    ribbons <- Filter(function(d) "ribbon_fill" %in% names(d), built$data)
    styles <- function(d) unname(vapply(split(d$linetype, d$group), unique, character(1)))
    expect_equal(styles(ribbons[[1]]), if (i == 2) c("solid", "dashed", "dashed") else rep("solid", 3))
    expect_equal(styles(ribbons[[2]]), rep("dotted", 3))
  }
})

test_that("empty primer results compose with circular layers", {
  hits <- find_primer_bindings(c(A = "AAAAAAAAAAAAAAAAAAAA"),
    primers = c(test = "ACGTTGCA"), min_annealed_bases = 8)
  p <- ggchord(data.frame(accver = "A", length = 20)) + geom_seq() +
    geom_primer(data = hits) + geom_primer_label_repel(data = hits) + coord_circular()
  out <- export_ggchord_layout(p, include = c("gene", "labels"))
  expect_equal(nrow(out$gene), 0L)
  expect_equal(nrow(out$labels), 0L)
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
