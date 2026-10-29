test_that("collinear coordinates preserve sequence direction and source identity", {
  p <- ggchord(smoke_sequences(), smoke_ribbon(), smoke_genes()) +
    geom_seq(seq_orientation = c(A = 1, B = -1)) + geom_link_ribbon() +
    geom_gene(position = position_strand()) + coord_collinear()
  out <- export_ggchord_layout(p, include = c("seq", "ribbon", "gene"))
  expect_identical(out$metadata$coordinate, "collinear")
  expect_setequal(out$gene$source_row, 1:2)
  expect_equal(out$sequence_layout$orientation, c(1, -1))
})

test_that("collinear labels preserve feature stacking and their own text", {
  seq <- data.frame(accver = "A", length = 1000)
  features <- data.frame(accver = "A", start = c(100, 150), end = c(600, 700),
    strand = "+", type = "CDS", label = c("gene A", "gene B"))
  p <- ggchord(seq, validate = "none") + geom_seq() +
    geom_feature(data = features,
      position = position_feature_stack(base_position = position_strand(.2))) +
    coord_collinear()
  before <- export_ggchord_layout(p, include = "feature")
  after <- export_ggchord_layout(p + geom_feature_label_repel(data = features),
    include = c("feature", "labels"))
  placement <- c("source_row", "lane", "normal_offset")
  expect_equal(unique(after$feature[placement]), unique(before$feature[placement]))
  text <- after$labels[after$labels$.component %in% "text", , drop = FALSE]
  expect_equal(as.character(text$label), features$label)
  expect_equal(lapply(unique(text[placement]), unname),
    lapply(unique(after$feature[placement]), unname))
})

test_that("ribbon outline defaults follow display adjacency without changing mapped layers", {
  seq <- data.frame(accver = LETTERS[1:4], length = 1000)
  ribbon <- data.frame(qaccver = "A", saccver = c("B", "C", "D"),
    qstart = 100, qend = 200, sstart = c(400, 300, 600),
    send = c(300, 400, 500), kind = c("a", "b", "c"))
  base <- ggchord(seq, ribbon, validate = "none") + geom_seq()
  styles <- function(data) {
    value <- if ("ribbon_linetype" %in% names(data)) {
      data$ribbon_linetype
    } else data$linetype
    unname(vapply(split(as.character(value), data$group), unique, character(1L)))
  }
  linear <- ggplot2::ggplot_build(base +
    geom_link_ribbon(adjacent_only = FALSE) + coord_collinear())
  expect_equal(styles(linear$data[[2L]]), c("solid", "dashed", "dashed"))
  mapped <- ggplot2::ggplot_build(base + aes(linetype = "chosen") +
    geom_link_ribbon(adjacent_only = FALSE) +
    geom_link_ribbon(adjacent_only = FALSE, inherit.aes = TRUE) +
    ggplot2::scale_linetype_manual(values = c(chosen = "dotdash")) +
    coord_collinear())
  expect_equal(styles(mapped$data[[2L]]), c("solid", "dashed", "dashed"))
  expect_equal(styles(mapped$data[[3L]]), rep("dotdash", 3L))
})


test_that("narrow linear channels retain ribbon endpoints between backbones", {
  p <- ggchord(smoke_sequences(), smoke_ribbon()) +
    geom_seq(seq_gap = .01) + geom_link_ribbon() + coord_collinear()
  out <- export_ggchord_layout(p, include = "ribbon")
  expect_true(all(out$ribbon$y >= min(out$sequence_layout$y) &
    out$ribbon$y <= max(out$sequence_layout$y)))
})

test_that("tiny linear arrows retain non-overlapping biological intervals", {
  features <- data.frame(accver = "A", start = c(100, 110), end = c(102, 112),
    strand = "+", type = "CDS")
  p <- ggchord(data.frame(accver = "A", length = 1e6)) + geom_seq() +
    geom_feature(data = features, position = position_feature_stack()) +
    coord_collinear()
  out <- export_ggchord_layout(p, include = "feature")$feature
  expect_false(any(out$display_extended))
  expect_lt(max(out$x[out$source_row == 1]), min(out$x[out$source_row == 2]))
})
