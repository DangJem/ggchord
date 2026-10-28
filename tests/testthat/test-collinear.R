test_that("collinear coordinates reuse public sequence, ribbon and gene layers", {
  seq <- data.frame(accver = c("A", "B"), length = c(1000, 800))
  ribbon <- data.frame(qaccver = "A", saccver = "B",
    qstart = 100, qend = 300, sstart = 400, send = 200)
  gene <- data.frame(accver = c("A", "B"), start = c(100, 200),
    end = c(300, 400), strand = c("+", "-"), anno = c("a", "b"))
  p <- ggchord(seq, ribbon, gene, validate = "none") +
    geom_seq(seq_orientation = c(A = 1, B = -1)) +
    geom_link_ribbon() + geom_gene(position = position_strand()) +
    geom_seq_label() + coord_collinear()
  old_dir <- setwd(tempdir())
  built <- tryCatch(ggplot2::ggplot_build(p), finally = setwd(old_dir))
  expect_s3_class(built, "ggplot_built")
  exported <- export_ggchord_layout(p, include = c("seq", "ribbon", "gene"))
  expect_identical(exported$metadata$coordinate, "collinear")
  expect_setequal(exported$sequence_layout$accver, seq$accver)
  expect_setequal(exported$gene$source_row, seq_len(nrow(gene)))
  expect_gt(nrow(exported$ribbon), 0L)
  expect_error(
    ggchord(seq) + geom_seq(seq_radius = 1) + coord_collinear(),
    "seq_radius"
  )
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
