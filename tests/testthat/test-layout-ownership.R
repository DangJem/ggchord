test_that("derived plots export their own layouts and retain their cached base", {
  seq <- data.frame(accver = c("A", "B", "C"), length = 1000)
  ribbon <- data.frame(qaccver = c("A", "A", "B"),
    saccver = c("B", "C", "C"), qstart = 100, qend = 200,
    sstart = 300, send = 400)
  base <- ggchord(seq, ribbon, validate = "none") + geom_seq()
  invisible(get_chord_layout(base))
  adjacent <- base + geom_link_ribbon() + coord_collinear()
  all <- base + geom_link_ribbon(adjacent_only = FALSE) + coord_collinear()

  first <- export_ggchord_layout(adjacent, include = "ribbon")
  second <- export_ggchord_layout(all, include = "ribbon")
  repeated <- export_ggchord_layout(adjacent, include = "ribbon")
  expect_setequal(first$ribbon$source_row, c(1L, 3L))
  expect_setequal(second$ribbon$source_row, 1:3)
  expect_setequal(repeated$ribbon$source_row, first$ribbon$source_row)
  expect_setequal(names(get_chord_layout(base, build = FALSE)$seq_arcs), seq$accver)
})

test_that("build then export retains independently reused layers", {
  s <- data.frame(accver = "A", length = 1000)
  d <- data.frame(accver = "A", start = 100, end = 300,
    strand = "+", anno = "one")
  shared <- geom_gene(data = d)
  p <- ggchord(s) + geom_seq() + shared + shared
  invisible(ggplot2::ggplot_build(p))
  out <- export_ggchord_layout(p, include = c("seq", "gene"))
  expect_gt(nrow(out$seq), 0L)
  expect_gt(nrow(out$gene), 0L)
  expect_length(unique(out$gene$layer_id), 2L)
  expect_setequal(out$gene$source_row, 1L)
})
