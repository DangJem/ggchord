test_that("bundled primer catalogue preserves sequence aliases", {
  catalogue <- primer_catalog("all")
  aliases <- primer_catalog("all", aliases = TRUE)
  expect_equal(nrow(catalogue), 137L)
  expect_equal(nrow(aliases), 163L)
  duplicated_sequence <- aliases$sequence == "GGGAAACGCCTGGTATCTTT"
  expect_setequal(aliases$name[duplicated_sequence], c("hU6-F", "pBR322ori-F"))
})

test_that("primer search derives strand and circular origin", {
  primers <- c(test = "ACGTTGCA")
  target <- c(circle = "TTGCAAAAACG")
  hits <- find_primer_bindings(
    target, primers = primers, circular = TRUE,
    min_annealed_bases = 8L
  )
  expect_equal(nrow(hits), 1L)
  expect_true(hits$crosses_origin)
  expect_equal(hits$strand, "+")
  expect_equal(hits$aliases[[1L]], "test")
})

test_that("primer geometry uses a backbone arc and dedicated unboxed labels", {
  sequence <- data.frame(accver = "circle", length = 1000)
  primers <- data.frame(
    accver = "circle", start = 100L, end = 119L,
    strand = "+", name = "test primer"
  )
  layer <- geom_primer(data = primers)
  expect_identical(layer$ggchord_params$feature_role, "primer")
  expect_true(layer$ggchord_params$is_primer)
  expect_identical(layer$ggchord_params$feature_shape, "primer_arc")
  expect_false(isTRUE(layer$ggchord_params$feature_position$ggchord_feature_stack))
  plot <- ggchord(sequence, validate = "none") + geom_seq() + layer +
    geom_primer_label_repel(data = primers) +
    coord_circular()
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  registry <- export_ggchord_layout(plot)$annotation_registry
  expect_true(any(registry$kind == "primer_label"))
  expect_true(all(registry$annotation_class[registry$kind == "primer_label"] ==
    "primer"))
})

test_that("reference primer coordinates are one-based inclusive", {
  data(plasmid_example_pSB1C3)
  primers <- find_primer_bindings(plasmid_example_pSB1C3, set = "reference")
  expect_equal(primers$start, c(155L, 1931L))
  expect_equal(primers$end, c(174L, 1950L))
})
