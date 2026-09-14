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

test_that("primer geometry uses the shared feature engine", {
  sequence <- data.frame(accver = "circle", length = 1000)
  primers <- data.frame(
    accver = "circle", start = 100L, end = 119L,
    strand = "+", name = "test primer"
  )
  layer <- geom_primer(data = primers)
  expect_identical(layer$ggchord_params$feature_role, "primer")
  expect_true(layer$ggchord_params$is_primer)
  plot <- ggchord(sequence, validate = "none") + geom_seq() + layer +
    coord_circular()
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})
