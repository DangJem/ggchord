# v0.7.0 tests: clean_ggchord_data()

test_that("clean_ggchord_data returns all components and never modifies input", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  seq_in <- seq_data_example
  rib_in <- ribbon_data_example
  gen_in <- gene_data_example
  out <- clean_ggchord_data(seq_in, rib_in, gen_in)
  expect_named(out, c("seq_data", "ribbon_data", "gene_data", "report"))
  expect_identical(out$seq_data, seq_in)
  expect_identical(out$ribbon_data, rib_in)
  expect_identical(out$gene_data, gen_in)
  expect_identical(rib_in, ribbon_data_example)  # input untouched
  expect_equal(nrow(out$report), 0)  # clean data -> empty report
})

test_that("unknown_id policy drops or errors with clear messages", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- ribbon_data_example
  bad$saccver[1] <- "UNKNOWN"
  out <- clean_ggchord_data(seq_data_example, bad, unknown_id = "drop")
  expect_equal(nrow(out$ribbon_data), nrow(bad) - 1)
  expect_true(all(c("ribbon", which(bad$saccver == "UNKNOWN"), "saccver",
                    "unknown sequence ID", "drop") %in%
                    unlist(out$report[1, ])))
  expect_error(clean_ggchord_data(seq_data_example, bad, unknown_id = "error"),
               "unknown_id = 'drop' or 'keep'")
  out2 <- clean_ggchord_data(seq_data_example, bad, unknown_id = "keep")
  expect_equal(nrow(out2$ribbon_data), nrow(bad))
  expect_true(any(out2$report$action == "keep"))
})

test_that("out_of_range clips coordinates to [1, sequence length]", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, qstart = 0L, qend = qend + 1e6)
  out <- clean_ggchord_data(seq_data_example, bad, out_of_range = "clip")
  lens <- stats::setNames(seq_data_example$length, seq_data_example$seq_id)
  expect_true(all(out$ribbon_data$qstart >= 1))
  expect_true(all(out$ribbon_data$qend <=
                    lens[out$ribbon_data$qaccver]))
  expect_true(any(out$report$action == "clip"))
  expect_error(clean_ggchord_data(seq_data_example, bad, out_of_range = "error"),
               "out_of_range")
})

test_that("degenerate intervals after clipping are dropped with a report entry", {
  data(seq_data_example)
  # interval entirely outside the sequence: clipping leaves zero length
  rb <- data.frame(
    qaccver = "MT108731.1", saccver = "MT118296.1",
    length = 100, pident = 90,
    qstart = 1e6, qend = 1e6 + 99, sstart = 1, send = 100
  )
  out <- clean_ggchord_data(seq_data_example, rb, out_of_range = "clip")
  expect_equal(nrow(out$ribbon_data), 0)
  expect_true(any(grepl("degenerate", out$report$reason)))
})

test_that("reversed_interval sorts endpoints and records the original direction", {
  data(seq_data_example)
  data(ribbon_data_example)
  rev <- data.frame(
    qaccver = ribbon_data_example$qaccver,
    saccver = ribbon_data_example$saccver,
    length = ribbon_data_example$length,
    pident = ribbon_data_example$pident,
    qstart = ribbon_data_example$qend,
    qend = ribbon_data_example$qstart,
    sstart = ribbon_data_example$sstart,
    send = ribbon_data_example$send
  )
  out <- clean_ggchord_data(seq_data_example, rev, reversed_interval = "sort")
  expect_true(all(out$ribbon_data$qstart <= out$ribbon_data$qend))
  expect_true(any(out$report$action == "sort"))
  expect_true(any(grepl("original direction", out$report$reason)))

  out2 <- clean_ggchord_data(seq_data_example, rev, reversed_interval = "drop")
  expect_equal(nrow(out2$ribbon_data), 0)
  expect_true(any(out2$report$action == "drop"))
  expect_error(clean_ggchord_data(seq_data_example, rev,
                                  reversed_interval = "error"),
               "reversed_interval")
})

test_that("invalid_pident clips to [0, 100]", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, pident = pident + 200)
  out <- clean_ggchord_data(seq_data_example, bad, invalid_pident = "clip")
  expect_true(all(out$ribbon_data$pident <= 100))
  expect_true(any(out$report$action == "clip"))
  out2 <- clean_ggchord_data(seq_data_example, bad, invalid_pident = "drop")
  expect_equal(nrow(out2$ribbon_data), 0)
})

test_that("empty_annotation replaces or drops", {
  data(seq_data_example)
  data(gene_data_example)
  bad <- transform(gene_data_example, anno = NA_character_)
  out <- clean_ggchord_data(seq_data_example, gene_data = bad,
                            empty_annotation = "replace",
                            replacement_annotation = "unannotated")
  expect_true(all(out$gene_data$anno == "unannotated"))
  expect_true(any(out$report$action == "replace"))
  out2 <- clean_ggchord_data(seq_data_example, gene_data = bad,
                             empty_annotation = "drop")
  expect_equal(nrow(out2$gene_data), 0)
  out3 <- clean_ggchord_data(seq_data_example, gene_data = bad,
                             empty_annotation = "keep")
  expect_equal(nrow(out3$gene_data), nrow(bad))
})

test_that("gene coordinates are clipped and reversed genes are sorted", {
  data(seq_data_example)
  data(gene_data_example)
  bad <- transform(gene_data_example, start = 0, end = end + 1e6)
  out <- clean_ggchord_data(seq_data_example, gene_data = bad,
                            out_of_range = "clip", reversed_interval = "sort")
  expect_true(all(out$gene_data$start >= 1))
  lens <- stats::setNames(seq_data_example$length, seq_data_example$seq_id)
  expect_true(all(out$gene_data$end <= lens[out$gene_data$seq_id]))
})

test_that("cleaned data can be plotted directly", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  bad_r <- transform(ribbon_data_example, qstart = 0L)
  bad_g <- transform(gene_data_example, anno = NA_character_)
  out <- clean_ggchord_data(seq_data_example, bad_r, bad_g)
  # default empty_annotation = "keep" leaves missing annos -> expected warning
  p <- suppressWarnings(ggchord(out$seq_data, out$ribbon_data, out$gene_data) +
                          geom_seq() + geom_ribbon() + geom_gene())
  expect_no_error(ggplot_build(p))
  # with explicit replacement the cleaned data is warning-free
  out2 <- clean_ggchord_data(seq_data_example, bad_r, bad_g,
                             empty_annotation = "replace")
  expect_no_warning(
    ggchord(out2$seq_data, out2$ribbon_data, out2$gene_data))
})

test_that("cleaned data has a print method and ggchord_clean class", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, qstart = 0L)
  out <- clean_ggchord_data(seq_data_example, bad, out_of_range = "clip")
  expect_s3_class(out, "ggchord_clean")
  expect_output(print(out), "ggchord cleaned data")
})
