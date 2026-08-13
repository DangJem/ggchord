# v0.7.0 tests: validate_ggchord_data() + print/summary methods

test_that("valid example data passes validation", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)
  v <- validate_ggchord_data(seq_data_example, ribbon_data_example,
                             gene_data_example)
  expect_s3_class(v, "ggchord_validation")
  expect_true(v$valid)
  expect_equal(nrow(v$errors), 0)
  expect_equal(nrow(v$warnings), 0)
  expect_equal(v$data_summary$n_seq, nrow(seq_data_example))
  expect_equal(v$data_summary$n_ribbon, nrow(ribbon_data_example))
  expect_equal(v$data_summary$n_gene, nrow(gene_data_example))
})

test_that("seq_data structural problems are severe errors", {
  expect_false(validate_ggchord_data("not a df")$valid)
  expect_false(validate_ggchord_data(data.frame())$valid)
  expect_false(validate_ggchord_data(data.frame(id = "a", len = 1))$valid)
  expect_false(validate_ggchord_data(
    data.frame(seq_id = c("a", NA), length = c(1, 2)))$valid)
  expect_false(validate_ggchord_data(
    data.frame(seq_id = c("a", ""), length = c(1, 2)))$valid)
  v <- validate_ggchord_data(data.frame(seq_id = c("a", "a"), length = c(1, 2)))
  expect_false(v$valid)
  expect_true("duplicate_id" %in% v$errors$category)
  expect_equal(sort(v$invalid_rows$duplicate_id), c(1, 2))
})

test_that("seq_data length problems are detected", {
  expect_false(validate_ggchord_data(
    data.frame(seq_id = "a", length = NA_real_))$valid)
  expect_false(validate_ggchord_data(
    data.frame(seq_id = "a", length = Inf))$valid)
  expect_false(validate_ggchord_data(
    data.frame(seq_id = "a", length = 0))$valid)
  expect_false(validate_ggchord_data(
    data.frame(seq_id = "a", length = -5))$valid)
  v <- validate_ggchord_data(data.frame(seq_id = "a", length = 10.5))
  expect_true(v$valid)  # fractional length is only a warning
  expect_true("fractional_length" %in% v$warnings$category)
})

test_that("extreme length differences are a warning", {
  v <- validate_ggchord_data(data.frame(
    seq_id = c("tiny", "huge"), length = c(1, 10000)))
  expect_true(v$valid)
  expect_true("extreme_length" %in% v$warnings$category)
})

test_that("ribbon_data required columns and numeric types are errors", {
  data(seq_data_example)
  expect_false(validate_ggchord_data(seq_data_example, "nope")$valid)
  bad <- data.frame(qaccver = "a", saccver = "b")
  v <- validate_ggchord_data(seq_data_example, bad)
  expect_false(v$valid)
  expect_true("missing_columns" %in% v$errors$category)
})

test_that("ribbon_data non-finite values and unknown IDs are errors with rows", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, qstart = NA_real_)
  v <- validate_ggchord_data(seq_data_example, bad)
  expect_false(v$valid)
  expect_true("non_finite" %in% v$errors$category)
  expect_equal(v$invalid_rows$non_finite, which(is.na(bad$qstart)))

  bad2 <- transform(ribbon_data_example, saccver = "UNKNOWN_SEQ")
  v2 <- validate_ggchord_data(seq_data_example, bad2)
  expect_false(v2$valid)
  expect_true("unknown_id" %in% v2$errors$category)
  expect_equal(v2$invalid_rows$unknown_id, which(bad2$saccver == "UNKNOWN_SEQ"))
  expect_gt(v2$data_summary$n_unknown_ribbon_ids, 0)
})

test_that("ribbon self-links are errors and can be disabled", {
  data(seq_data_example)
  rb <- data.frame(
    qaccver = c("A", "A"), saccver = c("A", "B"),
    length = c(100, 100), pident = c(90, 90),
    qstart = c(1, 1), qend = c(100, 100),
    sstart = c(1, 1), send = c(100, 100)
  )
  seqs <- data.frame(seq_id = c("A", "B"), length = c(500, 500))
  v <- validate_ggchord_data(seqs, rb)
  expect_false(v$valid)
  expect_true("self_link" %in% v$errors$category)
  expect_equal(v$invalid_rows$self_link, 1)
  v2 <- validate_ggchord_data(seqs, rb, check_self_links = FALSE)
  expect_true(v2$valid)
})

test_that("ribbon out-of-range coordinates are errors and can be disabled", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, qend = qend + 1e6)
  v <- validate_ggchord_data(seq_data_example, bad)
  expect_false(v$valid)
  expect_true("out_of_range" %in% v$errors$category)
  v2 <- validate_ggchord_data(seq_data_example, bad, check_coordinates = FALSE)
  expect_true(v2$valid)
})

test_that("reversed ribbon intervals and bad pident are warnings", {
  data(seq_data_example)
  data(ribbon_data_example)
  rev <- transform(ribbon_data_example, qstart = qend + 1, qend = qstart)
  # rebuild a proper frame (transform evaluates simultaneously)
  rev <- data.frame(
    qaccver = ribbon_data_example$qaccver,
    saccver = ribbon_data_example$saccver,
    length = ribbon_data_example$length,
    pident = ribbon_data_example$pident,
    qstart = ribbon_data_example$qend + 1,
    qend = ribbon_data_example$qstart,
    sstart = ribbon_data_example$sstart,
    send = ribbon_data_example$send
  )
  v <- validate_ggchord_data(seq_data_example, rev)
  expect_true(v$valid)
  expect_true("reversed_interval" %in% v$warnings$category)

  pid <- transform(ribbon_data_example, pident = pident + 200)
  v2 <- validate_ggchord_data(seq_data_example, pid)
  expect_true(v2$valid)
  expect_true("invalid_pident" %in% v2$warnings$category)
})

test_that("ribbon duplicates and overlap are warnings", {
  data(seq_data_example)
  data(ribbon_data_example)
  dup <- rbind(ribbon_data_example, ribbon_data_example[1, ])
  v <- validate_ggchord_data(seq_data_example, dup)
  expect_true(v$valid)
  expect_true("exact_duplicate" %in% v$warnings$category)

  near <- rbind(ribbon_data_example, transform(ribbon_data_example[1, ],
                                               qstart = qstart + 2,
                                               qend = qend + 2))
  v2 <- validate_ggchord_data(seq_data_example, near)
  expect_true("near_duplicate" %in% v2$warnings$category)

  v3 <- validate_ggchord_data(seq_data_example, dup, check_duplicates = FALSE)
  expect_false("exact_duplicate" %in% v3$warnings$category)
})

test_that("optional BLAST fields give diagnostics when invalid", {
  data(seq_data_example)
  rb <- data.frame(
    qaccver = "A", saccver = "B", length = 100, pident = 90,
    qstart = 1, qend = 100, sstart = 1, send = 100,
    evalue = -1, bitscore = -5, qcovs = 150, qlen = 0, slen = 0
  )
  seqs <- data.frame(seq_id = c("A", "B"), length = c(500, 500))
  v <- validate_ggchord_data(seqs, rb)
  expect_true(v$valid)
  cats <- unique(v$warnings$category)
  expect_true(all(c("evalue_diag", "bitscore_diag", "qcovs_diag",
                    "qlen_slen_diag") %in% cats))
})

test_that("gene_data problems are detected", {
  data(seq_data_example)
  data(gene_data_example)
  # missing columns
  v <- validate_ggchord_data(seq_data_example, gene_data =
                               gene_data_example[, c("seq_id", "start")])
  expect_false(v$valid)
  expect_true("missing_columns" %in% v$errors$category)

  # non-finite coordinates
  v2 <- validate_ggchord_data(seq_data_example, gene_data =
                                transform(gene_data_example, start = Inf))
  expect_false(v2$valid)

  # invalid strand
  v3 <- validate_ggchord_data(seq_data_example, gene_data =
                                transform(gene_data_example, strand = "x"))
  expect_false(v3$valid)
  expect_true("invalid_strand" %in% v3$errors$category)

  # unknown seq_id -> warning
  v4 <- validate_ggchord_data(seq_data_example, gene_data =
                                transform(gene_data_example,
                                          seq_id = "UNKNOWN_GENE"))
  expect_true(v4$valid)
  expect_true("unknown_id" %in% v4$warnings$category)
  expect_gt(v4$data_summary$n_unknown_gene_ids, 0)

  # out-of-range -> error
  v5 <- validate_ggchord_data(seq_data_example, gene_data =
                                transform(gene_data_example, end = end + 1e6))
  expect_false(v5$valid)
  expect_true("out_of_range" %in% v5$errors$category)

  # reversed -> warning
  v6 <- validate_ggchord_data(seq_data_example, gene_data =
                                transform(gene_data_example,
                                          start = end, end = start))
  expect_true(v6$valid)
  expect_true("reversed_interval" %in% v6$warnings$category)

  # empty annotation -> warning
  v7 <- validate_ggchord_data(seq_data_example, gene_data =
                                transform(gene_data_example,
                                          anno = NA_character_))
  expect_true(v7$valid)
  expect_true("missing_anno" %in% v7$warnings$category)

  # duplicates -> warning
  v8 <- validate_ggchord_data(seq_data_example, gene_data =
                                rbind(gene_data_example, gene_data_example[1, ]))
  expect_true(v8$valid)
  expect_true("exact_duplicate" %in% v8$warnings$category)
})

test_that("strict = TRUE stops on errors but not on warnings only", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, saccver = "UNKNOWN")
  expect_error(validate_ggchord_data(seq_data_example, bad, strict = TRUE),
               "strict = FALSE")

  # warnings only: no stop
  pid <- transform(ribbon_data_example, pident = pident + 200)
  expect_no_error(validate_ggchord_data(seq_data_example, pid, strict = TRUE))
})

test_that("print and summary methods work", {
  data(seq_data_example)
  v <- validate_ggchord_data(seq_data_example)
  expect_output(print(v), "VALID")
  s <- summary(v)
  expect_s3_class(s, "ggchord_validation_summary")
  expect_equal(nrow(s), 0)

  bad <- validate_ggchord_data(data.frame(seq_id = c("a", "a"),
                                          length = c(1, 2)))
  expect_output(print(bad), "INVALID")
  s2 <- summary(bad)
  expect_true("duplicate_id" %in% s2$category)
  expect_equal(s2$severity[s2$category == "duplicate_id"], "error")
  expect_output(print(s2), "Valid: FALSE")
})

test_that("cleanable lists the fixable issues", {
  data(seq_data_example)
  data(ribbon_data_example)
  bad <- transform(ribbon_data_example, saccver = "UNKNOWN")
  v <- validate_ggchord_data(seq_data_example, bad)
  expect_true("unknown_id" %in% v$cleanable$category)
  expect_true(grepl("clean_ggchord_data", v$cleanable$suggestion[1]))
})
