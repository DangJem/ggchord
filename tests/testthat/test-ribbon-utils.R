# v0.8.0 tests: filter_ggchord_ribbons(), deduplicate_ggchord_ribbons(),
#             merge_ggchord_ribbons()

test_that("filter_ggchord_ribbons applies filters and reports", {
  data(ribbon_data_example)
  out <- filter_ggchord_ribbons(ribbon_data_example, min_pident = 95)
  expect_true(all(out$data$pident >= 95))
  expect_equal(out$report$n_kept + out$report$n_removed,
               out$report$n_input)
  expect_true("min_pident" %in% out$report$removed_by_reason$reason)
  # source_rows attribute tracks original rows
  expect_equal(nrow(out$data), length(attr(out$data, "source_rows")))
  expect_true(all(attr(out$data, "source_rows") %in% seq_len(nrow(ribbon_data_example))))
})

test_that("filter_ggchord_ribbons filters by pairs undirected", {
  data(ribbon_data_example)
  pairs <- data.frame(
    q = c("MT108731.1", "OQ646790.1"),
    s = c("MT118296.1", "OR222515.1"))
  out <- filter_ggchord_ribbons(ribbon_data_example, keep_pairs = pairs)
  expect_true(all(
    (out$data$qaccver %in% pairs$q & out$data$saccver %in% pairs$s) |
      (out$data$qaccver %in% pairs$s & out$data$saccver %in% pairs$q)))
  # list form with reversed direction works too
  out2 <- filter_ggchord_ribbons(
    ribbon_data_example,
    keep_pairs = list(c("MT118296.1", "MT108731.1")))
  expect_true(all(
    (out2$data$qaccver == "MT108731.1" & out2$data$saccver == "MT118296.1") |
      (out2$data$qaccver == "MT118296.1" & out2$data$saccver == "MT108731.1")))
})

test_that("filter_ggchord_ribbons drops self links and sorts", {
  data(ribbon_data_example)
  with_self <- rbind(ribbon_data_example,
                     transform(ribbon_data_example[1, ],
                               saccver = qaccver))
  out <- filter_ggchord_ribbons(with_self)
  expect_true(all(out$data$qaccver != out$data$saccver))
  expect_true("self_link" %in% out$report$removed_by_reason$reason)

  out2 <- filter_ggchord_ribbons(
    ribbon_data_example, sort_by = c("pident", "-length"))
  expect_true(all(diff(out2$data$pident) >= 0))
  # within equal pident, length descending
  pid <- out2$data$pident
  for (v in unique(pid)) {
    lens <- out2$data$length[pid == v]
    expect_true(all(diff(lens) <= 0))
  }
})

test_that("filter_ggchord_ribbons errors on missing columns", {
  data(ribbon_data_example)
  no_bitscore <- ribbon_data_example[
    , setdiff(colnames(ribbon_data_example), "bitscore"), drop = FALSE]
  expect_error(
    filter_ggchord_ribbons(no_bitscore, min_bitscore = 10),
    "column 'bitscore' is required")
  no_slen <- ribbon_data_example[
    , setdiff(colnames(ribbon_data_example), "slen"), drop = FALSE]
  expect_error(
    filter_ggchord_ribbons(no_slen, min_subject_coverage = 50),
    "column 'slen' is required")
})

test_that("deduplicate_ggchord_ribbons removes exact duplicates", {
  data(ribbon_data_example)
  dup <- rbind(ribbon_data_example, ribbon_data_example[1, ])
  out <- deduplicate_ggchord_ribbons(dup, by = "exact")
  expect_equal(out$report$n_removed, 1)
  expect_equal(nrow(out$data), nrow(ribbon_data_example))
  expect_equal(out$report$removed$row, 32)
})

test_that("deduplicate_ggchord_ribbons handles coordinates and overlap", {
  rb <- data.frame(
    qaccver = c("A", "A", "A"), saccver = c("B", "B", "B"),
    length = c(100, 100, 60), pident = c(90, 95, 80),
    qstart = c(1, 3, 20), qend = c(100, 102, 79),
    sstart = c(1, 3, 20), send = c(100, 102, 79))
  out <- deduplicate_ggchord_ribbons(rb, by = "coordinates", tolerance = 5)
  expect_equal(out$report$n_removed, 1)  # rows 1 and 2 are within 5 bp
  # best_pident keeps row 2 (95) and row 3 (80) stays
  expect_equal(out$data$pident, c(95, 80))

  out2 <- deduplicate_ggchord_ribbons(rb, by = "overlap",
                                      min_reciprocal_overlap = 0.9)
  # row 3 overlaps row 1 by 60/60 on both sequences -> duplicate
  expect_true(any(out2$report$removed$row == 3))
})

test_that("deduplicate_ggchord_ribbons keep = 'first'/'longest' works", {
  rb <- data.frame(
    qaccver = c("A", "A"), saccver = c("B", "B"),
    length = c(100, 50), pident = c(90, 99),
    qstart = c(1, 1), qend = c(100, 100),
    sstart = c(1, 1), send = c(100, 100))
  out <- deduplicate_ggchord_ribbons(rb, by = "exact", keep = "first")
  expect_equal(nrow(out$data), 1)
  expect_equal(out$data$pident, 90)
  out2 <- deduplicate_ggchord_ribbons(rb, by = "exact", keep = "longest")
  expect_equal(out2$data$pident, 90)
  out3 <- deduplicate_ggchord_ribbons(rb, by = "exact", keep = "best_pident")
  expect_equal(out3$data$pident, 99)
})

test_that("merge_ggchord_ribbons merges adjacent collinear blocks", {
  rb <- data.frame(
    qaccver = c("A", "A"), saccver = c("B", "B"),
    length = c(100, 100), pident = c(95, 97),
    qstart = c(1, 101), qend = c(100, 200),
    sstart = c(501, 601), send = c(600, 700))
  out <- merge_ggchord_ribbons(rb, max_gap = 0)
  expect_equal(nrow(out$data), 1)
  expect_equal(out$data$qstart, 1)
  expect_equal(out$data$qend, 200)
  expect_equal(out$data$sstart, 501)
  expect_equal(out$data$send, 700)
  expect_equal(out$data$length, 200)
  expect_equal(out$data$pident, 96)  # length-weighted mean
  expect_equal(out$report$n_merged, 2)
  expect_equal(out$report$from_rows, "1,2")
})

test_that("merge_ggchord_ribbons respects gaps and orientation", {
  rb <- data.frame(
    qaccver = c("A", "A", "A"), saccver = c("B", "B", "B"),
    length = c(100, 100, 100), pident = c(95, 95, 95),
    qstart = c(1, 201, 1), qend = c(100, 300, 100),
    sstart = c(501, 701, 1), send = c(600, 800, 100))
  # row3 is inverted relative to row1 (same q interval, subject runs opposite)
  out <- merge_ggchord_ribbons(rb, max_gap = 0)
  expect_equal(nrow(out$data), 3)  # nothing merges: gap and orientation differ
  # with a gap allowance, rows 1-2 merge
  out2 <- merge_ggchord_ribbons(rb[1:2, ], max_gap = 100)
  expect_equal(nrow(out2$data), 1)
  # without orientation requirement, row1 and row3 would be considered same
  # orientation? no: they differ, so even require_same_orientation = FALSE they
  # are not adjacent on the subject (interval_gap) - assert unmerged
  out3 <- merge_ggchord_ribbons(rb, max_gap = 0,
                                require_same_orientation = FALSE)
  expect_equal(nrow(out3$data), 3)
})

test_that("merge_ggchord_ribbons honours min_pident_difference", {
  rb <- data.frame(
    qaccver = c("A", "A"), saccver = c("B", "B"),
    length = c(100, 100), pident = c(80, 99),
    qstart = c(1, 101), qend = c(100, 200),
    sstart = c(501, 601), send = c(600, 700))
  out <- merge_ggchord_ribbons(rb, max_gap = 0, min_pident_difference = 5)
  expect_equal(nrow(out$data), 2)  # pident difference 19 > 5 -> no merge
})

test_that("ribbon utility output can be plotted directly", {
  data(seq_data_example)
  data(ribbon_data_example)
  f <- filter_ggchord_ribbons(ribbon_data_example, min_pident = 90)
  d <- deduplicate_ggchord_ribbons(f$data, by = "exact")
  m <- merge_ggchord_ribbons(d$data)
  p <- ggchord(seq_data_example, m$data) + geom_seq() + geom_ribbon()
  expect_no_error(ggplot_build(p))
})
