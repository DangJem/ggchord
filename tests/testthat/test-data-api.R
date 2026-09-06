test_that("validation and cleaning return their public result objects", {
  data(seq_data_example)
  data(ribbon_data_example)
  data(gene_data_example)

  validation <- validate_ggchord_data(
    seq_data_example, ribbon_data_example, gene_data_example
  )
  expect_s3_class(validation, "ggchord_validation")

  cleaned <- clean_ggchord_data(
    seq_data_example, ribbon_data_example, gene_data_example
  )
  expect_s3_class(cleaned, "ggchord_clean")
  expect_true(all(c("seq_data", "ribbon_data", "gene_data", "report") %in%
                    names(cleaned)))

  malformed <- transform(
    ribbon_data_example[1, ], length = "bad", pident = "bad",
    qstart = "bad", qend = "bad", sstart = "bad", send = "bad"
  )
  malformed_result <- validate_ggchord_data(
    seq_data_example, malformed, check_duplicates = TRUE
  )
  expect_s3_class(malformed_result, "ggchord_validation")
  expect_true(any(malformed_result$errors$category == "non_numeric"))
  expect_error(
    ggchord(seq_data_example, malformed, validate = "none"),
    "ribbon_data\\$length must be numeric"
  )
})

test_that("packaged gene example is compact, distributed and strand-balanced", {
  data(seq_data_example)
  data(gene_data_example)
  expect_equal(
    as.integer(table(factor(
      gene_data_example$accver, levels = seq_data_example$accver
    ))),
    rep(8L, nrow(seq_data_example))
  )
  strand_counts <- table(gene_data_example$accver, gene_data_example$strand)
  expect_true(all(strand_counts[, "+"] == 4L))
  expect_true(all(strand_counts[, "-"] == 4L))
  expect_true("source_strand" %in% names(gene_data_example))
  expect_true(all(gene_data_example$start < gene_data_example$end))
})

test_that("packaged ribbon example is balanced across sequence pairs", {
  data(ribbon_data_example)
  pair <- paste(
    pmin(ribbon_data_example$qaccver, ribbon_data_example$saccver),
    pmax(ribbon_data_example$qaccver, ribbon_data_example$saccver),
    sep = "\r"
  )
  counts <- table(pair)
  expect_equal(nrow(ribbon_data_example), 13L)
  expect_equal(length(counts), 6L)
  expect_lte(max(counts), 3L)
  expect_true(all(ribbon_data_example$length >= 300))
})

test_that("the three import helpers parse minimal files", {
  fasta <- tempfile(fileext = ".fna")
  writeLines(c(">seqA", "ACGTACGT"), fasta)
  expect_equal(read_fasta_lengths(fasta)$length, 8)

  blast <- tempfile(fileext = ".tsv")
  writeLines(paste(
    "seqA", "seqB", "98", "100", "0", "0",
    "1", "100", "1", "100", "1e-20", "200",
    sep = "\t"
  ), blast)
  expect_equal(nrow(read_blast(blast, format = "outfmt6")), 1)

  gff <- tempfile(fileext = ".gff3")
  writeLines(paste(
    "seqA", "source", "CDS", "1", "100", ".", "+", "0",
    "ID=cds1;product=test%20protein",
    sep = "\t"
  ), gff)
  expect_equal(read_gff3(gff)$anno, "test protein")
})

test_that("ribbon preparation helpers run on simple inputs", {
  data(ribbon_data_example)

  filtered <- filter_ggchord_ribbons(
    ribbon_data_example,
    min_pident = 90
  )
  expect_true(all(filtered$data$pident >= 90))

  duplicated <- rbind(ribbon_data_example, ribbon_data_example[1, ])
  deduplicated <- deduplicate_ggchord_ribbons(duplicated)
  expect_equal(nrow(deduplicated$data), nrow(ribbon_data_example))

  blocks <- data.frame(
    qaccver = c("A", "A"),
    saccver = c("B", "B"),
    length = c(100, 100),
    pident = c(95, 97),
    qstart = c(1, 101),
    qend = c(100, 200),
    sstart = c(501, 601),
    send = c(600, 700)
  )
  expect_equal(nrow(merge_ggchord_ribbons(blocks)$data), 1)
})

test_that("dense ribbon helpers bundle explicitly and optimize deterministically", {
  seq <- data.frame(accver = c("A", "B", "C"), length = c(1000, 1000, 1000))
  ribbons <- data.frame(
    qaccver = c("A", "A", "A", "A", "B"),
    saccver = c("B", "B", "B", "B", "C"),
    length = rep(100, 5), pident = c(90, 94, 80, 84, 88),
    qstart = c(100, 120, 100, 120, 700),
    qend = c(199, 219, 199, 219, 799),
    sstart = c(100, 120, 219, 239, 100),
    send = c(199, 219, 120, 140, 199),
    note = c("same", "different", "reverse", "reverse", "single")
  )
  original <- ribbons

  bundled <- bundle_ggchord_ribbons(ribbons, seq, bins = 10)
  expect_identical(ribbons, original)
  expect_lt(nrow(bundled$data), nrow(ribbons))
  expect_equal(sort(bundled$data$.bundle_n), c(1L, 2L, 2L))
  expect_equal(sum(bundled$data$.bundle_weight), sum(ribbons$length))
  expect_true(is.na(bundled$data$note[bundled$data$pident == 92]))
  expect_equal(sort(bundled$report$n_input), c(1L, 2L, 2L))
  expect_setequal(bundled$report$source_rows, c("1,2", "3,4", "5"))

  identity_weighted <- bundle_ggchord_ribbons(
    ribbons, seq, bins = 10, weight = "pident"
  )
  expect_equal(
    sum(identity_weighted$data$.bundle_weight),
    sum(ribbons$length * ribbons$pident / 100)
  )

  optimized <- optimize_ggchord_layout(seq, ribbons)
  repeated <- optimize_ggchord_layout(seq, ribbons)
  expect_setequal(optimized$seq_order, seq$accver)
  expect_true(all(optimized$seq_orientation %in% c(-1, 1)))
  expect_lt(optimized$score_after, optimized$score_before)
  expect_identical(optimized, repeated)

  many <- ribbons[rep(seq_len(nrow(ribbons)), length.out = 2001), ]
  expect_true(optimize_ggchord_layout(seq, many)$report$approximate)
})

test_that("cleaning preserves kept unknown genes and ribbon direction", {
  seq <- data.frame(accver = c("A", "B"), length = c(100, 100))
  genes <- data.frame(
    accver = "unknown", start = 1, end = 10, strand = "+", anno = "x"
  )
  ribbons <- data.frame(
    qaccver = "A", saccver = "B", length = 10, pident = 90,
    qstart = 1, qend = 10, sstart = 20, send = 11
  )

  cleaned <- clean_ggchord_data(
    seq, ribbons, genes, unknown_id = "keep", reversed_interval = "sort"
  )
  expect_equal(nrow(cleaned$gene_data), 1)
  expect_equal(cleaned$ribbon_data$direction, "reverse")
  expect_lt(cleaned$ribbon_data$sstart, cleaned$ribbon_data$send)

  validation <- validate_ggchord_data(seq, gene_data = genes, strict = FALSE)
  expect_false(validation$valid)
  expect_true(any(validation$errors$category == "unknown_id"))
})

test_that("ribbon reports retain all reasons and input-first semantics", {
  ribbons <- data.frame(
    qaccver = c("A", "A"), saccver = c("A", "B"),
    length = c(5, 5), pident = c(10, 20),
    qstart = c(50, 1), qend = c(55, 5),
    sstart = c(50, 1), send = c(55, 5), score = c("x", "y")
  )
  filtered <- filter_ggchord_ribbons(
    ribbons, min_pident = 50, min_length = 10, drop_self_links = TRUE
  )
  expect_equal(nrow(filtered$report$removed_reasons), 5)

  dup <- ribbons[c(1, 1), ]
  dup$qstart <- c(50, 1)
  dup$qend <- c(55, 6)
  dup$sstart <- c(50, 1)
  dup$send <- c(55, 6)
  dedup <- deduplicate_ggchord_ribbons(
    dup, by = "coordinates", tolerance = 100, keep = "first"
  )
  expect_equal(attr(dedup$data, "source_rows"), 1L)

  merged <- merge_ggchord_ribbons(
    transform(ribbons, qaccver = "A", saccver = "B",
              qstart = c(1, 6), qend = c(5, 10),
              sstart = c(1, 6), send = c(5, 10), pident = c(90, 90))
  )
  expect_true(is.na(merged$data$score))
})

test_that("outfmt7 fields, GFF3 FASTA boundaries and source files are parsed", {
  blast <- tempfile(fileext = ".o7")
  writeLines(c(
    "# BLASTN 2.15.0+",
    "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, q. start, q. end, s. start, s. end",
    "A\tB\t99\t10\t1\t10\t20\t11"
  ), blast)
  parsed <- read_blast(blast, format = "outfmt7", source_file = TRUE)
  expect_true(all(c("qaccver", "saccver", ".source_file") %in% names(parsed)))

  gff <- tempfile(fileext = ".gff3")
  writeLines(c(
    "##gff-version 3",
    "A\tsrc\tCDS\t1\t10\t.\t+\t0\tID=x",
    "##FASTA", ">A", "ACGT"
  ), gff)
  expect_equal(nrow(read_gff3(gff)), 1)
})

test_that("focus_ggchord_data synchronizes loci, genes and ribbon direction", {
  seq <- data.frame(accver = c("A", "B"), length = c(100, 100))
  ribbons <- data.frame(
    qaccver = "A", saccver = "B", length = 81, pident = 95,
    qstart = 10, qend = 90, sstart = 90, send = 10
  )
  genes <- data.frame(
    accver = "A", start = c(5, 30), end = c(25, 50),
    strand = "+", anno = c("edge", "inside")
  )
  loci <- data.frame(
    accver = c("A", "B"), start = c(20, 20), end = c(80, 80)
  )

  focused <- focus_ggchord_data(seq, ribbons, genes, loci, boundary = "trim")
  expect_equal(focused$seq_data$length, c(61, 61))
  expect_equal(
    focused$ribbon_data[c("qstart", "qend", "sstart", "send")],
    data.frame(qstart = 1, qend = 61, sstart = 61, send = 1)
  )
  expect_equal(focused$gene_data$start, c(1, 11))
  expect_equal(focused$gene_data$end, c(6, 31))
  expect_equal(focused$ribbon_data$.source_row, 1L)

  dropped <- focus_ggchord_data(seq, ribbons, genes, loci, boundary = "drop")
  expect_equal(nrow(dropped$ribbon_data), 0)
  expect_equal(dropped$gene_data$anno, "inside")
})
