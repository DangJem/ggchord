# v0.8.0 tests: read_blast(), read_gff3(), read_fasta_lengths()

test_that("read_blast parses outfmt7 17-column files", {
  f <- tempfile(fileext = ".o7")
  writeLines(c(
    "# BLASTN 2.13.0+",
    "# Query: seqA",
    paste("seqA", "seqB", "98.5", "1200", "18", "0", "1", "1200",
          "1", "1200", "1e-180", "2400", "100", "5000", "4800",
          "plus", "seqB description", sep = "\t"),
    paste("seqA", "seqC", "95.0", "800", "40", "2", "1300", "2100",
          "50", "850", "1e-100", "1500", "80", "5000", "6000",
          "minus", "seqC description", sep = "\t")
  ), f)
  rb <- read_blast(f)
  expect_equal(
    colnames(rb),
    c("qaccver", "saccver", "length", "pident", "qstart", "qend",
      "sstart", "send", "mismatch", "gapopen", "evalue", "bitscore",
      "qcovs", "qlen", "slen", "sstrand", "stitle"))
  expect_equal(nrow(rb), 2)
  expect_equal(rb$pident, c(98.5, 95.0))
  expect_true(is.numeric(rb$evalue))
})

test_that("read_blast parses outfmt6 12-column files", {
  f <- tempfile(fileext = ".o6")
  writeLines(paste("seqA", "seqB", "98.5", "1200", "18", "0",
                   "1", "1200", "1", "1200", "1e-180", "2400",
                   sep = "\t"), f)
  rb <- read_blast(f, format = "outfmt6")
  expect_equal(ncol(rb), 12)
  expect_equal(names(rb)[1:8],
               c("qaccver", "saccver", "length", "pident",
                 "qstart", "qend", "sstart", "send"))
})

test_that("read_blast reports missing files and unknown layouts clearly", {
  expect_error(read_blast(tempfile()), "file not found")
  f <- tempfile()
  writeLines("a\tb\tc", f)
  expect_error(read_blast(f), "auto-detect")
})

test_that("read_gff3 produces gene_data with extracted annotations", {
  f <- tempfile(fileext = ".gff3")
  writeLines(c(
    "##gff-version 3",
    paste("seqA", "source", "CDS", "101", "500", ".", "+", "0",
          "ID=cds1;product=hypothetical%20protein", sep = "\t"),
    paste("seqA", "source", "CDS", "600", "900", ".", "-", "0",
          "ID=cds2;Name=integrase", sep = "\t"),
    paste("seqA", "source", "tRNA", "1000", "1080", ".", "+", "0",
          "ID=trna1", sep = "\t")
  ), f)
  gd <- read_gff3(f)
  expect_equal(colnames(gd)[1:5],
               c("seq_id", "start", "end", "strand", "anno"))
  expect_equal(nrow(gd), 2)  # tRNA filtered out by default
  expect_equal(gd$anno, c("hypothetical protein", "integrase"))
  expect_equal(gd$strand, c("+", "-"))
  expect_true("type" %in% colnames(gd))
  expect_true("attributes" %in% colnames(gd))
})

test_that("read_gff3 handles unstranded features and feature_types", {
  f <- tempfile(fileext = ".gff3")
  writeLines(c(
    paste("seqA", "s", "CDS", "1", "100", ".", ".", "0", "ID=a", sep = "\t"),
    paste("seqA", "s", "tRNA", "200", "300", ".", "?", "0", "ID=b", sep = "\t"),
    paste("seqA", "s", "rRNA", "400", "500", ".", "+", "0", "ID=c", sep = "\t")
  ), f)
  gd <- read_gff3(f, feature_types = c("CDS", "tRNA", "rRNA"))
  expect_equal(gd$strand, c("+", "+", "+"))  # "." and "?" -> "+"
  gd2 <- read_gff3(f, feature_types = c("CDS", "tRNA", "rRNA"),
                   unstranded = "drop")
  expect_equal(nrow(gd2), 1)
  expect_error(read_gff3(f, feature_types = "gene"), "no features")
})

test_that("read_fasta_lengths computes lengths and parses IDs", {
  f <- tempfile(fileext = ".fna")
  writeLines(c(">seqA some description", "ACGTACGTACGTACGT", "ACGTACGT",
               ">seqB", "TTTTGGGG"), f)
  fl <- read_fasta_lengths(f)
  expect_equal(fl$seq_id, c("seqA", "seqB"))
  expect_equal(fl$length, c(24, 8))

  f2 <- tempfile(fileext = ".fna")
  writeLines(c(">NC_000001.1|cds|geneA", "ACGTACGT",
               ">NC_000002.1|cds|geneB", "TTTT"), f2)
  fl2 <- read_fasta_lengths(f2, header_delim = "|")
  expect_equal(fl2$seq_id, c("NC_000001.1", "NC_000002.1"))
})

test_that("read_fasta_lengths reports empty files and missing headers", {
  f <- tempfile()
  writeLines(character(0), f)
  expect_error(read_fasta_lengths(f), "empty")
  writeLines(c("ACGT", "TGCA"), f)
  expect_error(read_fasta_lengths(f), "no FASTA headers")
})

test_that("read_* functions accept `files` and combine multiple files", {
  # FASTA: one sequence per file
  f1 <- tempfile(fileext = ".fna")
  writeLines(c(">seqA", "ACGT"), f1)
  f2 <- tempfile(fileext = ".fna")
  writeLines(c(">seqB", "ACGTACGT"), f2)
  fl <- read_fasta_lengths(files = c(f1, f2))
  expect_equal(fl$seq_id, c("seqA", "seqB"))
  expect_equal(fl$length, c(4, 8))

  # GFF3: combine features from two files
  g1 <- tempfile(fileext = ".gff3")
  writeLines(paste("seqA", "s", "CDS", "1", "100", ".", "+", "0",
                   "ID=a;product=alpha", sep = "\t"), g1)
  g2 <- tempfile(fileext = ".gff3")
  writeLines(paste("seqA", "s", "CDS", "200", "300", ".", "-", "0",
                   "ID=b;product=beta", sep = "\t"), g2)
  gd <- read_gff3(files = c(g1, g2))
  expect_equal(nrow(gd), 2)
  expect_equal(gd$anno, c("alpha", "beta"))

  # BLAST: combine files with different column layouts
  b1 <- tempfile(fileext = ".o6")
  writeLines(paste("seqA", "seqB", "98", "10", "0", "0",
                   "1", "10", "1", "10", "1e-5", "100", sep = "\t"), b1)
  b2 <- tempfile(fileext = ".o7")
  writeLines(paste("seqA", "seqC", "99", "20", "0", "0",
                   "1", "20", "1", "20", "1e-8", "200", "80", "5000",
                   "5000", "plus", "seqC", sep = "\t"), b2)
  rb <- read_blast(files = c(b1, b2))
  expect_equal(nrow(rb), 2)
  expect_true(all(c("qaccver", "qcovs", "stitle") %in% names(rb)))
})

test_that("read_* `files` supports wildcards and reports no-match clearly", {
  d <- tempfile(pattern = "ggchord_fasta_")
  dir.create(d)
  writeLines(c(">seqA", "ACGT"), file.path(d, "a.fna"))
  writeLines(c(">seqB", "ACGTACGT"), file.path(d, "b.fna"))
  fl <- read_fasta_lengths(files = file.path(d, "*.fna"))
  expect_equal(fl$seq_id, c("seqA", "seqB"))
  expect_error(read_fasta_lengths(files = file.path(d, "*.gff3")),
               "no files match pattern")
})
