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
