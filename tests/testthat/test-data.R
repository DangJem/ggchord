test_that("validation and cleaning preserve usable input", {
  seqs <- data.frame(accver = c("A", "B"), length = 1000)
  links <- as_ggchord_link(data.frame(qaccver = "A", saccver = "B",
    qstart = 100, qend = 200, sstart = 400, send = 300))
  data <- ggchord_data(seqs, links = list(hsps = links))
  expect_true(validate_ggchord_data(data)$valid)
  expect_equal(clean_ggchord_data(data)$links$hsps$sstrand, "-")
  expect_error(ggchord(transform(seqs, length = "bad")), "length")
})

test_that("basic readers retain identifiers and coordinate conventions", {
  files <- setNames(vapply(c("fa", "tsv", "gff", "bed", "paf", "bedpe", "bedGraph"),
    function(ext) tempfile(fileext = paste0(".", ext)), character(1)),
    c("fa", "tsv", "gff", "bed", "paf", "bedpe", "bedGraph"))
  on.exit(unlink(files), add = TRUE)
  writeLines(c(">A", "ACgtACGT"), files["fa"])
  writeLines(paste("001", "NA", 98, 100, 0, 0, 1, 100, 1, 100,
    "1e-20", 200, sep = "\t"), files["tsv"])
  writeLines(c("##gff-version 3",
    "A\tsource\tCDS\t1\t8\t.\t+\t0\tID=cds1;product=test%20protein"), files["gff"])
  writeLines("A\t0\t1\tpoint\t0\t.", files["bed"])
  writeLines("A\t1000\t99\t200\t-\tB\t1000\t299\t400\t90\t101\t255", files["paf"])
  writeLines("browser\t0\t1\tB\t2\t4\tpair\t.\t+\t-", files["bedpe"])
  writeLines("track\t0\t1\t-0.5", files["bedGraph"])
  expect_equal(read_fasta_lengths(files["fa"], sequence_unit = "bp")$length, 8L)
  hit <- read_blast(files["tsv"], format = "outfmt6")
  expect_equal(c(hit$qaccver, hit$saccver), c("001", "NA"))
  expect_equal(read_gff3(files["gff"])$label, "test protein")
  expect_equal(read_bed(files["bed"])$start, 1)
  expect_equal(read_paf(files["paf"])$sstart, 300)
  expect_equal(read_bedpe(files["bedpe"])$sstart, 3)
  expect_equal(read_bedgraph(files["bedGraph"])$value, -.5)
})

test_that("link utilities each retain a basic success path", {
  link <- data.frame(qaccver = "A", saccver = "B", qstart = 100, qend = 200,
    sstart = 300, send = 400, qstrand = "+", sstrand = "-", pident = 95, length = 101)
  expect_equal(nrow(filter_ggchord_links(link, min_pident = 90)), 1L)
  expect_equal(nrow(deduplicate_ggchord_links(rbind(link, link))), 1L)
  blocks <- rbind(link, transform(link, qstart = 201, qend = 301,
    sstart = 199, send = 299))
  expect_equal(nrow(merge_ggchord_links(blocks)), 1L)
})
