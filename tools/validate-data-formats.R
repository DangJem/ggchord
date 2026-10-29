# Explicit lightweight format checks; no testthat, plotting or package check.
pkgload::load_all(quiet = TRUE)
local({
  directory <- tempfile("ggchord-formats-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  bed <- file.path(directory, "blocks.bed")
  writeLines(c("track name=example", paste("contig.2", 10, 40, "gene", 500,
    "-", 10, 40, "0", 2, "10,10,", "0,20,", sep = "\t")), bed)
  blocks <- read_bed(bed, source_file = TRUE)
  stopifnot(blocks$accver == "contig.2", blocks$start == 11, blocks$end == 40,
    blocks$strand == "-", blocks$chromStart == 10, blocks$blockStarts == "0,20,",
    identical(blocks$segments[[1L]]$start, c(11, 31)),
    identical(blocks$segments[[1L]]$end, c(20, 40)),
    blocks$.source_line == 2L, blocks$.source_file == bed)
  one_bp <- file.path(directory, "one.bed")
  writeLines("contig.2\t0\t1", one_bp)
  intervals <- read_bed(files = c(bed, one_bp))
  stopifnot(nrow(intervals) == 2L, intervals$start[2L] == 1,
    intervals$end[2L] == 1, intervals$strand[2L] == ".",
    intervals$directionality[2L] == "none", is.list(intervals$segments),
    attr(intervals, "ggchord_import")$input_coordinates == "0-based-half-open")

  paf <- file.path(directory, "links.paf")
  writeLines(c("q\t100\t0\t90\t-\ts\t200\t20\t110\t85\t90\t255\tcg:Z:90M\txx:Z:keep",
    "q\t100\t10\t20\t+\ts\t200\t0\t10\t10\t10\t60"), paf)
  links <- read_paf(paf)
  stopifnot(identical(links$qstart, c(1, 11)), identical(links$qend, c(90, 20)),
    identical(links$sstart, c(110, 1)), identical(links$send, c(21, 10)),
    links$direction[1L] == "reverse", is.na(links$mapq[1L]),
    links$source_mapq[1L] == 255, links$mapq[2L] == 60,
    links$length[1L] == 90, links$nmatch[1L] == 85,
    links$pident[2L] == 100,
    identical(links$tags[[1L]], c("cg:Z:90M", "xx:Z:keep")),
    length(links$tags[[2L]]) == 0L)
  # Two invalid topology cases: empty interval and overlapping BED12 blocks.
  invalid <- file.path(directory, "invalid.bed")
  for (line in c("q\t5\t5", "q\t0\t10\tx\t0\t+\t0\t10\t0\t2\t8,5\t0,5")) {
    writeLines(line, invalid)
    stopifnot(inherits(tryCatch(read_bed(invalid), error = identity), "error"))
  }
})
cat("BED/PAF coordinate and metadata checks passed.\n")
