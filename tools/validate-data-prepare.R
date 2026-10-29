# Explicit development checks, not a package test/check run.
pkgload::load_all(quiet = TRUE)
seqs <- data.frame(accver = c("A", "B"), length = c(100, 100),
  genome = c("isolate1", "isolate2"))
rb <- data.frame(qaccver = "A", saccver = "B", qstart = c(1, 1, 11),
  qend = c(10, 10, 20), sstart = c(20, 20, 30), send = c(29, 29, 39),
  length = 10, pident = 95)
sites <- data.frame(accver = "A", start = 30, end = 40, feature_id = "site1")
original <- rb
out <- prepare_ggchord_data(seqs, rb, features = list(sites = sites),
  steps = list(deduplicate = list(keep = "first"), merge = list()))
stopifnot(out$report$valid, identical(rb, original), nrow(out$ribbon_data) == 1L,
  identical(out$report$source_rows[[1L]], c(1L, 3L)),
  out$features$sites$strand == ".", out$seq_data$genome[1L] == "isolate1")

unchanged <- prepare_ggchord_data(seqs, rb)
stopifnot(nrow(unchanged$ribbon_data) == nrow(rb),
  identical(unchanged$ribbon_data$qstart, rb$qstart),
  length(unchanged$report$steps) == 0L)
invalid <- prepare_ggchord_data(seqs, transform(rb, qend = 101),
  features = list(sites = transform(sites, accver = "missing")))
stopifnot(!invalid$report$valid, nrow(invalid$ribbon_data) == nrow(rb),
  "unknown_id" %in% invalid$report$feature_issues$category)

directional <- rb[1:2, ]
directional$sstart <- c(20, 29)
directional$send <- c(29, 20)
kept <- prepare_ggchord_data(seqs, directional,
  steps = list(clean = list(), deduplicate = list(keep = "first"),
    bundle = list(bins = 1, weight = "count")))
stopifnot(nrow(kept$ribbon_data) == 2L,
  all(c("same", "reverse") %in% kept$ribbon_data$direction))
rejected <- tryCatch(prepare_ggchord_data(seqs, rb,
  steps = list(clean = list(reversed_interval = "sort"), bundle = list())), error = identity)
stopifnot(inherits(rejected, "error"), grepl("endpoint direction", conditionMessage(rejected)))
abbreviated <- tryCatch(prepare_ggchord_data(seqs, rb,
  steps = list(clean = list(reversed_interval = "s"), bundle = list())), error = identity)
stopifnot(inherits(abbreviated, "error"), grepl("endpoint direction", conditionMessage(abbreviated)))
sorted <- clean_ggchord_data(seqs, directional)$ribbon_data
conflict <- tryCatch(prepare_ggchord_data(seqs, sorted,
  steps = list(bundle = list(weight = "count"))), error = identity)
stopifnot(inherits(conflict, "error"), grepl("direction contradicts", conditionMessage(conflict)))
conflict_clean <- tryCatch(prepare_ggchord_data(seqs, sorted,
  steps = list(clean = list(), bundle = list(weight = "count"))), error = identity)
stopifnot(inherits(conflict_clean, "error"))

# GFF3: split multiple parents before decoding; retain repeated IDs and source strand.
gff <- tempfile(fileext = ".gff3")
writeLines(c("##gff-version 3",
  "A\tsrc\tCDS\t1\t10\t.\t.\t0\tID=part;Parent=one%2Ctwo,three",
  "A\tsrc\tCDS\t20\t30\t.\t-\t0\tID=part;Parent=three"), gff)
imported <- read_gff3(gff, unstranded = "keep")
stopifnot(identical(imported$feature_id, c("part", "part")),
  identical(imported$parent_ids[[1L]], c("one,two", "three")),
  imported$strand[1L] == ".", imported$source_strand[1L] == ".")
converted <- prepare_ggchord_data(seqs, features = list(annotation = imported))
stopifnot(identical(converted$report$conversions$features$annotation$source_import$files, gff),
  converted$report$conversions$features$annotation$source_import$format == "GFF3")
unlink(gff)
plot <- ggchord(out$seq_data, out$ribbon_data) + geom_seq() +
  geom_link_ribbon() + geom_feature(data = out$features$sites) + coord_collinear()
built <- ggplot2::ggplot_build(plot)
stopifnot(inherits(built, "ggplot_built"), nrow(built$data[[3L]]) > 0L)
cat("Explicit preparation, provenance and GFF3 ID checks passed.\n")
