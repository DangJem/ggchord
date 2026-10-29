# Explicit local development check; does not run testthat or package checks.
pkgload::load_all(quiet = TRUE)
local({
  file <- tempfile(fileext = ".gff3")
  on.exit(unlink(file), add = TRUE)
  writeLines(c("##gff-version 3",
    "A\tsrc\tgene\t100\t204\t.\t+\t.\tID=gene",
    "A\tsrc\tCDS\t100\t103\t.\t+\t0\tID=coding",
    "A\tsrc\tCDS\t200\t204\t.\t+\t2\tID=coding"), file)
  model <- build_ggchord_cds(read_gff3(file))
  stopifnot(model$models$codon_count == 3)
})
cds_data <- data.frame(feature_id = c("plus", "plus", "minus", "minus"),
  accver = c("A", "A", "B", "B"), start = c(100, 200, 400, 300),
  end = c(103, 204, 403, 304), strand = c("+", "+", "-", "-"), phase = c(0, 2, 0, 2))
cds <- build_ggchord_cds(cds_data, genetic_code = c(plus = 11, minus = 1))
stopifnot(identical(cds$models$codon_count, c(3, 3)),
  identical(cds$models$genetic_code, c(11, 1)))
domains <- data.frame(feature_id = c("domainA", "domainB"), parent_id = c("plus", "minus"),
  start = 2, end = 3, label = c("A-domain", "B-domain"))
original <- domains
projected <- project_ggchord_features(domains, cds)
stopifnot(identical(domains, original), identical(projected$strand, c("+", "-")),
  identical(projected$segments[[1]]$start, c(103, 200)),
  identical(projected$segments[[1]]$end, c(103, 204)),
  identical(projected$segments[[2]]$start, c(400, 300)),
  identical(projected$segments[[2]]$end, c(400, 304)),
  identical(projected$parent_ids[[1]], "plus"), projected$label[1] == "A-domain")
roundtrip <- project_ggchord_features(projected, cds, from = "genome")
stopifnot(identical(roundtrip$start, c(2, 2)), identical(roundtrip$end, c(3, 3)),
  all(roundtrip$excluded_bases == 0), all(!roundtrip$codon_expanded))
# A whole envelope explicitly excludes intronic sequence; coding boundaries remain exact.
envelope <- projected; envelope$segments <- NULL
inverse <- project_ggchord_features(envelope, cds, from = "genome")
stopifnot(all(inverse$excluded_bases > 0), identical(inverse$start, c(2, 2)))
links <- data.frame(qfeature_id = "plus", sfeature_id = "minus", qstart = 1,
  qend = 3, sstart = 1, send = 3, pident = 80, length = 3)
paired <- project_ggchord_links(links, cds, ungapped = TRUE)
stopifnot(nrow(paired) == 2L, sum(paired$projected_bases) == 9,
  all(paired$qstart <= paired$qend), all(paired$sstart >= paired$send),
  all(paired$direction == "reverse"), all(paired$source_pident == 80),
  !any(c("pident", "length") %in% names(paired)))
# Partial CDS starts at the first complete observed codon, without inventing an upstream residue.
partial <- data.frame(feature_id = "p", accver = "C", start = c(10, 30),
  end = c(13, 36), strand = "+", phase = c(1, 0), partial5 = TRUE)
trimmed <- build_ggchord_cds(partial, partial = "trim")
stopifnot(trimmed$models$codon_count == 3, trimmed$models$terminal_bases == 1,
  trimmed$models$partial3, trimmed$segments$p$start[1] == 11)
# Reject bad phase, declared exceptions, missing parent and partial-codon inverse.
error <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
stopifnot(error(build_ggchord_cds(transform(cds_data, phase = 0))),
  error(build_ggchord_cds(partial)),
  error(build_ggchord_cds(transform(cds_data, frameshift = TRUE))),
  error(project_ggchord_features(transform(domains, parent_id = "missing"), cds)),
  error(project_ggchord_links(links, cds)))
split_codon <- data.frame(feature_id = "split", parent_id = "plus", accver = "A", start = 101, end = 102)
stopifnot(error(project_ggchord_features(split_codon, cds, from = "genome")))
expanded <- project_ggchord_features(split_codon, cds, from = "genome", codon_boundary = "expand")
stopifnot(expanded$start == 1, expanded$end == 1, expanded$codon_expanded)
# Nested source pieces use explicit transcription order and derive phase from codon_start.
nested <- data.frame(feature_id = "nested", accver = "N", start = 10, end = 103,
  strand = "-", codon_start = 1, genetic_code = 11)
nested$segments <- I(list(data.frame(start = c(100, 10), end = c(103, 14),
  strand = "-", segment_index = 1:2)))
nested_map <- build_ggchord_cds(nested)
stopifnot(identical(nested_map$segments$nested$phase, c(0, 2)),
  identical(nested_map$segments$nested$start, c(100, 10)),
  all(nested_map$segments$nested$source_row == 1L))
# A generic BED-style segment_index remains genomic indexing, not transcript order.
bed_order <- nested
bed_order$segments <- I(list(nested$segments[[1L]][2:1, ]))
bed_order$segments[[1L]]$segment_index <- 1:2
bed_map <- build_ggchord_cds(bed_order)
stopifnot(identical(bed_map$segments$nested$start, c(100, 10)),
  identical(bed_map$segments$nested$phase, c(0, 2)))
# A declared circular join retains its stated order rather than genomic sorting.
circular <- nested
circular$strand <- "+"
circular$segments[[1L]]$strand <- "+"
circular$segments[[1L]]$segment_order <- 1:2
circular_map <- build_ggchord_cds(circular)
stopifnot(identical(circular_map$segments$nested$start, c(100, 10)),
  identical(circular_map$segments$nested$phase, c(0, 2)))
# Disjoint domains remain disjoint across both projection directions.
simple <- build_ggchord_cds(data.frame(feature_id = "simple", accver = "S",
  start = 101, end = 109, strand = "+", phase = 0))
disjoint <- data.frame(feature_id = "d", parent_id = "simple", accver = "S", start = 101, end = 109)
disjoint$segments <- I(list(data.frame(start = c(101, 107), end = c(103, 109))))
aa <- project_ggchord_features(disjoint, simple, from = "genome")
back <- project_ggchord_features(aa, simple)
stopifnot(identical(back$segments[[1L]]$start, c(101, 107)),
  identical(back$segments[[1L]]$end, c(103, 109)))
# Projected domain segments work directly with the existing biological feature layer.
plot <- ggchord(data.frame(accver = c("A", "B"), length = 1000)) + geom_seq() +
  geom_feature(data = projected) + coord_collinear()
geometry <- export_ggchord_layout(plot, include = "feature")
stopifnot(setequal(geometry$feature$source_row, 1:2))
cat("CDS phase, bidirectional projection, parent domains and protein links passed.\n")
# Reuse an explicitly selected optional parser environment for one joined CDS.
if (nzchar(Sys.getenv("GGCHORD_PYTHON"))) local({
  file <- tempfile(fileext = ".gbff")
  on.exit(unlink(file), add = TRUE)
  writeLines(c("LOCUS       project                   40 bp    DNA     linear   BCT 01-JAN-2000",
    "DEFINITION  Synthetic projection fixture.", "ACCESSION   project", "VERSION     project.1",
    "KEYWORDS    .", "SOURCE      synthetic construct", "  ORGANISM  synthetic construct",
    "            other sequences; artificial sequences.", "FEATURES             Location/Qualifiers",
    "     CDS             complement(join(1..5,21..24))",
    "                     /protein_id=\"coding.1\"", "                     /codon_start=1",
    "                     /transl_table=11", "ORIGIN",
    "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt", "//"), file)
  imported <- read_genbank(file, python = Sys.getenv("GGCHORD_PYTHON"), feature_types = "CDS")
  model <- build_ggchord_cds(imported$features$annotation)
  mapped <- project_ggchord_features(data.frame(feature_id = "domain", parent_id = "coding.1",
    start = 2, end = 3), model)
  stopifnot(identical(model$segments$coding.1$phase, c(0, 2)), model$models$genetic_code == 11,
    mapped$strand == "-", identical(mapped$segments[[1L]]$start, c(21, 1)),
    identical(mapped$segments[[1L]]$end, c(21, 5)))
  cat("Real GenBank joined-CDS to domain projection passed.\n")
})
