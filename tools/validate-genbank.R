# Explicit local check. Select an existing Python environment with Biopython.
# Example: GGCHORD_PYTHON=/path/to/python Rscript tools/validate-genbank.R
pkgload::load_all(quiet = TRUE)
local({
  path <- tempfile(fileext = ".gbff")
  on.exit(unlink(path), add = TRUE)
  header <- function(id, size, topology = "linear") c(
    sprintf("LOCUS       %-16s %11d bp    DNA     %-8s BCT 01-JAN-2000", id, size, topology),
    paste("DEFINITION  Synthetic local adapter validation", id, "."),
    paste("ACCESSION  ", id), paste("VERSION    ", paste0(id, ".2")),
    "KEYWORDS    .", "SOURCE      synthetic construct", "  ORGANISM  synthetic construct",
    "            other sequences; artificial sequences.",
    "FEATURES             Location/Qualifiers")
  first <- c(header("recordA", 60, "circular"),
    "     source          1..60", "                     /organism=\"synthetic construct\"",
    "     CDS             complement(join(1..9,31..39))",
    "                     /locus_tag=\"geneA\"", "                     /protein_id=\"proteinA.1\"",
    "                     /codon_start=1", "                     /transl_table=11",
    "                     /product=\"joined coding feature\"",
    "     misc_feature    <45..>55", "                     /note=\"retain fuzzy bounds\"",
    "ORIGIN", paste("        1", paste(rep("acgtacgtac", 6), collapse = " ")), "//")
  # No ORIGIN bases: retain the declared full length, not observed endpoints.
  second <- c(header("recordB", 90),
    "     CDS             10..18", "                     /protein_id=\"proteinB.1\"",
    "                     /codon_start=1", "                     /transl_table=11", "ORIGIN", "//")
  writeLines(c(first, second), path)
  result <- read_genbank(path, python = Sys.getenv("GGCHORD_PYTHON", "python3"), source_file = TRUE)
  stopifnot(identical(result$seq_data$accver, c("recordA.2", "recordB.2")),
    identical(as.numeric(result$seq_data$length), c(60, 90)),
    result$seq_data$circular[1L], !result$seq_data$sequence_defined[2L])
  coding <- result$features$annotation[result$features$annotation$type == "CDS", ]
  stopifnot(coding$feature_id[1L] == "proteinA.1", coding$strand[1L] == "-",
    identical(as.numeric(coding$segments[[1L]]$start), c(31, 1)),
    identical(as.numeric(coding$segments[[1L]]$end), c(39, 9)),
    coding$genetic_code[1L] == 11, coding$codon_start[1L] == 1,
    coding$qualifiers[[1L]]$locus_tag[[1L]] == "geneA",
    result$features$annotation$partial_start[3L],
    result$features$annotation$partial_end[3L],
    length(result$report$raw_files[[1L]]) == length(c(first, second)))
  prepared <- prepare_ggchord_data(result$seq_data, features = result$features)
  stopifnot(prepared$report$valid,
    prepared$report$conversions$features$annotation$source_import$format == "GenBank")
  # Preserve a declared origin join, and distinguish incomplete codons from
  # uncertain genomic boundaries rather than calling exact bounds complete.
  joined <- c(header("origin", 60, "circular"),
    "     CDS             join(51..60,1..10)",
    "                     /protein_id=\"origin_protein\"",
    "                     /codon_start=2", "ORIGIN", "//")
  writeLines(joined, path)
  joined <- read_genbank(path, python = Sys.getenv("GGCHORD_PYTHON", "python3"), feature_types = "CDS")
  coding <- joined$features$annotation
  model <- build_ggchord_cds(coding, partial = "trim")
  domain <- project_ggchord_features(data.frame(feature_id = "first", parent_id = "origin_protein",
    start = 1, end = 1), model)
  stopifnot(identical(coding$segments[[1L]]$segment_order, 1:2),
    all(coding$segments[[1L]]$start == c(51, 1)),
    coding$partial5, coding$partial3, !coding$bound_partial5,
    model$models$codon_count == 6, domain$start == 52, domain$end == 54)
  default <- c(header("default", 60), "     CDS             10..18",
    "                     /protein_id=\"default_protein\"", "ORIGIN", "//")
  writeLines(default, path)
  default <- read_genbank(path, python = Sys.getenv("GGCHORD_PYTHON", "python3"), feature_types = "CDS")
  stopifnot(default$features$annotation$codon_start == 1L,
    default$features$annotation$codon_start_source == "format-default",
    build_ggchord_cds(default$features$annotation)$models$codon_count == 3)
  # Missing/remote topology is retained as an explicit unplaced source record.
  remote <- c(header("recordC", 90), "     misc_feature    REMOTE.1:10..20", "ORIGIN", "//")
  writeLines(remote, path)
  reported <- read_genbank(path, python = Sys.getenv("GGCHORD_PYTHON", "python3"),
    unsupported = "report")
  stopifnot(length(reported$report$unplaced_features) == 1L,
    nrow(reported$features$annotation) == 0L,
    grepl("remote", reported$report$unplaced_features[[1L]]$reason))
})
if (requireNamespace("Biostrings", quietly = TRUE)) {
  sequences <- Biostrings::DNAStringSet(c("A.2" = "ACGT", "B" = "AT"))
  converted <- as_ggchord_seq(sequences)
  single <- as_ggchord_seq(Biostrings::DNAString("ACGT"), ids = "explicit")
  stopifnot(identical(converted$accver, c("A.2", "B")),
    identical(converted$length, c(4, 2)), single$accver == "explicit", single$length == 4)
} else cat("Biostrings unavailable: real S4-object checks remain a candidate-release validation item.\n")
cat("Real GenBank parser, joined locations and provenance checks passed.\n")
