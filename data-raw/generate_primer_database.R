# Generate the internal sequencing-primer catalogue from the reviewed workbook.
# Keep this script deterministic so a later workbook release can be audited and
# rebuilt by updating the expected digest and row-count assertions together.

source_path <- "examples/Addgene_Sequencing_Primers.xlsx"
source_url <- "https://www.addgene.org/mol-bio-reference/sequencing-primers/"
expected_sha256 <-
  "11ea7c65ba07d3d4ab450900ece6ebc254e83dca7f54037fe6ad2f92b5ff5b09"

stopifnot(
  file.exists(source_path),
  requireNamespace("readxl", quietly = TRUE),
  requireNamespace("digest", quietly = TRUE)
)
actual_sha256 <- digest::digest(source_path, algo = "sha256", file = TRUE)
stopifnot(identical(actual_sha256, expected_sha256))

read_primer_sheet <- function(sheet, is_universal) {
  out <- as.data.frame(
    readxl::read_excel(source_path, sheet = sheet),
    stringsAsFactors = FALSE
  )
  expected <- c("Name", "Sequence (5' to 3')", "Description", "Direction")
  stopifnot(identical(names(out), expected))
  names(out) <- c("name", "sequence", "description", "direction_hint")
  out[] <- lapply(out, function(x) trimws(gsub("\\u00a0", " ", as.character(x))))
  out$sequence <- toupper(gsub("[[:space:]]", "", out$sequence))
  out$is_universal <- is_universal
  out$source_sheet <- sheet
  out$source_row <- seq_len(nrow(out)) + 1L
  out
}

universal <- read_primer_sheet("Universal Sequencing Primers", TRUE)
common <- read_primer_sheet("Other Common Sequencing Primers", FALSE)
stopifnot(nrow(universal) == 23L, nrow(common) == 140L)
aliases <- rbind(universal, common)
stopifnot(
  !anyNA(aliases),
  all(nzchar(aliases$name)),
  all(grepl("^[ACGTRYSWKMBDHVN]+$", aliases$sequence)),
  all(aliases$direction_hint %in% c("Forward", "Reverse"))
)

primer_sequences <- sort(unique(aliases$sequence))
primer_ids <- paste0(
  "primer_",
  vapply(primer_sequences, digest::digest, character(1L), algo = "sha256",
    serialize = FALSE),
  USE.NAMES = FALSE
)
primer_ids <- substr(primer_ids, 1L, 23L)
id_map <- stats::setNames(primer_ids, primer_sequences)
aliases$primer_id <- unname(id_map[aliases$sequence])

# Preserve every source row as an alias record. A repeated oligo may have
# biologically distinct names; the matching engine searches it once and keeps
# all aliases instead of selecting a misleading global canonical name.
aliases <- aliases[, c(
  "primer_id", "name", "sequence", "description", "direction_hint",
  "is_universal", "source_sheet", "source_row"
)]
aliases <- aliases[order(aliases$primer_id, !aliases$is_universal,
  aliases$source_sheet, aliases$source_row), , drop = FALSE]
rownames(aliases) <- NULL

preferred_name <- vapply(primer_sequences, function(sequence) {
  rows <- aliases[aliases$sequence == sequence, , drop = FALSE]
  rows$name[order(!rows$is_universal, nchar(rows$name), rows$name)][1L]
}, character(1L))
is_universal <- vapply(primer_sequences, function(sequence) {
  any(aliases$is_universal[aliases$sequence == sequence])
}, logical(1L))

catalog <- data.frame(
  primer_id = primer_ids,
  preferred_name = preferred_name,
  sequence = primer_sequences,
  length = nchar(primer_sequences),
  is_universal = is_universal,
  source = "Addgene sequencing primer reference",
  stringsAsFactors = FALSE
)

ggchord_primer_database <- list(
  metadata = list(
    database_id = "ggchord-sequencing-primers",
    database_version = "2025-10-22",
    source = "Addgene Sequencing Primers",
    source_url = source_url,
    content_last_reviewed = "2025-10-22",
    source_workbook = basename(source_path),
    source_workbook_sha256 = actual_sha256,
    record_counts = c(catalog = nrow(catalog), aliases = nrow(aliases)),
    generating_version = "0.13.0"
  ),
  catalog = catalog,
  aliases = aliases
)

stopifnot(
  nrow(catalog) == 137L,
  nrow(aliases) == 163L,
  nrow(unique(aliases[c("name", "sequence", "description", "direction_hint")])) == 145L,
  !anyDuplicated(catalog$primer_id),
  all(aliases$primer_id %in% catalog$primer_id)
)

old <- new.env(parent = emptyenv())
if (file.exists("R/sysdata.rda")) load("R/sysdata.rda", envir = old)
objects <- list(ggchord_primer_database = ggchord_primer_database)
for (name in c("ggchord_common_feature_database", "ggchord_rebase_database")) {
  if (exists(name, envir = old, inherits = FALSE)) {
    objects[[name]] <- get(name, envir = old, inherits = FALSE)
  }
}
list2env(objects, envir = environment())
save(list = names(objects), file = "R/sysdata.rda", compress = "xz")
