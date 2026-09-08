# Build the internal normalized common-feature database from the audited
# workbook. readxl and digest are development-only dependencies.

workbook <- "examples/standardCommonFeatures_tables.xlsx"
expected_sha256 <- "3d9e7a78c705de2beeb9d8af3c74dede714faa468cd6eee61f78e0f784b65312"
if (!requireNamespace("readxl", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE)) {
  stop("Install readxl and digest to regenerate the common-feature database",
    call. = FALSE)
}
actual_sha256 <- digest::digest(workbook, algo = "sha256", file = TRUE,
  serialize = FALSE)
stopifnot(identical(actual_sha256, expected_sha256))

read_sheet <- function(name) {
  as.data.frame(readxl::read_excel(workbook, sheet = name),
    stringsAsFactors = FALSE, check.names = FALSE)
}
metadata_table <- read_sheet("FileMetadata")
features <- read_sheet("Features")
segments <- read_sheet("Segments")
qualifiers <- read_sheet("Qualifiers")
qualifier_links <- read_sheet("QualifierLinks")
feature_type_summary <- read_sheet("FeatureTypes")

required <- list(
  Features = c("feature_id", "name", "type", "directionality_label",
    "detectionMode", "geneticCode", "segment_count",
    "reference_dna_feature_5to3", "reference_protein"),
  Segments = c("feature_id", "segment_index", "segment_type", "start", "end",
    "length_bp", "color", "translated", "segment_name",
    "dna_sequence_top_strand"),
  Qualifiers = c("feature_id", "qualifier_name", "value_index",
    "value_display"),
  QualifierLinks = c("feature_id", "qualifier_name", "value_index", "url"),
  FeatureTypes = c("feature_type", "feature_count")
)
tables <- list(Features = features, Segments = segments,
  Qualifiers = qualifiers, QualifierLinks = qualifier_links,
  FeatureTypes = feature_type_summary)
for (name in names(required)) {
  missing <- setdiff(required[[name]], names(tables[[name]]))
  if (length(missing)) stop(name, " is missing columns: ",
    paste(missing, collapse = ", "), call. = FALSE)
}
stopifnot(
  nrow(features) == 1273L,
  nrow(segments) == 1730L,
  nrow(qualifiers) == 4551L,
  !anyDuplicated(features$feature_id),
  all(segments$feature_id %in% features$feature_id),
  all(qualifiers$feature_id %in% features$feature_id),
  all(qualifier_links$feature_id %in% features$feature_id),
  all(segments$segment_type %in% c("standard", "gap")),
  all(segments$start >= 1L),
  all(segments$end >= segments$start),
  all(segments$length_bp == segments$end - segments$start + 1L)
)
standard_dna <- toupper(gsub("[[:space:]]", "",
  segments$dna_sequence_top_strand[segments$segment_type == "standard"]))
stopifnot(all(is.na(standard_dna) | grepl("^[ACGTRYSWKMBDHVN]+$", standard_dna)))

features$common_feature_id <- as.character(features$feature_id)
source_metadata <- stats::setNames(metadata_table$value, metadata_table$key)
ggchord_common_feature_database <- list(
  metadata = list(
    database_id = "ggchord-common-features",
    database_version = paste0("export-", source_metadata[["exportVersion"]]),
    source = as.character(source_metadata[["Source file"]]),
    source_workbook_sha256 = actual_sha256,
    source_file_sha256 = as.character(source_metadata[["SHA256"]]),
    export_version = as.integer(source_metadata[["exportVersion"]]),
    import_version = as.integer(source_metadata[["importVersion"]]),
    record_counts = c(features = nrow(features), segments = nrow(segments),
      qualifiers = nrow(qualifiers), qualifier_links = nrow(qualifier_links)),
    generating_version = "0.13.0"
  ),
  features = features,
  segments = segments,
  qualifiers = qualifiers,
  qualifier_links = qualifier_links,
  feature_type_summary = feature_type_summary,
  search_indexes = list(
    dna_feature_ids = sort(unique(segments$feature_id[
      segments$segment_type == "standard" &
        !is.na(segments$dna_sequence_top_strand) &
        nzchar(segments$dna_sequence_top_strand)])),
    protein_feature_ids = sort(features$feature_id[
      !is.na(features$reference_protein) & nzchar(features$reference_protein)]),
    feature_name = split(features$feature_id,
      factor(tolower(features$name), levels = unique(tolower(features$name)))),
    feature_type = split(features$feature_id,
      factor(features$type, levels = unique(features$type)))
  )
)

old <- new.env(parent = emptyenv())
if (file.exists("R/sysdata.rda")) load("R/sysdata.rda", envir = old)
objects <- list(ggchord_common_feature_database = ggchord_common_feature_database)
if (exists("ggchord_rebase_database", envir = old, inherits = FALSE)) {
  objects$ggchord_rebase_database <- get("ggchord_rebase_database", envir = old)
}
list2env(objects, envir = environment())
save(list = names(objects), file = "R/sysdata.rda", compress = "xz")
