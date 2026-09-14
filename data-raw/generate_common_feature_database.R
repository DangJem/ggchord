# Build the internal normalized common-feature database from the audited
# workbook and the audited annotations embedded in the example SnapGene files.
# readxl, digest, and xml2 are development-only dependencies.

workbook <- "examples/standardCommonFeatures_8.2.3_tables.xlsx"
legacy_provenance_path <- "data-raw/common-feature-legacy-provenance.rds"
expected_sha256 <- "7d9294df1ba89e6d9437f54d2aed5681d86ff4663f098fcc3d866b9ec58db69b"
if (!requireNamespace("readxl", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE) ||
    !requireNamespace("xml2", quietly = TRUE)) {
  stop("Install readxl, digest, and xml2 to regenerate the common-feature database",
    call. = FALSE)
}
actual_sha256 <- digest::digest(workbook, algo = "sha256", file = TRUE,
  serialize = FALSE)
stopifnot(identical(actual_sha256, expected_sha256))
legacy_provenance <- readRDS(legacy_provenance_path)
stopifnot(
  identical(legacy_provenance$schema_version, 1L),
  nrow(legacy_provenance$legacy_release) == 1L,
  nrow(legacy_provenance$legacy_feature_records) == 1273L,
  nrow(legacy_provenance$feature_correspondence) == 1454L
)

read_sheet <- function(name, path = workbook) {
  as.data.frame(readxl::read_excel(path, sheet = name, guess_max = 100000L),
    stringsAsFactors = FALSE, check.names = FALSE)
}
metadata_table <- read_sheet("README")
raw_features <- read_sheet("Features")
raw_segments <- read_sheet("Segments")
raw_qualifiers <- read_sheet("Qualifiers")
raw_qualifier_links <- read_sheet("QualifierLinks")
feature_type_summary <- read_sheet("FeatureTypes")
raw_source_sequence <- read_sheet("SourceSequence")

required <- list(
  Features = c("feature_id", "recentID", "name", "type",
    "directionality_label", "detectionMode", "geneticCode", "segment_count",
    "start_min_1based_inclusive", "end_max_1based_inclusive", "prioritize"),
  Segments = c("feature_id", "segment_index", "segment_type",
    "start_1based_inclusive", "end_1based_inclusive", "length_bp", "color",
    "translated", "segment_name", "translationNumberingStartsFrom"),
  Qualifiers = c("feature_id", "qualifier_name", "value_index", "text",
    "int", "predef", "bool", "value_display", "url_occurrence_count"),
  QualifierLinks = c("feature_id", "qualifier_name", "value_index",
    "link_index", "anchor_text", "url"),
  FeatureTypes = c("feature_type", "feature_count")
)
tables <- list(Features = raw_features, Segments = raw_segments,
  Qualifiers = raw_qualifiers, QualifierLinks = raw_qualifier_links,
  FeatureTypes = feature_type_summary)
for (name in names(required)) {
  missing <- setdiff(required[[name]], names(tables[[name]]))
  if (length(missing)) stop(name, " is missing columns: ",
    paste(missing, collapse = ", "), call. = FALSE)
}
stopifnot(
  nrow(raw_features) == 1454L,
  nrow(raw_segments) == 1970L,
  nrow(raw_qualifiers) == 5216L,
  nrow(raw_qualifier_links) == 548L,
  !anyDuplicated(raw_features$feature_id),
  all(raw_segments$feature_id %in% raw_features$feature_id),
  all(raw_qualifiers$feature_id %in% raw_features$feature_id),
  all(raw_qualifier_links$feature_id %in% raw_features$feature_id),
  all(raw_segments$segment_type %in% c("standard", "gap")),
  all(raw_segments$start_1based_inclusive >= 1L),
  all(raw_segments$end_1based_inclusive >=
    raw_segments$start_1based_inclusive),
  all(raw_segments$length_bp == raw_segments$end_1based_inclusive -
    raw_segments$start_1based_inclusive + 1L),
  sum(raw_features$detectionMode == "exactProteinMatch", na.rm = TRUE) == 250L,
  sum(raw_features$prioritize %in% TRUE, na.rm = TRUE) == 1L
)
source_metadata <- stats::setNames(metadata_table$value, metadata_table$key)
source_release_id <- "snapgene-common-features-8.2.3-export15"
legacy_source_release_id <- "snapgene-common-features-legacy-export13"
source_sequence_id <- paste0(source_release_id, ":backing-sequence")
current_sequence <- toupper(paste0(raw_source_sequence$sequence,
  collapse = ""))
stopifnot(
  nchar(current_sequence) == 903397L,
  grepl("^[ACGTRYSWKMBDHVN]+$", current_sequence),
  all(raw_source_sequence$start_1based_inclusive ==
    c(1L, utils::head(raw_source_sequence$end_1based_inclusive, -1L) + 1L)),
  all(raw_source_sequence$end_1based_inclusive -
    raw_source_sequence$start_1based_inclusive + 1L ==
    nchar(raw_source_sequence$sequence))
)

common_feature_ids <- sprintf("ggcf_%06d", seq_len(nrow(raw_features)))
id_map <- stats::setNames(common_feature_ids, raw_features$feature_id)

features <- raw_features
features$common_feature_id <- unname(id_map[features$feature_id])
features$source_release_id <- source_release_id
features$source_sequence_id <- source_sequence_id
features$source_feature_id <- as.character(features$feature_id)
features$source_recent_id <- as.integer(features$recentID)
features$feature_id <- NULL
features$recentID <- NULL
features <- features[c("common_feature_id", "source_release_id",
  "source_sequence_id", "source_feature_id", "source_recent_id",
  setdiff(names(features), c("common_feature_id", "source_release_id",
    "source_sequence_id", "source_feature_id", "source_recent_id")))]

segments <- raw_segments
segments$common_feature_id <- unname(id_map[segments$feature_id])
segments$source_sequence_id <- source_sequence_id
segments <- segments[c("common_feature_id", "source_sequence_id",
  "segment_index", "segment_type", "start_1based_inclusive",
  "end_1based_inclusive", "length_bp", "color", "translated",
  "segment_name", "translationNumberingStartsFrom")]

qualifiers <- raw_qualifiers
qualifiers$common_feature_id <- unname(id_map[qualifiers$feature_id])
qualifiers <- qualifiers[c("common_feature_id", "qualifier_name",
  "value_index", "text", "int", "predef", "bool", "value_display",
  "url_occurrence_count")]

qualifier_links <- raw_qualifier_links
qualifier_links$common_feature_id <- unname(id_map[qualifier_links$feature_id])
qualifier_links <- qualifier_links[c("common_feature_id", "qualifier_name",
  "value_index", "link_index", "anchor_text", "url")]

source_sequences <- data.frame(
  source_sequence_id = source_sequence_id,
  source_release_id = source_release_id,
  length_bp = nchar(current_sequence),
  sequence_sha256 = digest::digest(current_sequence, algo = "sha256",
    serialize = FALSE),
  sequence = current_sequence,
  stringsAsFactors = FALSE
)

current_release <- data.frame(
  source_release_id = source_release_id, snapgene_version = "8.2.3",
  source_workbook = workbook, source_workbook_sha256 = actual_sha256,
  source_file = as.character(source_metadata[["Source file"]]),
  source_file_sha256 = as.character(source_metadata[["SHA256"]]),
  source_file_bytes = as.integer(source_metadata[["file bytes"]]),
  export_version = as.integer(source_metadata[["export-format version"]]),
  import_version = as.integer(source_metadata[["import-format version"]]),
  sequence_length_bp = nchar(current_sequence), feature_count = nrow(features),
  segment_count = nrow(segments), qualifier_count = nrow(qualifiers),
  qualifier_link_count = nrow(qualifier_links), primer_count = 0L,
  primer_binding_count = 0L, stringsAsFactors = FALSE
)
source_releases <- rbind(legacy_provenance$legacy_release, current_release)

segment_counts <- table(segments$common_feature_id)
stopifnot(all(as.integer(segment_counts[features$common_feature_id]) ==
  features$segment_count))
segment_groups <- split(segments,
  factor(segments$common_feature_id, levels = features$common_feature_id))
segment_min <- vapply(segment_groups, function(x)
  min(x$start_1based_inclusive), numeric(1))
segment_max <- vapply(segment_groups, function(x)
  max(x$end_1based_inclusive), numeric(1))
standard_length <- vapply(segment_groups, function(x)
  sum(x$length_bp[x$segment_type == "standard"]), numeric(1))
gap_length <- vapply(segment_groups, function(x)
  sum(x$length_bp[x$segment_type == "gap"]), numeric(1))
stopifnot(
  all(segment_min == features$start_min_1based_inclusive),
  all(segment_max == features$end_max_1based_inclusive),
  all(standard_length == features$standard_length_bp),
  all(gap_length == features$gap_length_bp),
  max(segments$end_1based_inclusive) <= nchar(current_sequence)
)
qualifier_keys <- paste(qualifiers$common_feature_id,
  qualifiers$qualifier_name, qualifiers$value_index, sep = "\r")
link_keys <- paste(qualifier_links$common_feature_id,
  qualifier_links$qualifier_name, qualifier_links$value_index, sep = "\r")
stopifnot(
  !anyDuplicated(qualifier_keys),
  all(link_keys %in% qualifier_keys),
  !anyDuplicated(paste(link_keys, qualifier_links$link_index, sep = "\r"))
)
observed_link_counts <- table(factor(link_keys, levels = qualifier_keys))
stopifnot(all(as.integer(observed_link_counts) ==
  qualifiers$url_occurrence_count))
exact_protein_ids <- features$common_feature_id[
  features$detectionMode == "exactProteinMatch" &
    !is.na(features$detectionMode)]
translation_ids <- unique(qualifiers$common_feature_id[
  qualifiers$qualifier_name == "translation"])
stopifnot(length(exact_protein_ids) == 250L,
  all(exact_protein_ids %in% translation_ids))

# The legacy wide workbook is deliberately absent from the repository. Its
# audited, compact migration result is frozen in legacy_provenance above.
# The former workbook-to-correspondence derivation is retained here only as
# historical specification for reproducing that audit from an external copy.
if (FALSE) {
value_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "<NA>"
  x
}

derive_source_dna <- function(feature_ids, segment_table, sequence,
                              start_col, end_col) {
  split_segments <- split(segment_table,
    factor(segment_table$feature_id, levels = feature_ids))
  stats::setNames(vapply(split_segments, function(rows) {
    rows <- rows[order(rows$segment_index), , drop = FALSE]
    rows <- rows[rows$segment_type == "standard", , drop = FALSE]
    paste0(mapply(function(start, end) substr(sequence, start, end),
      as.integer(rows[[start_col]]), as.integer(rows[[end_col]]),
      USE.NAMES = FALSE), collapse = "")
  }, character(1)), feature_ids)
}

record_index <- function(feature_table, segment_table, qualifier_table,
                         sequence, legacy = FALSE) {
  feature_ids <- as.character(feature_table$feature_id)
  start_col <- if (legacy) "start" else "start_1based_inclusive"
  end_col <- if (legacy) "end" else "end_1based_inclusive"
  dna <- derive_source_dna(feature_ids, segment_table, sequence,
    start_col, end_col)
  if (legacy) {
    stopifnot(identical(unname(dna),
      toupper(as.character(feature_table$reference_dna_top_strand))))
  }
  segment_groups <- split(segment_table,
    factor(segment_table$feature_id, levels = feature_ids))
  structure <- vapply(segment_groups, function(rows) {
    rows <- rows[order(rows$segment_index), , drop = FALSE]
    paste(value_text(rows$segment_type), value_text(rows$length_bp),
      value_text(rows$color), value_text(rows$translated),
      value_text(rows$segment_name),
      value_text(rows$translationNumberingStartsFrom), sep = ":",
      collapse = "|")
  }, character(1))
  q <- qualifier_table
  if (legacy) {
    names(q)[match(c("value_text", "value_int", "value_predef",
      "value_bool"), names(q))] <- c("text", "int", "predef", "bool")
  }
  qualifier_groups <- split(q,
    factor(q$feature_id, levels = feature_ids))
  qualifier_signature <- vapply(qualifier_groups, function(rows) {
    if (!nrow(rows)) return("")
    rows <- rows[order(rows$qualifier_name, rows$value_index), , drop = FALSE]
    paste(value_text(rows$qualifier_name), value_text(rows$value_index),
      value_text(rows$text), value_text(rows$int), value_text(rows$predef),
      value_text(rows$bool), value_text(rows$value_display), sep = ":",
      collapse = "|")
  }, character(1))
  protein <- if (legacy) {
    as.character(feature_table$reference_protein)
  } else {
    vapply(qualifier_groups, function(rows) {
      hit <- rows$qualifier_name == "translation"
      if (!any(hit)) return(NA_character_)
      rows <- rows[hit, , drop = FALSE]
      rows <- rows[order(rows$value_index), , drop = FALSE]
      value <- as.character(rows$text[1L])
      if (is.na(value) || !nzchar(value)) as.character(rows$value_display[1L])
      else value
    }, character(1))
  }
  # maxRunOn/maxFusedRunOn are coordinates in the release-local backing
  # sequence. They remain stored source attributes, but moving a record inside
  # that backing sequence is not a biological feature change.
  comparable_fields <- c("directionality_raw", "directionality_label",
    "strand_symbol", "detectionMode", "translationMW", "readingFrame",
    "cleavageArrows", "allowSegmentOverlaps",
    "consecutiveTranslationNumbering", "consecutiveNumberingStartsFrom",
    "swappedSegmentNumbering", "translateFirstCodonAsMet", "hitsStopCodon",
    "originalName", "originalSequence", "isFavorite", "geneticCode",
    "prioritize")
  comparable <- vapply(seq_len(nrow(feature_table)), function(i) {
    paste(vapply(comparable_fields, function(field) {
      if (field %in% names(feature_table)) value_text(feature_table[[field]][i])
      else "<NA>"
    }, character(1)), collapse = "|")
  }, character(1))
  full_signature <- vapply(seq_len(nrow(feature_table)), function(i) {
    digest::digest(paste(feature_table$name[i], feature_table$type[i], dna[i],
      value_text(protein[i]), structure[i], qualifier_signature[i],
      comparable[i], sep = "\r"), algo = "sha256", serialize = FALSE)
  }, character(1))
  data.frame(
    source_feature_id = feature_ids,
    source_recent_id = as.integer(feature_table$recentID),
    name = as.character(feature_table$name),
    type = as.character(feature_table$type),
    directionality_label = as.character(feature_table$directionality_label),
    detectionMode = as.character(feature_table$detectionMode),
    segment_count = as.integer(feature_table$segment_count),
    dna_sha256 = vapply(dna, digest::digest, character(1), algo = "sha256",
      serialize = FALSE),
    protein_sha256 = vapply(protein, function(x) {
      if (is.na(x) || !nzchar(x)) NA_character_ else
        digest::digest(x, algo = "sha256", serialize = FALSE)
    }, character(1)),
    structure_signature = structure,
    qualifier_signature = qualifier_signature,
    feature_attribute_signature = comparable,
    full_signature = full_signature,
    stringsAsFactors = FALSE
  )
}

legacy_index <- record_index(legacy_features, legacy_segments,
  legacy_qualifiers, legacy_sequence, legacy = TRUE)
current_index <- record_index(raw_features, raw_segments, raw_qualifiers,
  current_sequence)
current_index$common_feature_id <- unname(id_map[current_index$source_feature_id])

old_matched <- rep(FALSE, nrow(legacy_index))
new_matched <- rep(FALSE, nrow(current_index))
correspondence_rows <- list()

add_correspondence <- function(old_i = NA_integer_, new_i = NA_integer_,
                               relationship, confidence, evidence) {
  changed <- function(field) {
    if (is.na(old_i) || is.na(new_i)) return(NA)
    old <- legacy_index[[field]][old_i]
    new <- current_index[[field]][new_i]
    if (is.na(old) && is.na(new)) return(FALSE)
    !isTRUE(old == new)
  }
  correspondence_rows[[length(correspondence_rows) + 1L]] <<- data.frame(
    old_source_release_id = if (is.na(old_i)) NA_character_ else
      legacy_source_release_id,
    old_source_feature_id = if (is.na(old_i)) NA_character_ else
      legacy_index$source_feature_id[old_i],
    old_name = if (is.na(old_i)) NA_character_ else legacy_index$name[old_i],
    old_type = if (is.na(old_i)) NA_character_ else legacy_index$type[old_i],
    new_source_release_id = if (is.na(new_i)) NA_character_ else
      source_release_id,
    new_source_feature_id = if (is.na(new_i)) NA_character_ else
      current_index$source_feature_id[new_i],
    common_feature_id = if (is.na(new_i)) NA_character_ else
      current_index$common_feature_id[new_i],
    new_name = if (is.na(new_i)) NA_character_ else current_index$name[new_i],
    new_type = if (is.na(new_i)) NA_character_ else current_index$type[new_i],
    relationship = relationship, confidence = confidence,
    evidence = evidence,
    name_changed = changed("name"),
    type_changed = changed("type"),
    directionality_changed = changed("directionality_label"),
    detection_mode_changed = changed("detectionMode"),
    dna_changed = changed("dna_sha256"),
    protein_changed = changed("protein_sha256"),
    segment_structure_changed = changed("structure_signature"),
    qualifiers_changed = changed("qualifier_signature"),
    feature_attributes_changed = changed("feature_attribute_signature"),
    stringsAsFactors = FALSE
  )
}

pair_groups <- function(old_groups, new_groups, renamed = FALSE) {
  shared <- intersect(names(old_groups), names(new_groups))
  for (group in shared) {
    old_ids <- old_groups[[group]]
    new_ids <- new_groups[[group]]
    old_ids <- old_ids[!old_matched[old_ids]]
    new_ids <- new_ids[!new_matched[new_ids]]
    while (length(old_ids) && length(new_ids)) {
      candidates <- expand.grid(old_i = old_ids, new_i = new_ids)
      candidates$score <- mapply(function(old_i, new_i) {
        old <- legacy_index[old_i, , drop = FALSE]
        new <- current_index[new_i, , drop = FALSE]
        200L * isTRUE(old$full_signature == new$full_signature) +
          100L * isTRUE(old$dna_sha256 == new$dna_sha256) +
          30L * isTRUE(old$type == new$type) +
          20L * (!is.na(old$protein_sha256) &&
            isTRUE(old$protein_sha256 == new$protein_sha256)) +
          10L * isTRUE(old$structure_signature == new$structure_signature) +
          10L * isTRUE(old$qualifier_signature == new$qualifier_signature) +
          5L * identical(old$directionality_label, new$directionality_label) +
          5L * identical(old$detectionMode, new$detectionMode)
      }, candidates$old_i, candidates$new_i)
      candidates <- candidates[order(-candidates$score,
        legacy_index$source_recent_id[candidates$old_i],
        current_index$source_recent_id[candidates$new_i]), , drop = FALSE]
      chosen <- candidates[1L, , drop = FALSE]
      old_i <- chosen$old_i; new_i <- chosen$new_i
      old <- legacy_index[old_i, , drop = FALSE]
      new <- current_index[new_i, , drop = FALSE]
      relationship <- if (!isTRUE(old$type == new$type)) "reclassified" else
        if (renamed && !isTRUE(old$name == new$name)) "renamed" else
        if (isTRUE(old$full_signature == new$full_signature)) "unchanged" else
        "modified"
      tied <- sum(candidates$score == chosen$score) > 1L
      if (tied) relationship <- "unresolved"
      confidence <- if (tied) "low" else if (
        isTRUE(old$dna_sha256 == new$dna_sha256)) "high" else "medium"
      evidence <- paste(c(
        if (isTRUE(old$dna_sha256 == new$dna_sha256))
          "exact_standard_segment_dna",
        if (!is.na(old$protein_sha256) &&
            isTRUE(old$protein_sha256 == new$protein_sha256))
          "exact_translation",
        if (isTRUE(old$structure_signature == new$structure_signature))
          "exact_segment_structure",
        if (isTRUE(old$qualifier_signature == new$qualifier_signature))
          "exact_qualifiers",
        if (tied) "ambiguous_same_name_group"), collapse = ";")
      add_correspondence(old_i, new_i, relationship, confidence, evidence)
      old_matched[old_i] <<- TRUE; new_matched[new_i] <<- TRUE
      old_ids <- old_ids[old_ids != old_i]
      new_ids <- new_ids[new_ids != new_i]
    }
  }
}

old_name_groups <- split(seq_len(nrow(legacy_index)), legacy_index$name)
new_name_groups <- split(seq_len(nrow(current_index)), current_index$name)
pair_groups(old_name_groups, new_name_groups)
old_case_groups <- split(seq_len(nrow(legacy_index)),
  tolower(legacy_index$name))
new_case_groups <- split(seq_len(nrow(current_index)),
  tolower(current_index$name))
pair_groups(old_case_groups, new_case_groups, renamed = TRUE)

split_old <- which(!old_matched & legacy_index$name ==
  "tPA signal/pro sequence")
split_new <- which(!new_matched & current_index$name %in%
  c("tPA signal sequence", "tPA propeptide sequence"))
stopifnot(length(split_old) == 1L, length(split_new) == 2L)
for (new_i in split_new) {
  add_correspondence(split_old, new_i, "split", "high",
    "curated_tPA_signal_propeptide_split")
  new_matched[new_i] <- TRUE
}
old_matched[split_old] <- TRUE

for (new_i in which(!new_matched)) {
  add_correspondence(new_i = new_i, relationship = "added",
    confidence = "high", evidence = "no_name_or_curated_predecessor")
}
for (old_i in which(!old_matched)) {
  add_correspondence(old_i = old_i, relationship = "removed",
    confidence = "high", evidence = "no_name_or_curated_successor")
}
feature_correspondence <- do.call(rbind, correspondence_rows)
stopifnot(
  all(current_index$source_feature_id %in%
    feature_correspondence$new_source_feature_id),
  all(legacy_index$source_feature_id %in%
    feature_correspondence$old_source_feature_id)
)

inherited_identity_ids <- unique(feature_correspondence$common_feature_id[
  !is.na(feature_correspondence$old_source_feature_id) &
    feature_correspondence$relationship != "unresolved"])
feature_identities <- data.frame(
  common_feature_id = current_index$common_feature_id,
  canonical_name = current_index$name,
  canonical_type = current_index$type,
  created_in_release = ifelse(current_index$common_feature_id %in%
    inherited_identity_ids, legacy_source_release_id, source_release_id),
  status = "active",
  stringsAsFactors = FALSE
)

legacy_feature_records <- legacy_index[c("source_feature_id",
  "source_recent_id", "name", "type", "directionality_label",
  "detectionMode", "segment_count", "dna_sha256", "protein_sha256")]
legacy_feature_records$source_release_id <- legacy_source_release_id
legacy_feature_records <- legacy_feature_records[c("source_release_id",
  setdiff(names(legacy_feature_records), "source_release_id"))]
}

legacy_feature_records <- legacy_provenance$legacy_feature_records
feature_correspondence <- legacy_provenance$feature_correspondence
stopifnot(
  setequal(as.character(raw_features$feature_id),
    as.character(na.omit(feature_correspondence$new_source_feature_id))),
  all(legacy_feature_records$source_release_id == legacy_source_release_id)
)
inherited_identity_ids <- unique(feature_correspondence$common_feature_id[
  !is.na(feature_correspondence$old_source_feature_id) &
    feature_correspondence$relationship != "unresolved"])
feature_identities <- data.frame(
  common_feature_id = features$common_feature_id,
  canonical_name = features$name,
  canonical_type = features$type,
  created_in_release = ifelse(features$common_feature_id %in%
    inherited_identity_ids, legacy_source_release_id, source_release_id),
  status = "active",
  stringsAsFactors = FALSE
)

read_uint32_be <- function(bytes) {
  values <- as.numeric(bytes)
  values[1L] * 256^3 + values[2L] * 256^2 + values[3L] * 256 + values[4L]
}

value_or <- function(x, fallback) {
  if (length(x) == 0L || is.na(x) || !nzchar(x)) fallback else x
}

read_binary_plasmid_reference <- function(path) {
  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  position <- 1L
  sequence <- NULL
  feature_xml <- NULL
  primer_xml <- NULL
  enzyme_display_packet <- NULL
  custom_enzyme_xml <- NULL
  while (position + 4L <= length(bytes)) {
    block_id <- as.integer(bytes[position])
    block_length <- read_uint32_be(bytes[(position + 1L):(position + 4L)])
    first <- position + 5L
    last <- first + block_length - 1L
    if (last > length(bytes)) stop("Invalid SnapGene block in ", path,
      call. = FALSE)
    payload <- bytes[first:last]
    if (block_id == 0L) sequence <- rawToChar(payload[-1L])
    if (block_id == 10L) feature_xml <- rawToChar(payload)
    if (block_id == 5L) primer_xml <- rawToChar(payload)
    if (block_id == 13L) enzyme_display_packet <- payload
    if (block_id == 14L) custom_enzyme_xml <- rawToChar(payload)
    position <- last + 1L
  }
  if (is.null(sequence) || is.null(feature_xml)) {
    stop("SnapGene reference lacks sequence or feature data: ", path,
      call. = FALSE)
  }
  reference_id <- make.names(tools::file_path_sans_ext(basename(path)))
  document <- xml2::read_xml(feature_xml)
  nodes <- xml2::xml_find_all(document, ".//Feature")
  strand_map <- c(`0` = ".", `1` = "+", `2` = "-", `3` = "+/-")
  sequence_length <- nchar(sequence)
  circular_envelope <- function(segment_rows) {
    starts <- segment_rows$start
    ends <- segment_rows$end
    candidates <- lapply(starts, function(origin) {
      shifted_start <- starts
      shifted_start[shifted_start < origin] <-
        shifted_start[shifted_start < origin] + sequence_length
      width <- ifelse(starts <= ends, ends - starts,
        sequence_length - starts + ends)
      shifted_end <- shifted_start + width
      c(start = origin, end = max(shifted_end), span = max(shifted_end) - origin)
    })
    chosen <- candidates[[which.min(vapply(candidates, `[[`, numeric(1L), "span"))]]
    c(start = ((chosen[["start"]] - 1L) %% sequence_length) + 1L,
      end = ((chosen[["end"]] - 1L) %% sequence_length) + 1L)
  }
  rows <- lapply(seq_along(nodes), function(i) {
    node <- nodes[[i]]
    attrs <- xml2::xml_attrs(node)
    segment_nodes <- xml2::xml_find_all(node, "./Segment")
    segment_rows <- lapply(seq_along(segment_nodes), function(j) {
      segment_attrs <- xml2::xml_attrs(segment_nodes[[j]])
      bounds <- as.integer(strsplit(segment_attrs[["range"]], "-",
        fixed = TRUE)[[1L]])
      data.frame(
        segment_index = j,
        segment_type = value_or(segment_attrs["type"], "standard"),
        start = bounds[1L], end = bounds[2L],
        length_bp = if (bounds[1L] <= bounds[2L]) {
          bounds[2L] - bounds[1L] + 1L
        } else sequence_length - bounds[1L] + bounds[2L] + 1L,
        translated = FALSE, segment_name = NA_character_,
        color = value_or(segment_attrs["color"], "#B8BDC3"),
        stringsAsFactors = FALSE
      )
    })
    reference_segments <- do.call(rbind, segment_rows)
    directionality <- value_or(attrs["directionality"], "0")
    envelope <- circular_envelope(reference_segments)
    cleavage_text <- value_or(attrs["cleavageArrows"], "")
    cleavage_arrows <- if (nzchar(cleavage_text)) {
      as.numeric(strsplit(cleavage_text, "[,; ]+")[[1L]])
    } else numeric()
    cleavage_arrows <- cleavage_arrows[is.finite(cleavage_arrows)]
    data.frame(
      reference_id = reference_id,
      common_feature_id = paste0("reference:", reference_id, ":",
        value_or(attrs["recentID"], sprintf("%04d", i))),
      name = unname(attrs["name"]), type = unname(attrs["type"]),
      start = envelope[["start"]], end = envelope[["end"]],
      strand = unname(value_or(strand_map[[directionality]], ".")),
      color = reference_segments$color[1L],
      cleavage_arrows = I(list(cleavage_arrows)),
      segments = I(list(reference_segments)), stringsAsFactors = FALSE
    )
  })
  empty_primers <- data.frame(
    reference_id = character(), primer_id = character(), name = character(),
    sequence = character(), start = integer(), end = integer(),
    strand = character(), annealed_bases = character(),
    melting_temperature = numeric(), stringsAsFactors = FALSE
  )
  primer_rows <- if (is.null(primer_xml)) empty_primers else {
    primer_document <- xml2::read_xml(primer_xml)
    primer_nodes <- xml2::xml_find_all(primer_document, ".//Primer")
    parsed <- unlist(lapply(seq_along(primer_nodes), function(i) {
      primer <- primer_nodes[[i]]
      primer_attrs <- xml2::xml_attrs(primer)
      # Simplified binding sites duplicate the canonical records in these
      # SnapGene files and are not independent biological annotations.
      sites <- xml2::xml_find_all(primer,
        "./BindingSite[not(@simplified='1')]")
      lapply(seq_along(sites), function(j) {
        site_attrs <- xml2::xml_attrs(sites[[j]])
        bounds <- as.integer(strsplit(site_attrs[["location"]], "-",
          fixed = TRUE)[[1L]])
        data.frame(
          reference_id = reference_id,
          primer_id = paste0("reference:", reference_id, ":primer:",
            value_or(primer_attrs["recentID"], sprintf("%04d", i)), ":", j),
          name = unname(primer_attrs["name"]),
          sequence = toupper(unname(primer_attrs["sequence"])),
          start = bounds[1L], end = bounds[2L],
          strand = if (value_or(site_attrs["boundStrand"], "0") == "0")
            "+" else "-",
          annealed_bases = toupper(unname(value_or(
            site_attrs["annealedBases"], ""))),
          melting_temperature = as.numeric(value_or(
            site_attrs["meltingTemperature"], NA_character_)),
          stringsAsFactors = FALSE
        )
      })
    }), recursive = FALSE)
    if (length(parsed)) do.call(rbind, parsed) else empty_primers
  }
  printable_runs <- if (is.null(enzyme_display_packet)) character() else {
    printable <- as.integer(enzyme_display_packet)
    printable[printable < 32L | printable > 126L] <- 32L
    runs <- strsplit(trimws(rawToChar(as.raw(printable))),
      "[[:space:]]{2,}")[[1L]]
    trimws(runs[nchar(trimws(runs)) >= 2L])
  }
  selected_enzyme_set <- if (length(printable_runs)) {
    printable_runs[which.max(nchar(printable_runs))]
  } else NA_character_
  custom_sets <- data.frame(name = character(), enzymes = I(list()))
  if (!is.null(custom_enzyme_xml)) {
    custom_document <- xml2::read_xml(custom_enzyme_xml)
    custom_nodes <- xml2::xml_find_all(custom_document, ".//CustomEnzymeSet")
    if (length(custom_nodes)) {
      custom_sets <- do.call(rbind, lapply(custom_nodes, function(node) {
        attrs <- xml2::xml_attrs(node)
        enzyme_text <- value_or(attrs["enzymeNames"], "")
        enzymes <- if (nzchar(enzyme_text)) {
          strsplit(trimws(enzyme_text), "[[:space:]]+")[[1L]]
        } else character()
        data.frame(name = value_or(attrs["name"], ""),
          enzymes = I(list(enzymes)), stringsAsFactors = FALSE)
      }))
    }
  }
  selected_custom <- match(selected_enzyme_set, custom_sets$name)
  enzyme_profile <- data.frame(
    reference_id = reference_id,
    set_name = selected_enzyme_set,
    profile_type = if (is.na(selected_enzyme_set)) "unknown" else if (
      identical(selected_enzyme_set, "None")) "none" else if (
      !is.na(selected_custom)) "custom" else "preset",
    enzymes = I(list(if (is.na(selected_custom)) character() else
      custom_sets$enzymes[[selected_custom]])),
    stringsAsFactors = FALSE
  )
  list(
    sequence = data.frame(
      reference_id = reference_id,
      label = tools::file_path_sans_ext(basename(path)),
      sequence = toupper(sequence),
      source_path = path,
      source_sha256 = digest::digest(path, algo = "sha256", file = TRUE,
        serialize = FALSE),
      stringsAsFactors = FALSE
    ),
    features = do.call(rbind, rows), primers = primer_rows,
    enzyme_profile = enzyme_profile
  )
}

binary_reference_paths <- file.path("examples", "plasmid", c(
  "pBR322.dna", "pUC19.dna", "pBluescript II SK(+).dna", "pSB1C3.dna",
  "pET-28a(+).dna", "pETDuet-1.dna", "pcDNA3.1(+).dna",
  "pTRE-Tight-BI.dna", "pSpCas9(BB)-2A-GFP (PX458).dna", "pDONR221.dna",
  "pCAMBIA1300.dna", "pEarleyGate 201.dna", "pTRIPZ.dna"
))
stopifnot(all(file.exists(binary_reference_paths)))
binary_references <- lapply(binary_reference_paths,
  read_binary_plasmid_reference)
reference_sequences <- do.call(rbind,
  lapply(binary_references, `[[`, "sequence"))
reference_features <- do.call(rbind,
  lapply(binary_references, `[[`, "features"))
reference_primers <- do.call(rbind,
  lapply(binary_references, `[[`, "primers"))
reference_enzyme_profiles <- do.call(rbind,
  lapply(binary_references, `[[`, "enzyme_profile"))
stopifnot(
  nrow(reference_sequences) == 13L,
  nrow(reference_features) == 201L,
  sum(vapply(reference_features$segments, nrow, integer(1L))) == 218L,
  nrow(reference_primers) == 7L,
  nrow(reference_enzyme_profiles) == 13L,
  sum(reference_enzyme_profiles$profile_type == "custom") == 2L,
  sum(reference_enzyme_profiles$profile_type == "none") == 1L,
  identical(reference_enzyme_profiles$set_name[
    reference_enzyme_profiles$reference_id == make.names(
      "pSpCas9(BB)-2A-GFP (PX458)")], "BbsI + EcoRI"),
  identical(reference_enzyme_profiles$set_name[
    reference_enzyme_profiles$reference_id == make.names("pTRIPZ")],
    "Unique Cutters + BamHI"),
  identical(sort(unique(reference_primers$reference_id)),
    sort(make.names(c("pSB1C3", "pETDuet-1"))))
)

ggchord_common_feature_database <- list(
  metadata = list(
    database_id = "ggchord-common-features",
    database_version = "snapgene-8.2.3-export15",
    source_release_id = source_release_id,
    snapgene_version = "8.2.3",
    source = as.character(source_metadata[["Source file"]]),
    source_workbook_sha256 = actual_sha256,
    source_file_sha256 = as.character(source_metadata[["SHA256"]]),
    export_version = as.integer(source_metadata[["export-format version"]]),
    import_version = as.integer(source_metadata[["import-format version"]]),
    record_counts = c(features = nrow(features), segments = nrow(segments),
      qualifiers = nrow(qualifiers), qualifier_links = nrow(qualifier_links),
      source_sequence_bp = nchar(current_sequence), primers = 0L,
      primer_bindings = 0L),
    reference_record_counts = c(
      plasmids = nrow(reference_sequences),
      features = nrow(reference_features),
      segments = sum(vapply(reference_features$segments, nrow, integer(1L))),
      primers = nrow(reference_primers),
      enzyme_profiles = nrow(reference_enzyme_profiles)
    ),
    generating_version = "0.13.0"
  ),
  source_releases = source_releases,
  source_sequences = source_sequences,
  feature_identities = feature_identities,
  features = features,
  segments = segments,
  qualifiers = qualifiers,
  qualifier_links = qualifier_links,
  feature_type_summary = feature_type_summary,
  legacy_feature_records = legacy_feature_records,
  feature_correspondence = feature_correspondence,
  source_audit = list(
    primer_settings = read_sheet("PrimerSettings"),
    alignment_settings = read_sheet("AlignmentSettings"),
    enzyme_settings = read_sheet("EnzymeSettings"),
    root_metadata = read_sheet("RootMetadata"),
    additional_properties = read_sheet("AdditionalProps"),
    notes = read_sheet("Notes"),
    packet_summary = read_sheet("PacketSummary"),
    completeness_audit = read_sheet("CompletenessAudit"),
    cross_validation = read_sheet("CrossValidation"),
    coordinate_conventions = read_sheet("CoordinateConventions")
  ),
  reference_sequences = reference_sequences,
  reference_features = reference_features,
  reference_primers = reference_primers,
  reference_enzyme_profiles = reference_enzyme_profiles,
  search_indexes = list(
    dna_feature_ids = sort(unique(segments$common_feature_id[
      segments$segment_type == "standard"])),
    protein_feature_ids = sort(unique(qualifiers$common_feature_id[
      qualifiers$qualifier_name == "translation"])),
    translated_feature_ids = sort(unique(segments$common_feature_id[
      !is.na(segments$translated) & segments$translated])),
    feature_name = split(features$common_feature_id,
      factor(tolower(features$name), levels = unique(tolower(features$name)))),
    feature_type = split(features$common_feature_id,
      factor(features$type, levels = unique(features$type)))
  )
)

old <- new.env(parent = emptyenv())
if (file.exists("R/sysdata.rda")) load("R/sysdata.rda", envir = old)
objects <- list(ggchord_common_feature_database = ggchord_common_feature_database)
if (exists("ggchord_rebase_database", envir = old, inherits = FALSE)) {
  objects$ggchord_rebase_database <- get("ggchord_rebase_database", envir = old)
}
if (exists("ggchord_primer_database", envir = old, inherits = FALSE)) {
  objects$ggchord_primer_database <- get("ggchord_primer_database", envir = old)
}
list2env(objects, envir = environment())
save(list = names(objects), file = "R/sysdata.rda", compress = "xz")
