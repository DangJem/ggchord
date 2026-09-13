# Build the internal normalized common-feature database from the audited
# workbook and the audited annotations embedded in the example SnapGene files.
# readxl, digest, and xml2 are development-only dependencies.

workbook <- "examples/standardCommonFeatures_tables.xlsx"
expected_sha256 <- "3d9e7a78c705de2beeb9d8af3c74dede714faa468cd6eee61f78e0f784b65312"
if (!requireNamespace("readxl", quietly = TRUE) ||
    !requireNamespace("digest", quietly = TRUE) ||
    !requireNamespace("xml2", quietly = TRUE)) {
  stop("Install readxl, digest, and xml2 to regenerate the common-feature database",
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
    features = do.call(rbind, rows), primers = primer_rows
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
stopifnot(
  nrow(reference_sequences) == 13L,
  nrow(reference_features) == 201L,
  sum(vapply(reference_features$segments, nrow, integer(1L))) == 218L,
  nrow(reference_primers) == 7L,
  identical(sort(unique(reference_primers$reference_id)),
    sort(make.names(c("pSB1C3", "pETDuet-1"))))
)

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
    reference_record_counts = c(
      plasmids = nrow(reference_sequences),
      features = nrow(reference_features),
      segments = sum(vapply(reference_features$segments, nrow, integer(1L))),
      primers = nrow(reference_primers)
    ),
    generating_version = "0.13.0"
  ),
  features = features,
  segments = segments,
  qualifiers = qualifiers,
  qualifier_links = qualifier_links,
  feature_type_summary = feature_type_summary,
  reference_sequences = reference_sequences,
  reference_features = reference_features,
  reference_primers = reference_primers,
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
