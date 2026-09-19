# Sequencing-primer catalogue and deterministic binding search.

ggchord_builtin_primers <- function() {
  if (!exists("ggchord_primer_database", inherits = TRUE)) {
    ggchord_stop("The internal primer database is missing; reinstall ggchord")
  }
  get("ggchord_primer_database", inherits = TRUE)
}

#' Inspect the bundled sequencing-primer catalogue
#'
#' The catalogue stores one searchable record per distinct oligonucleotide.
#' Set `aliases = TRUE` to obtain the source rows, including multiple biological
#' names that share one sequence. The source direction is retained only as a
#' descriptive hint; binding strand is determined from the target sequence.
#'
#' @param set Return all records or only the commonly used universal subset.
#' @param aliases Return source alias records instead of the sequence catalogue.
#' @return A data frame.
#' @export
#' @examples
#' head(primer_catalog("universal"))
#' head(primer_catalog("universal", aliases = TRUE))
primer_catalog <- function(set = c("all", "universal"), aliases = FALSE) {
  set <- match.arg(set)
  if (!is.logical(aliases) || length(aliases) != 1L || is.na(aliases)) {
    ggchord_stop("primer_catalog(): aliases must be TRUE or FALSE")
  }
  database <- ggchord_builtin_primers()
  out <- if (aliases) database$aliases else database$catalog
  if (set == "universal") out <- out[out$is_universal, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "primer_database_metadata") <- database$metadata
  out
}

ggchord_primer_definitions <- function(primers, set) {
  if (is.null(primers)) {
    database <- ggchord_builtin_primers()
    catalog <- database$catalog
    aliases <- database$aliases
    if (set == "universal") {
      catalog <- catalog[catalog$is_universal, , drop = FALSE]
      aliases <- aliases[aliases$primer_id %in% catalog$primer_id, , drop = FALSE]
    }
    return(list(catalog = catalog, aliases = aliases))
  }
  if (is.character(primers)) {
    primers <- data.frame(
      name = names(primers) %||% paste0("primer_", seq_along(primers)),
      sequence = unname(primers), stringsAsFactors = FALSE
    )
  }
  if (!is.data.frame(primers)) {
    ggchord_stop("find_primer_bindings(): primers must be a data frame or named character vector")
  }
  ggchord_require_columns(primers, "sequence", "find_primer_bindings()")
  primers <- as.data.frame(primers, stringsAsFactors = FALSE)
  primers$sequence <- toupper(gsub("[[:space:]]", "", as.character(primers$sequence)))
  if (anyNA(primers$sequence) || any(!grepl("^[ACGTRYSWKMBDHVN]+$", primers$sequence))) {
    ggchord_stop("find_primer_bindings(): invalid primer DNA symbols")
  }
  if (!"name" %in% names(primers)) primers$name <- paste0("primer_", seq_len(nrow(primers)))
  if (!"primer_id" %in% names(primers)) {
    unique_sequences <- unique(primers$sequence)
    ids <- paste0("custom_primer_", seq_along(unique_sequences))
    primers$primer_id <- unname(stats::setNames(ids, unique_sequences)[primers$sequence])
  }
  aliases <- data.frame(
    primer_id = as.character(primers$primer_id),
    name = as.character(primers$name),
    sequence = primers$sequence,
    description = if ("description" %in% names(primers)) as.character(primers$description) else NA_character_,
    direction_hint = if ("direction_hint" %in% names(primers)) as.character(primers$direction_hint) else NA_character_,
    is_universal = FALSE, source_sheet = "custom", source_row = seq_len(nrow(primers)),
    stringsAsFactors = FALSE
  )
  first <- !duplicated(aliases$primer_id)
  catalog <- data.frame(
    primer_id = aliases$primer_id[first], preferred_name = aliases$name[first],
    sequence = aliases$sequence[first], length = nchar(aliases$sequence[first]),
    is_universal = FALSE, source = "custom", stringsAsFactors = FALSE
  )
  if (anyDuplicated(catalog$primer_id) || anyDuplicated(catalog$sequence)) {
    ggchord_stop("find_primer_bindings(): each primer_id must identify one distinct sequence")
  }
  list(catalog = catalog, aliases = aliases)
}

ggchord_primer_empty <- function() {
  data.frame(
    accver = character(), primer_id = character(), name = character(),
    sequence = character(), start = integer(), end = integer(),
    crosses_origin = logical(), strand = character(),
    annealed_bases = character(), annealed_length = integer(),
    mismatch_count = integer(), melting_temperature = numeric(),
    tm_method = character(), match_type = character(),
    binding_count = integer(), is_unique = logical(),
    stringsAsFactors = FALSE
  )
}

ggchord_primer_tm_wallace <- function(sequence) {
  chars <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  2 * sum(chars %in% c("A", "T")) + 4 * sum(chars %in% c("G", "C"))
}

ggchord_reference_primer_bindings <- function(ids, sequences) {
  database <- ggchord_builtin_common_features()
  if (!all(c("reference_sequences", "reference_primers") %in% names(database))) {
    return(ggchord_primer_empty())
  }
  rows <- list()
  for (i in seq_along(sequences)) {
    reference_id <- database$reference_sequences$reference_id[
      database$reference_sequences$sequence == sequences[i]
    ]
    if (!length(reference_id)) next
    hit <- database$reference_primers[
      database$reference_primers$reference_id %in% reference_id, , drop = FALSE
    ]
    if (!nrow(hit)) next
    # Databases generated before the v0.13 primer-geometry revision retained
    # SnapGene's zero-based binding-site coordinates.  Convert at the runtime
    # boundary as well so installed development data remains correct until the
    # next deterministic data regeneration.
    metadata_version <- database$metadata$primer_coordinate_system %||%
      "snapgene_zero_based_inclusive"
    if (identical(metadata_version, "snapgene_zero_based_inclusive")) {
      hit$start <- hit$start + 1L
      hit$end <- hit$end + 1L
    }
    hit$accver <- ids[i]
    hit$crosses_origin <- hit$end < hit$start
    hit$annealed_length <- nchar(hit$annealed_bases)
    hit$mismatch_count <- 0L
    hit$tm_method <- "source"
    hit$match_type <- ifelse(hit$annealed_length == nchar(hit$sequence),
      "reference_exact", "reference_partial")
    hit$binding_count <- ave(hit$primer_id, hit$primer_id, FUN = length)
    hit$is_unique <- hit$binding_count == 1L
    rows[[length(rows) + 1L]] <- hit
  }
  if (!length(rows)) return(ggchord_primer_empty())
  out <- do.call(rbind, rows)
  out <- out[, names(ggchord_primer_empty()), drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Find sequencing-primer binding sites
#'
#' Searches each distinct oligonucleotide once on both target strands and then
#' attaches all catalogue aliases. Exact matching is the default. Optional
#' `three_prime` matching uses the longest exact 3-prime suffix available, with
#' no mismatches, and is intended for inspecting partial annealing rather than
#' silently treating it as a full match.
#'
#' @param sequence DNA text, a named character vector, or a data frame with
#'   `accver` and `sequence`.
#' @param primers Optional custom primers with `sequence` and optional `name`
#'   and `primer_id`, or a named character vector.
#' @param set Bundled catalogue subset. `reference` returns authoritative primer
#'   annotations from a recognized bundled reference plasmid.
#' @param match Full exact matching or longest exact 3-prime annealing.
#' @param min_annealed_bases Minimum suffix length for `three_prime` matching.
#' @param circular Allow matches to cross the sequence origin.
#' @param unique_only Keep only primer sequences with one binding site per target.
#' @return One row per binding site. The `aliases` list-column preserves every
#'   catalogue name associated with the matched oligonucleotide.
#' @export
#' @examples
#' target <- c(circle = "TTGCAAAAACG")
#' find_primer_bindings(
#'   target, primers = c(test = "ACGTTGCA"),
#'   min_annealed_bases = 8
#' )
find_primer_bindings <- function(sequence, primers = NULL,
                                 set = c("universal", "all", "reference"),
                                 match = c("exact", "three_prime"),
                                 min_annealed_bases = 15L,
                                 circular = TRUE, unique_only = FALSE) {
  set <- match.arg(set)
  match <- match.arg(match)
  normalized <- ggchord_normalize_sequence_input(sequence, "find_primer_bindings()")
  if (!is.logical(circular) || length(circular) != 1L || is.na(circular) ||
      !is.logical(unique_only) || length(unique_only) != 1L || is.na(unique_only)) {
    ggchord_stop("find_primer_bindings(): circular and unique_only must be TRUE or FALSE")
  }
  min_annealed_bases <- as.integer(min_annealed_bases)
  if (length(min_annealed_bases) != 1L || is.na(min_annealed_bases) || min_annealed_bases < 8L) {
    ggchord_stop("find_primer_bindings(): min_annealed_bases must be an integer of at least 8")
  }
  if (set == "reference") {
    if (!is.null(primers)) ggchord_stop("find_primer_bindings(): primers cannot be combined with set = 'reference'")
    out <- ggchord_reference_primer_bindings(normalized$ids, normalized$sequences)
    out$aliases <- I(lapply(out$name, function(x) x))
    return(if (unique_only) out[out$is_unique, , drop = FALSE] else out)
  }
  definitions <- ggchord_primer_definitions(primers, set)
  rows <- list()
  for (s in seq_along(normalized$sequences)) {
    target <- normalized$sequences[s]
    target_length <- nchar(target)
    for (d in seq_len(nrow(definitions$catalog))) {
      full <- definitions$catalog$sequence[d]
      lengths <- if (match == "exact") nchar(full) else
        seq.int(nchar(full), min(min_annealed_bases, nchar(full)), by = -1L)
      primer_hits <- list()
      for (annealed_length in lengths) {
        query <- substr(full, nchar(full) - annealed_length + 1L, nchar(full))
        found <- list()
        for (strand in c("+", "-")) {
          pattern <- if (strand == "+") query else ggchord_reverse_complement(query)
          search <- if (circular && annealed_length > 1L)
            paste0(target, substr(target, 1L, annealed_length - 1L)) else target
          starts <- ggchord_restriction_match_starts(search, pattern, target_length)
          if (length(starts)) {
            found[[length(found) + 1L]] <- data.frame(
              accver = normalized$ids[s], primer_id = definitions$catalog$primer_id[d],
              name = definitions$catalog$preferred_name[d], sequence = full,
              start = starts, end = ((starts + annealed_length - 2L) %% target_length) + 1L,
              crosses_origin = starts + annealed_length - 1L > target_length,
              strand = strand, annealed_bases = query,
              annealed_length = annealed_length, mismatch_count = 0L,
              melting_temperature = ggchord_primer_tm_wallace(query),
              tm_method = "Wallace", match_type = if (annealed_length == nchar(full)) "exact" else "three_prime",
              stringsAsFactors = FALSE
            )
          }
        }
        if (length(found)) {
          primer_hits <- found
          break
        }
      }
      if (length(primer_hits)) rows <- c(rows, primer_hits)
    }
  }
  if (!length(rows)) {
    out <- ggchord_primer_empty()
    out$aliases <- I(list())
    return(out)
  }
  out <- do.call(rbind, rows)
  key <- paste(out$accver, out$primer_id, sep = "\r")
  out$binding_count <- ave(seq_along(key), key, FUN = length)
  out$is_unique <- out$binding_count == 1L
  alias_map <- split(definitions$aliases$name, definitions$aliases$primer_id)
  out$aliases <- I(lapply(out$primer_id, function(id) unique(as.character(alias_map[[id]]))))
  if (unique_only) out <- out[out$is_unique, , drop = FALSE]
  rownames(out) <- NULL
  out
}
