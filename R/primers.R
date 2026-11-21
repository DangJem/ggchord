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
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' catalog <- primer_catalog(set="universal")
#' catalog[c("primer_id", "name", "sequence")]
#' primer_catalog(set="universal", aliases=TRUE)
primer_catalog <- function(set = "all", aliases = FALSE) {
  set <- match.arg(set, c("all", "universal"))
  if (!is.logical(aliases) || length(aliases) != 1L || is.na(aliases)) {
    ggchord_stop("primer_catalog(): aliases must be TRUE or FALSE")
  }
  database <- ggchord_builtin_primers()
  out <- if (aliases) database$aliases else database$catalog
  if (set == "universal") out <- out[out$is_universal, , drop = FALSE]
  rownames(out) <- NULL
  out$name <- out$name %||% out$preferred_name
  if (!aliases) {
    names_by_id <- split(database$aliases$name, database$aliases$primer_id)
    out$aliases <- I(lapply(out$primer_id, function(id) unique(names_by_id[[id]])))
  }
  ggchord_attach_data_report(out, "primer_catalog", metadata = database$metadata)
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
    primers$primer_id <- unname(setNames(ids, unique_sequences)[primers$sequence])
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
  if (anyDuplicated(unique(aliases[c("primer_id", "sequence")])$primer_id) || anyDuplicated(catalog$sequence)) {
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
    # Legacy zero-based binding-site coordinates. Convert at the runtime
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
#' Binding hits are sequence-search results, not experimental validation.
#' `melting_temperature` uses the Wallace composition estimate unless source
#' annotations supply it; ambiguous bases do not contribute to that estimate.
#'
#' @param sequence DNA text, a named character vector, or a data frame with
#'   `accver` and `sequence`.
#' @param primers Optional custom primers with `sequence` and optional `name`
#'   and `primer_id` (one ID per distinct sequence), or a named character vector.
#' @param set Bundled catalogue subset. `reference` returns authoritative primer
#'   annotations from a recognized bundled reference plasmid.
#' @param match Full exact matching or longest exact 3-prime annealing.
#' @param min_annealed_bases One finite integer of at least 8, the minimum
#'   suffix length for `three_prime` matching. Fractional values are rejected.
#' @param circular Allow matches to cross the sequence origin.
#' @param unique_only Keep only primer sequences with one binding site per target.
#' @return One row per binding site. The `aliases` list-column preserves every
#'   catalogue name associated with the matched oligonucleotide.
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' target <- c(A="TTGCAAAAACG")
#' bindings <- find_primer_bindings(target, primers=c(test="ACGTTGCA"),
#'                                 circular=TRUE, min_annealed_bases=8)
#' ref <- as_ggchord_seq(data.frame(accver="A", length=nchar(target)), sequence_unit="bp", topology="circular")
#' ggchord(ref) + geom_primer(data=bindings) + coord_circular()
find_primer_bindings <- function(sequence, primers = NULL, set = "universal",
    match = "exact", min_annealed_bases = 15L, circular = NULL, unique_only = FALSE) {
  caller <- "find_primer_bindings()"
  set <- match.arg(set, c("universal", "all", "reference"))
  match <- match.arg(match, c("exact", "three_prime"))
  normalized <- ggchord_normalize_sequence_input(sequence, caller)
  topology <- ggchord_search_topology(sequence, circular, normalized$ids, caller)
  if (!is.logical(unique_only) || length(unique_only) != 1L || is.na(unique_only)) {
    ggchord_stop(caller, ": unique_only must be TRUE or FALSE")
  }
  finish <- function(out) {
    if (unique_only) out <- out[out$is_unique, , drop = FALSE]
    out$feature_id <- paste(out$accver, out$primer_id, out$strand, out$start, sep = ":")[seq_len(nrow(out))]
    out$label <- out$name
    out$matched_bases <- out$annealed_length
    out$segments <- I(ggchord_segment_search_matches(out$start, out$end,
      nchar(normalized$sequences[match(out$accver, normalized$ids)])))
    across_origin <- which(out$crosses_origin)
    out$start[across_origin] <- 1
    out$end[across_origin] <- nchar(normalized$sequences[match(out$accver[across_origin], normalized$ids)])
    ggchord_attach_data_report(out, "find_primer_bindings", inputs = normalized$ids,
      metadata = list(set = set, match = match, min_annealed_bases = min_annealed_bases,
        topology = topology, temperature_method = "Wallace or source record"))
  }
  if (!is.numeric(min_annealed_bases) || length(min_annealed_bases) != 1L ||
      !is.finite(min_annealed_bases) || min_annealed_bases < 8L ||
      min_annealed_bases != floor(min_annealed_bases) || min_annealed_bases > .Machine$integer.max) {
    ggchord_stop("find_primer_bindings(): min_annealed_bases must be an integer of at least 8")
  }
  min_annealed_bases <- as.integer(min_annealed_bases)
  if (set == "reference") {
    if (!is.null(primers)) ggchord_stop("find_primer_bindings(): primers cannot be combined with set = 'reference'")
    references <- ggchord_builtin_common_features()$reference_sequences$sequence
    if (any(!normalized$sequences %in% references))
      ggchord_stop(caller, ": reference mode requires matching complete source characters")
    out <- ggchord_reference_primer_bindings(normalized$ids, normalized$sequences)
    out <- out[!out$crosses_origin | topology[match(out$accver, normalized$ids)], , drop = FALSE]
    out$aliases <- I(lapply(out$name, function(x) x))
    return(finish(out))
  }
  definitions <- ggchord_primer_definitions(primers, set)
  rows <- list()
  for (s in seq_along(normalized$sequences)) {
    target <- normalized$sequences[s]
    circular <- topology[s]
    target_length <- nchar(target)
    for (d in seq_len(nrow(definitions$catalog))) {
      full <- definitions$catalog$sequence[d]
      if (match == "three_prime" && nchar(full) < min_annealed_bases) next
      lengths <- if (match == "exact") nchar(full) else
        seq.int(nchar(full), min_annealed_bases, by = -1L)
      lengths <- lengths[lengths <= target_length]
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
    return(finish(out))
  }
  out <- do.call(rbind, rows)
  key <- paste(out$accver, out$primer_id, sep = "\r")
  out$binding_count <- ave(seq_along(key), key, FUN = length)
  out$is_unique <- out$binding_count == 1L
  alias_map <- split(definitions$aliases$name, definitions$aliases$primer_id)
  out$aliases <- I(lapply(out$primer_id, function(id) unique(as.character(alias_map[[id]]))))
  rownames(out) <- NULL
  finish(out)
}
