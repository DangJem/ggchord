# CDS coordinate maps use the GFF3 phase definition, not genomic frame.
# Specification: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

ggchord_project_integer <- function(x, name, caller, minimum = 1) {
  if (!is.numeric(x) || any(!is.finite(x)) || any(x < minimum | x != floor(x)))
    ggchord_stop(caller, ": ", name, " must contain finite integer values >= ", minimum)
  out <- as.numeric(x)
  names(out) <- names(x)
  out
}

ggchord_project_model <- function(cds, caller) {
  if (!inherits(cds, "ggchord_cds") || !is.data.frame(cds$models) || !is.list(cds$segments))
    ggchord_stop(caller, ": cds must be produced by build_ggchord_cds()")
  cds
}

ggchord_project_parent <- function(data, caller) {
  if ("parent_id" %in% names(data)) {
    out <- as.character(data$parent_id)
    if ("parent_ids" %in% names(data)) {
      if (!is.list(data$parent_ids) || any(lengths(data$parent_ids) != 1L) ||
          !identical(out, unname(vapply(data$parent_ids, as.character, character(1)))))
        ggchord_stop(caller, ": parent_id and single-valued parent_ids must agree")
    }
  } else if ("parent_ids" %in% names(data) && is.list(data$parent_ids) &&
      all(lengths(data$parent_ids) == 1L)) {
    out <- unname(vapply(data$parent_ids, as.character, character(1)))
  } else ggchord_stop(caller, ": each domain requires one explicit parent_id or parent_ids entry")
  ggchord_conversion_ids(out, caller)
}

ggchord_project_phase <- function(x, caller) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    if (any(!is.na(x) & !x %in% c("0", "1", "2", "."))) ggchord_stop(
      caller, ": GFF3 phase must be 0, 1, 2 or unknown '.'")
    x <- suppressWarnings(as.numeric(x))
  }
  x
}

# Expand an explicitly ordered nested CDS, never its intronic envelope.
ggchord_project_expand_cds <- function(input, caller) {
  pieces <- list()
  for (i in seq_len(nrow(input))) {
    row <- input[i, , drop = FALSE]
    if ("phase" %in% names(row)) row$phase <- ggchord_project_phase(row$phase, caller)
    nested <- if ("segments" %in% names(row)) row$segments[[1L]] else NULL
    if (is.null(nested)) {
      if (!"phase" %in% names(row)) {
        if (!"codon_start" %in% names(row)) ggchord_stop(caller, ": provide phase or codon_start explicitly")
        row$phase <- row$codon_start - 1
      }
      row$.cds_source_row <- i
      row$.cds_source_segment <- 1L
      pieces[[length(pieces) + 1L]] <- row
      next
    }
    if (!is.data.frame(nested) || !nrow(nested)) ggchord_stop(caller, ": CDS segments must be nonempty data frames")
    ggchord_require_columns(nested, c("start", "end"), caller)
    for (column in c("start", "end")) ggchord_project_integer(nested[[column]], column, caller)
    for (column in c("start", "end")) ggchord_project_integer(row[[column]], column, caller)
    if (row$start > row$end || any(nested$start > nested$end) ||
        any(nested$start < row$start | nested$end > row$end))
      ggchord_stop(caller, ": nested CDS segments must have ordered endpoints within the parent envelope")
    if ("segment_type" %in% names(nested) && any(nested$segment_type %in% "gap"))
      ggchord_stop(caller, ": CDS segments contain an explicit gap; supply only coding segments")
    if ("strand" %in% names(nested) && anyNA(nested$strand)) ggchord_stop(caller, ": nested CDS strands must be known")
    if ("strand" %in% names(nested) && any(as.character(nested$strand) != as.character(row$strand)))
      ggchord_stop(caller, ": nested CDS strand conflicts with parent")
    if ("segment_order" %in% names(nested)) {
      index <- ggchord_project_integer(nested$segment_order, "segment_order", caller)
      if (!setequal(index, seq_len(nrow(nested)))) ggchord_stop(caller, ": nested segment_order must identify each transcription-order segment once")
      source_index <- order(index)
    } else {
      source_index <- order(nested$start, nested$end)
      if (as.character(row$strand) == "-") source_index <- rev(source_index)
    }
    nested <- nested[source_index, , drop = FALSE]
    top_phase <- if ("phase" %in% names(row) && !is.na(row$phase)) row$phase else NULL
    if ("codon_start" %in% names(row) && !is.na(row$codon_start)) {
      start <- ggchord_project_integer(row$codon_start, "codon_start", caller)
      if (start > 3 || (!is.null(top_phase) && top_phase != start - 1))
        ggchord_stop(caller, ": codon_start must be 1/2/3 and agree with initial phase")
      top_phase <- start - 1
    }
    phase <- if ("phase" %in% names(nested)) ggchord_project_phase(nested$phase, caller) else rep(NA_real_, nrow(nested))
    if (!is.numeric(phase)) ggchord_stop(caller, ": nested phase must be numeric when supplied")
    if (is.na(phase[1L])) {
      if (is.null(top_phase)) ggchord_stop(caller, ": nested CDS needs explicit initial phase or codon_start")
      phase[1L] <- top_phase
    } else if (!is.null(top_phase) && phase[1L] != top_phase) ggchord_stop(caller, ": parent and initial segment phase disagree")
    width <- nested$end - nested$start + 1
    before <- c(0, utils::head(cumsum(width), -1L))
    expected <- (3 - ((before - phase[1L]) %% 3)) %% 3
    expected[1L] <- phase[1L]
    phase[is.na(phase)] <- expected[is.na(phase)]
    for (j in seq_len(nrow(nested))) {
      out <- row; out$segments <- NULL
      out$start <- nested$start[j]; out$end <- nested$end[j]; out$phase <- phase[j]
      out$segment_order <- j
      out$.cds_source_row <- i; out$.cds_source_segment <- source_index[j]
      pieces[[length(pieces) + 1L]] <- out
    }
  }
  ggchord_rbind_fill(pieces)
}

# Source intervals are explicit domain pieces, not a bounding envelope.
ggchord_project_intervals <- function(row, caller) {
  out <- if ("segments" %in% names(row) && !is.null(row$segments[[1L]])) row$segments[[1L]] else
    data.frame(start = row$start, end = row$end)
  if (!is.data.frame(out) || !nrow(out)) ggchord_stop(caller, ": segments must be a nonempty data.frame")
  ggchord_require_columns(out, c("start", "end"), caller)
  if ("segment_type" %in% names(out)) out <- out[!is.na(out$segment_type) & out$segment_type != "gap", , drop = FALSE]
  if (!nrow(out)) ggchord_stop(caller, ": feature contains no domain segments")
  for (column in c("start", "end")) ggchord_project_integer(out[[column]], column, caller)
  out <- out[order(out$start), , drop = FALSE]
  if (any(out$start > out$end) || any(out$start < row$start | out$end > row$end) ||
      (nrow(out) > 1L && any(out$start[-1L] <= utils::head(out$end, -1L))))
    ggchord_stop(caller, ": source segments must be disjoint and within the source interval")
  out
}

#' Build an explicit spliced CDS coordinate map
#'
#' Builds a mapping from observed protein residues to genomic coding bases.
#' CDS rows sharing `feature_id` form one coding model. Their genomic intervals
#' must be non-overlapping and have one sequence and known strand. Segments are
#' ordered in transcription direction, increasing on `+` and decreasing on `-`.
#' An explicit `segment_order` overrides genomic sorting, allowing a declared
#' circular-origin join. A `segments` list-column can declare `segment_order`
#' for transcription order, including reverse-strand and origin joins; otherwise
#' its pieces follow genomic strand order. Generic `segment_index` identifiers
#' are preserved as source metadata and do not define transcription order. Nested
#' phases may be supplied or derived from an explicit first phase or
#' `codon_start - 1`; no phase is guessed from the genomic envelope.
#' Trans-spliced models spanning different sequence IDs are rejected.
#'
#' GFF3 `phase` counts bases from a segment's biological 5-prime end to the
#' next codon boundary. Only the first phase is skipped; later phases validate
#' the accumulated coding length, retaining bases of codons split across exons.
#' `partial = "error"` rejects partial ends. Explicit `"trim"` permits a
#' nonzero first phase, partial-end flags and incomplete terminal codons, and
#' excludes bases without complete codons. Residue one then means the first
#' complete observed codon, not an inferred position in an unobserved protein.
#'
#' Optional `partial5`, `partial3`, `genetic_code`, and `protein_id` columns
#' must agree within each model. A nonzero initial phase requires partial5;
#' if absent it is recorded as partial5. A terminal remainder similarly records
#' partial3 unless a conflicting explicit flag is supplied. Genetic-code IDs
#' are retained as metadata: coordinate projection does not translate sequence,
#' identify stop codons, or infer a genetic code. `codon_count` includes any
#' source terminal stop codon. Frameshifts, ribosomal slippage and translation
#' exceptions or non-join location operators cannot be represented by a simple
#' triplet map and are rejected
#' when declared. No full protein length or identity is inferred.
#'
#' @param cds_data Data frame with `feature_id`, `accver`, ordered one-based
#'   inclusive `start`/`end`, `strand` (`+`/`-`) and GFF3 `phase` (0–2,
#'   numeric or literal text as returned by [read_gff3()]). `.` remains unknown;
#'   unknown initial phase requires explicit `codon_start` for nested CDS.
#'   If `type` is present, select CDS rows explicitly before building a map.
#'   Alternatively supply one row per CDS with a `segments` list-column and
#'   explicit initial `phase` or `codon_start` (1–3). Nested genomic pieces
#'   retain explicit `segment_order` when supplied and are never replaced by the
#'   parent's outer envelope.
#' @param genetic_code Optional positive integer code ID, scalar or named by
#'   CDS `feature_id`. Must agree with a supplied `genetic_code` column. Missing
#'   IDs remain unknown; they do not prevent coordinate-only projection.
#' @param partial `"error"` (default) or explicit `"trim"`.
#' @return A `ggchord_cds` list containing a `models` summary, transcription-
#'   ordered `segments` per feature, original `source_data` and provenance.
#'   Segment `cds_start`/`cds_end` identify complete-codon nucleotide positions;
#'   genomic coordinates and original row numbers remain available.
#' @seealso <https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md>
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' coding <- data.frame(feature_id = "cdsA", accver = "contig",
#'   start = c(100, 200), end = c(103, 204), strand = "+", phase = c(0, 2))
#' cds <- build_ggchord_cds(coding, genetic_code = 11)
#' cds$models
build_ggchord_cds <- function(cds_data, genetic_code = NULL,
                              partial = c("error", "trim")) {
  caller <- "build_ggchord_cds()"
  partial <- match.arg(partial)
  if (!is.data.frame(cds_data) || !nrow(cds_data)) ggchord_stop(caller, ": cds_data must be a nonempty data.frame")
  ggchord_require_columns(cds_data, c("feature_id", "accver", "start", "end", "strand"), caller)
  if ("type" %in% names(cds_data) && (anyNA(cds_data$type) || any(cds_data$type != "CDS")))
    ggchord_stop(caller, ": select CDS rows explicitly; noncoding feature envelopes are not coding models")
  if (anyNA(cds_data$strand) || any(!cds_data$strand %in% c("+", "-")))
    ggchord_stop(caller, ": CDS models require a known '+' or '-' strand")
  if ("location_operator" %in% names(cds_data) &&
      any(!is.na(cds_data$location_operator) & !cds_data$location_operator %in% c("join", "")))
    ggchord_stop(caller, ": only declared joined coding segments have a continuous observed-protein map")
  data <- ggchord_project_expand_cds(cds_data, caller)
  data$strand <- as.character(data$strand)
  data$feature_id <- ggchord_conversion_ids(data$feature_id, caller)
  data$accver <- ggchord_conversion_ids(data$accver, caller)
  for (name in c("start", "end")) data[[name]] <- ggchord_project_integer(data[[name]], name, caller)
  data$phase <- ggchord_project_integer(data$phase, "phase", caller, 0)
  if (any(data$phase > 2) || any(data$start > data$end) || anyNA(data$strand) ||
      any(!data$strand %in% c("+", "-"))) ggchord_stop(caller, ": ordered intervals, phase 0/1/2 and known strand are required")
  for (column in intersect(c("partial5", "partial3"), names(data))) {
    if (!is.logical(data[[column]]) || anyNA(data[[column]])) ggchord_stop(caller, ": ", column, " must be non-missing logical values")
  }
  for (column in intersect(c("transl_except", "translation_exception", "exception", "ribosomal_slippage", "frameshift"), names(data))) {
    declared <- vapply(data[[column]], function(value) {
      if (is.logical(value)) return(any(value, na.rm = TRUE))
      any(!is.na(value) & nzchar(as.character(value)) & !as.character(value) %in% c(".", "FALSE", "false", "0"))
    }, logical(1))
    if (any(declared)) ggchord_stop(caller, ": declared ", column, " cannot be represented by an ordinary triplet map")
  }
  ids <- unique(data$feature_id)
  if (!is.null(genetic_code)) {
    genetic_code <- ggchord_project_integer(genetic_code, "genetic_code", caller)
    if (is.null(names(genetic_code))) {
      if (length(genetic_code) != 1L) ggchord_stop(caller, ": genetic_code must be scalar or named for every CDS ID")
    } else if (anyNA(names(genetic_code)) || any(!nzchar(names(genetic_code))) ||
        anyDuplicated(names(genetic_code)) || any(!ids %in% names(genetic_code)))
      ggchord_stop(caller, ": genetic_code must be scalar or named for every CDS ID")
  }
  segments <- list(); models <- list()
  for (id in ids) {
    rows <- which(data$feature_id == id)
    x <- data[rows, , drop = FALSE]
    if (length(unique(x$accver)) != 1L || length(unique(x$strand)) != 1L)
      ggchord_stop(caller, ": CDS '", id, "' must have one sequence and strand")
    genomic <- order(x$start, x$end)
    if (length(rows) > 1L && any(x$start[genomic][-1L] <= utils::head(x$end[genomic], -1L)))
      ggchord_stop(caller, ": overlapping CDS segments for '", id, "'")
    if ("segment_order" %in% names(x) && !all(is.na(x$segment_order))) {
      explicit <- ggchord_project_integer(x$segment_order, "segment_order", caller)
      if (!setequal(explicit, seq_len(nrow(x)))) ggchord_stop(caller, ": segment_order must identify every CDS segment once")
      order <- order(explicit)
    } else order <- if (x$strand[1L] == "+") genomic else rev(genomic)
    x <- x[order, , drop = FALSE]; rows <- rows[order]
    width <- x$end - x$start + 1
    before <- c(0, utils::head(cumsum(width), -1L))
    initial <- x$phase[1L]
    expected <- (3 - ((before[-1L] - initial) %% 3)) %% 3
    if (length(rows) > 1L && any(x$phase[-1L] != expected)) ggchord_stop(caller,
      ": inconsistent phase continuity for CDS '", id, "'")
    total <- sum(width) - initial
    if (total < 3) ggchord_stop(caller, ": CDS '", id, "' has no complete observed codon")
    tail_bases <- total %% 3
    value <- function(column, fallback) {
      if (!column %in% names(x)) return(fallback)
      result <- unique(x[[column]])
      if (length(result) != 1L) ggchord_stop(caller, ": ", column, " must agree within CDS '", id, "'")
      result[[1L]]
    }
    partial5 <- value("partial5", initial != 0)
    partial3 <- value("partial3", tail_bases != 0)
    if ((initial != 0 && !partial5) || (tail_bases != 0 && !partial3))
      ggchord_stop(caller, ": explicit complete-end flag conflicts with phase or remainder for '", id, "'")
    if (partial == "error" && (partial5 || partial3 || initial != 0 || tail_bases != 0))
      ggchord_stop(caller, ": partial CDS '", id, "'; use partial = 'trim' with observed-protein coordinates")
    code <- value("genetic_code", NA_real_)
    if (length(code) != 1L || (!is.na(code) && (!is.numeric(code) || !is.finite(code) || code < 1 || code != floor(code))))
      ggchord_stop(caller, ": genetic_code must be a positive integer or unknown")
    if (!is.null(genetic_code)) {
      supplied <- if (is.null(names(genetic_code))) unname(genetic_code) else unname(genetic_code[id])
      if (!is.na(code) && code != supplied) ggchord_stop(caller, ": genetic_code conflicts with source for '", id, "'")
      code <- supplied
    }
    first <- pmax(before + 1, initial + 1)
    last <- pmin(before + width, initial + total - tail_bases)
    keep <- first <= last
    local_first <- first[keep] - before[keep] - 1
    local_last <- last[keep] - before[keep] - 1
    gs <- if (x$strand[1L] == "+") x$start[keep] + local_first else x$end[keep] - local_last
    ge <- if (x$strand[1L] == "+") x$start[keep] + local_last else x$end[keep] - local_first
    segments[[id]] <- data.frame(feature_id = id, accver = x$accver[1L], strand = x$strand[1L],
      segment_index = which(keep), start = gs, end = ge,
      cds_start = first[keep] - initial, cds_end = last[keep] - initial,
      phase = x$phase[keep], source_row = x$.cds_source_row[keep],
      source_segment_index = x$.cds_source_segment[keep], stringsAsFactors = FALSE)
    summary <- data.frame(feature_id = id, accver = x$accver[1L], strand = x$strand[1L],
      genetic_code = code, initial_phase = initial, partial5 = partial5, partial3 = partial3,
      terminal_bases = tail_bases, codon_count = (total - tail_bases) / 3,
      stringsAsFactors = FALSE)
    summary$source_rows <- I(list(unique(x$.cds_source_row)))
    if ("protein_id" %in% names(x)) summary$protein_id <- as.character(value("protein_id", NA_character_))
    models[[id]] <- summary
  }
  structure(list(models = do.call(rbind, models), segments = segments, source_data = cds_data, expanded_data = data,
    provenance = list(coordinate_system = "1-based-inclusive", partial = partial,
      source_import = attr(cds_data, "ggchord_import"), source_conversion = attr(cds_data, "ggchord_conversion"))),
    class = "ggchord_cds")
}

# Convert an inclusive observed coding-base interval to genomic pieces.
ggchord_project_genomic <- function(segment, lower, upper) {
  lo <- pmax(segment$cds_start, lower); hi <- pmin(segment$cds_end, upper)
  keep <- lo <= hi
  x <- segment[keep, , drop = FALSE]; lo <- lo[keep]; hi <- hi[keep]
  from <- lo - x$cds_start; to <- hi - x$cds_start
  first <- ifelse(x$strand == "+", x$start + from, x$end - from)
  last <- ifelse(x$strand == "+", x$start + to, x$end - to)
  data.frame(start = pmin(first, last), end = pmax(first, last),
    segment_index = seq_len(nrow(x)), segment_type = "standard",
    cds_start = lo, cds_end = hi, source_cds_row = x$source_row,
    source_cds_segment = x$source_segment_index,
    stringsAsFactors = FALSE)
}

#' Project protein domains or genomic features through a spliced CDS map
#'
#' Protein input uses one-based inclusive observed-residue coordinates. Genomic
#' output preserves exon gaps in `segments` and uses the parent's biological
#' strand; its outer envelope is not a claim that introns encode the domain.
#' The inverse projects only intersections with coding bases, reports excluded
#' bases, and uses the CDS feature ID as its protein coordinate namespace.
#' This namespace is not an inferred protein accession. Inputs may supply a
#' genomic `segments` list-column to describe disjoint genomic pieces.
#'
#' Each domain retains its explicit `feature_id`, single CDS parent, source
#' interval and source row. Source metadata survives with `source_` prefixes
#' on name conflicts. Multiple biological parents must be expanded explicitly
#' before projection. Unknown parent, intronic-only, out-of-range, reversed,
#' inconsistent-strand and ambiguous codon-boundary input produce errors.
#' `codon_boundary = "expand"` explicitly permits an inverse interval to cover
#' partial codons and reports this expansion; it does not modify source bases.
#'
#' @param feature_data Domain/feature data frame with `feature_id`, `parent_id`
#'   (or a single-valued `parent_ids` list-column), `start` and `end`.
#' @param cds A map from [build_ggchord_cds()].
#' @param from `"protein"` (default) or `"genome"`.
#' @param codon_boundary Inverse projection policy: `"error"` or `"expand"`.
#' @return Ordinary feature table with projected `accver`, `start`, `end`,
#'   `strand`, `directionality`, `segments`, parent IDs, genetic-code metadata,
#'   source coordinates and a `ggchord_projection` provenance attribute.
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' coding <- data.frame(feature_id = "cdsA", accver = "contig",
#'   start = c(100, 200), end = c(103, 204), strand = "+", phase = c(0, 2))
#' cds <- build_ggchord_cds(coding, genetic_code = 11)
#' domains <- data.frame(feature_id = "domainA", parent_id = "cdsA", start = 2, end = 3)
#' genomic <- project_ggchord_features(domains, cds)
#' genomic$segments[[1]]
#' project_ggchord_features(genomic, cds, from = "genome")
project_ggchord_features <- function(feature_data, cds, from = c("protein", "genome"),
                                     codon_boundary = c("error", "expand")) {
  caller <- "project_ggchord_features()"
  from <- match.arg(from); codon_boundary <- match.arg(codon_boundary)
  cds <- ggchord_project_model(cds, caller)
  if (!is.data.frame(feature_data)) ggchord_stop(caller, ": feature_data must be a data.frame")
  ggchord_require_columns(feature_data, c("feature_id", "start", "end"), caller)
  ids <- ggchord_conversion_ids(feature_data$feature_id, caller)
  parents <- ggchord_project_parent(feature_data, caller)
  for (name in c("start", "end")) ggchord_project_integer(feature_data[[name]], name, caller)
  if (any(feature_data$start > feature_data$end)) ggchord_stop(caller, ": feature intervals must be ordered")
  index <- match(parents, cds$models$feature_id)
  if (anyNA(index)) ggchord_stop(caller, ": unknown CDS parent: ", paste(unique(parents[is.na(index)]), collapse = ", "))
  result <- list()
  for (i in seq_len(nrow(feature_data))) {
    model <- cds$models[index[i], , drop = FALSE]; segment <- cds$segments[[parents[i]]]
    original <- feature_data[i, , drop = FALSE]
    excluded <- 0; expanded <- FALSE
    intervals <- ggchord_project_intervals(original, caller)
    if (from == "protein") {
      if (original$end > model$codon_count) ggchord_stop(caller, ": protein interval exceeds complete observed codons for '", parents[i], "'")
      pieces <- ggchord_rbind_fill(lapply(seq_len(nrow(intervals)), function(j)
        ggchord_project_genomic(segment, (intervals$start[j] - 1) * 3 + 1, intervals$end[j] * 3)))
      pieces$segment_index <- seq_len(nrow(pieces))
      accver <- model$accver; strand <- model$strand
      lower <- min(pieces$start); upper <- max(pieces$end)
    } else {
      source_model <- cds$expanded_data[cds$expanded_data$feature_id == parents[i], , drop = FALSE]
      if (original$start < min(source_model$start) || original$end > max(source_model$end))
        ggchord_stop(caller, ": genomic interval lies outside its CDS parent envelope")
      if ("accver" %in% names(original) && (is.na(original$accver) || as.character(original$accver) != model$accver))
        ggchord_stop(caller, ": genomic sequence does not match its CDS parent")
      if ("strand" %in% names(original) && !as.character(original$strand) %in% c(model$strand, ".", "*"))
        ggchord_stop(caller, ": genomic strand conflicts with its CDS parent")
      coding <- list()
      for (j in seq_len(nrow(intervals))) {
        lo <- pmax(segment$start, intervals$start[j]); hi <- pmin(segment$end, intervals$end[j])
        keep <- lo <= hi
        if (!any(keep)) next
        x <- segment[keep, , drop = FALSE]
        a <- ifelse(x$strand == "+", x$cds_start + lo[keep] - x$start, x$cds_start + x$end - hi[keep])
        b <- ifelse(x$strand == "+", x$cds_start + hi[keep] - x$start, x$cds_start + x$end - lo[keep])
        coding[[length(coding) + 1L]] <- data.frame(start = a, end = b)
      }
      if (!length(coding)) ggchord_stop(caller, ": genomic feature has no complete-codon CDS overlap")
      coding <- do.call(rbind, coding); coding <- coding[order(coding$start), , drop = FALSE]
      # Disjoint coding blocks can describe disjoint protein domains; preserve them.
      breaks <- c(TRUE, coding$start[-1L] > utils::head(coding$end, -1L) + 1)
      groups <- split(seq_len(nrow(coding)), cumsum(breaks))
      lower_nt <- vapply(groups, function(g) min(coding$start[g]), numeric(1))
      upper_nt <- vapply(groups, function(g) max(coding$end[g]), numeric(1))
      expanded <- any((lower_nt - 1) %% 3 != 0 | upper_nt %% 3 != 0)
      if (expanded && codon_boundary == "error") ggchord_stop(caller, ": genomic interval splits a codon; choose codon_boundary = 'expand' explicitly")
      lower_aa <- floor((lower_nt - 1) / 3) + 1; upper_aa <- ceiling(upper_nt / 3)
      # Codon expansion can make separate base pieces overlap one residue.
      aa_groups <- cumsum(c(TRUE, lower_aa[-1L] > utils::head(cummax(upper_aa), -1L) + 1))
      aa_indices <- split(seq_along(lower_aa), aa_groups)
      lower_aa <- vapply(aa_indices, function(g) min(lower_aa[g]), numeric(1))
      upper_aa <- vapply(aa_indices, function(g) max(upper_aa[g]), numeric(1))
      pieces <- data.frame(start = unname(lower_aa), end = unname(upper_aa), segment_index = seq_along(aa_indices),
        segment_type = "standard", stringsAsFactors = FALSE)
      excluded <- sum(intervals$end - intervals$start + 1) - sum(coding$end - coding$start + 1)
      accver <- parents[i]; strand <- "."; lower <- min(lower_aa); upper <- max(upper_aa)
    }
    row <- data.frame(feature_id = ids[i], parent_id = parents[i], accver = accver,
      start = lower, end = upper, strand = strand,
      directionality = if (strand == "+") "forward" else if (strand == "-") "reverse" else "none",
      genetic_code = model$genetic_code, source_start = original$start, source_end = original$end,
      source_coordinate_system = if (from == "protein") "observed-protein-1-based-inclusive" else "genome-1-based-inclusive",
      source_row = i, excluded_bases = excluded, codon_expanded = expanded,
      stringsAsFactors = FALSE)
    row$parent_ids <- I(list(parents[i])); row$segments <- I(list(pieces))
    extras <- original[setdiff(names(original), c("feature_id", "parent_id", "parent_ids"))]
    row <- ggchord_conversion_append(row, extras, prefix = "source_")$data
    result[[i]] <- row
  }
  if (!length(result)) {
    out <- data.frame(feature_id = character(), parent_id = character(), accver = character(),
      start = numeric(), end = numeric(), strand = character(), directionality = character())
    out$parent_ids <- I(list()); out$segments <- I(list())
  } else out <- ggchord_rbind_fill(result)
  rownames(out) <- NULL
  attr(out, "ggchord_projection") <- list(from = from, coordinate_system = "1-based-inclusive",
    models = cds$models, cds_provenance = cds$provenance, source_import = attr(feature_data, "ggchord_import"),
    source_conversion = attr(feature_data, "ggchord_conversion"), source_projection = attr(feature_data, "ggchord_projection"))
  out
}

#' Project explicitly ungapped protein links to genomic coding segments
#'
#' Splits paired protein intervals wherever either CDS crosses a segment
#' boundary, retaining query/subject roles and genomic endpoint directions.
#' This is a codon-coordinate correspondence, not a new nucleotide alignment.
#' Source protein identity and length remain prefixed metadata and are never
#' emitted as nucleotide `pident` or alignment `length`. General gapped
#' alignments require explicit per-column alignment information and cannot be
#' projected by this endpoint-only helper.
#'
#' @param link_data Data frame with `qfeature_id`, `sfeature_id` identifying CDS
#'   models and ordered inclusive protein `qstart`, `qend`, `sstart`, `send`.
#'   Query and subject residue spans must be equal. If present, `gapopen` must
#'   be zero and source `length` must equal this residue span.
#' @param cds A map from [build_ggchord_cds()].
#' @param ungapped Must explicitly be `TRUE`, confirming an ungapped residue
#'   correspondence rather than inferring it from equal endpoint spans.
#' @return A paired genomic interval table with CDS parent IDs, source row,
#'   piece number, genomic strands, `projected_bases`, source-prefixed input
#'   metadata and projection provenance. No alignment identity is computed.
#' @md
#' @seealso \code{\link{ggchord_data_contract}}
#' @export
#' @examples
#' coding <- data.frame(feature_id = c("a", "b"), accver = c("A", "B"),
#'   start = c(10, 100), end = c(18, 108), strand = c("+", "-"), phase = 0)
#' cds <- build_ggchord_cds(coding)
#' link <- data.frame(qfeature_id = "a", sfeature_id = "b",
#'   qstart = 1, qend = 2, sstart = 2, send = 3)
#' project_ggchord_links(link, cds, ungapped = TRUE)
project_ggchord_links <- function(link_data, cds, ungapped = FALSE) {
  caller <- "project_ggchord_links()"
  cds <- ggchord_project_model(cds, caller)
  if (!isTRUE(ungapped)) ggchord_stop(caller, ": set ungapped = TRUE only for a known ungapped residue correspondence")
  if (!is.data.frame(link_data)) ggchord_stop(caller, ": link_data must be a data.frame")
  required <- c("qfeature_id", "sfeature_id", "qstart", "qend", "sstart", "send")
  ggchord_require_columns(link_data, required, caller)
  for (name in c("qstart", "qend", "sstart", "send")) ggchord_project_integer(link_data[[name]], name, caller)
  qids <- ggchord_conversion_ids(link_data$qfeature_id, caller); sids <- ggchord_conversion_ids(link_data$sfeature_id, caller)
  qi <- match(qids, cds$models$feature_id); si <- match(sids, cds$models$feature_id)
  if (anyNA(qi) || anyNA(si)) ggchord_stop(caller, ": unknown CDS parent in qfeature_id/sfeature_id")
  width <- link_data$qend - link_data$qstart + 1
  if (any(width < 1 | link_data$send < link_data$sstart | width != link_data$send - link_data$sstart + 1) ||
      any(link_data$qend > cds$models$codon_count[qi] | link_data$send > cds$models$codon_count[si]))
    ggchord_stop(caller, ": ordered equal residue spans within observed CDS codons are required")
  for (name in intersect(c("length", "gapopen"), names(link_data))) {
    values <- link_data[[name]]
    if (!is.numeric(values) || any(!is.finite(values)) ||
        any(values != if (name == "length") width else 0)) ggchord_stop(caller, ": source ", name, " conflicts with ungapped residue spans")
  }
  output <- list()
  position <- function(segment, offset) {
    index <- which(segment$cds_start <= offset & segment$cds_end >= offset)
    x <- segment[index, , drop = FALSE]
    if (nrow(x) != 1L) ggchord_stop(caller, ": coding position is not represented uniquely")
    if (x$strand == "+") x$start + offset - x$cds_start else x$end - offset + x$cds_start
  }
  for (i in seq_len(nrow(link_data))) {
    q <- cds$segments[[qids[i]]]; s <- cds$segments[[sids[i]]]
    q0 <- (link_data$qstart[i] - 1) * 3 + 1; s0 <- (link_data$sstart[i] - 1) * 3 + 1
    span <- width[i] * 3
    cuts <- sort(unique(c(0, span, q$cds_end - q0 + 1, s$cds_end - s0 + 1)))
    cuts <- cuts[cuts >= 0 & cuts <= span]
    for (j in seq_len(length(cuts) - 1L)) {
      lower <- cuts[j]; upper <- cuts[j + 1L] - 1
      row <- data.frame(qaccver = cds$models$accver[qi[i]], saccver = cds$models$accver[si[i]],
        qstart = position(q, q0 + lower), qend = position(q, q0 + upper),
        sstart = position(s, s0 + lower), send = position(s, s0 + upper),
        qfeature_id = qids[i], sfeature_id = sids[i],
        qstrand = cds$models$strand[qi[i]], sstrand = cds$models$strand[si[i]],
        source_row = i, projection_piece = j, projected_bases = upper - lower + 1,
        orientation_known = TRUE, direction = if (cds$models$strand[qi[i]] == cds$models$strand[si[i]]) "same" else "reverse",
        stringsAsFactors = FALSE)
      extras <- link_data[i, setdiff(names(link_data), c("qfeature_id", "sfeature_id")), drop = FALSE]
      names(extras) <- paste0("source_", names(extras))
      output[[length(output) + 1L]] <- ggchord_conversion_append(row, extras, prefix = "source_")$data
    }
  }
  if (length(output)) out <- ggchord_rbind_fill(output) else out <- data.frame(
    qaccver = character(), saccver = character(), qstart = numeric(), qend = numeric(), sstart = numeric(), send = numeric())
  rownames(out) <- NULL
  attr(out, "ggchord_projection") <- list(from = "ungapped-protein-correspondence",
    coordinate_system = "1-based-inclusive", models = cds$models, cds_provenance = cds$provenance,
    source_import = attr(link_data, "ggchord_import"), source_conversion = attr(link_data, "ggchord_conversion"))
  out
}
