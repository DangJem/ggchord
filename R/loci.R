# Explicit target selection produces independent genomic windows for focus.

#' Find genomic loci around explicitly identified features
#'
#' Selects targets by `feature_id` and expands their genomic bounds using
#' biological upstream/downstream directions. For a negative-strand target,
#' upstream extends towards larger genomic coordinates. Unknown or
#' bidirectional strands require equal upstream and downstream distances.
#' Windows are clipped to the known full sequence length; clipping is reported.
#'
#' Repeated rows with the same `feature_id` represent segments of one target
#' and must belong to the same sequence and strand. Their outer bounds define
#' the target window. Different targets always produce separate windows,
#' including overlapping targets on the same sequence. No intervening sequence
#' is removed, concatenated or locally reversed. Ordered one-based inclusive
#' genomic coordinates are required; origin-spanning circular features and
#' protein coordinates need an explicit conversion before this operation.
#'
#' The returned `locus_id` equals `target_id`. Pass the returned table as
#' `loci` to [focus_ggchord_data()] to create independent focused sequences
#' named by those IDs, including multiple windows from one source sequence.
#'
#' @param seq_data Sequence table containing unique `accver` and positive
#'   integer `length` values. Optional `genome` and `genome_id` columns are
#'   retained in the returned windows.
#' @param feature_data Feature table containing `feature_id`, `accver`,
#'   `start`, `end`, and optionally `strand`. Missing strand means unknown.
#' @param targets Non-empty character vector of distinct `feature_id` values.
#'   Selection is exact; no annotation text or accession versions are guessed.
#' @param upstream,downstream Single non-negative integer numbers of bases
#'   added before and after each target in its biological direction.
#' @return An ordinary data frame in target order with `target_id`, `locus_id`,
#'   `accver`, `start`, `end`, `strand`, unexpanded `target_start` and
#'   `target_end`, requested extensions, and `clipped_left`/`clipped_right`.
#'   The `target_source_rows` list column identifies feature input rows used
#'   for each target. Coordinates refer to the supplied sequence table; no
#'   implicit change from focused coordinates to an earlier source is made.
#' @md
#' @export
#' @examples
#' sequences <- data.frame(accver = "contig", length = 1000, genome = "isolate")
#' features <- data.frame(feature_id = c("geneA", "geneB"), accver = "contig",
#'   start = c(100, 700), end = c(200, 800), strand = c("+", "-"))
#' loci <- find_ggchord_loci(sequences, features, c("geneA", "geneB"),
#'   upstream = 50, downstream = 20)
#' loci
#' focused <- focus_ggchord_data(sequences, gene_data = features, loci = loci)
#' focused$seq_data
find_ggchord_loci <- function(seq_data, feature_data, targets,
                             upstream = 0, downstream = 0) {
  caller <- "find_ggchord_loci()"
  seq_data <- as_ggchord_seq(seq_data)
  if (!is.data.frame(feature_data)) ggchord_stop(
    caller, ": feature_data must be a data.frame with explicit feature_id values")
  ggchord_require_columns(feature_data, "feature_id", caller)
  ids <- as.character(feature_data$feature_id)
  if (anyNA(ids) || any(!nzchar(ids))) ggchord_stop(
    caller, ": feature_id must be non-missing and non-empty")
  features <- as_ggchord_feature(feature_data)
  if (!is.character(targets) || !length(targets) || anyNA(targets) ||
      any(!nzchar(targets)) || anyDuplicated(targets)) ggchord_stop(
    caller, ": targets must be distinct non-empty feature_id values")
  missing <- setdiff(targets, ids)
  if (length(missing)) ggchord_stop(caller, ": target feature_id not found: ",
    paste(missing, collapse = ", "))
  extension <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value != floor(value)) ggchord_stop(caller, ": ", name,
      " must be one non-negative integer number of bases")
    as.numeric(value)
  }
  upstream <- extension(upstream, "upstream")
  downstream <- extension(downstream, "downstream")
  sequence_rows <- match(features$accver, seq_data$accver)
  if (anyNA(sequence_rows)) ggchord_stop(caller,
    ": feature_data contains unknown sequence IDs")
  if (any(features$start > features$end) ||
      any(features$end > seq_data$length[sequence_rows])) ggchord_stop(caller,
    ": feature intervals must be ordered and within [1, sequence length]")

  windows <- lapply(targets, function(target) {
    rows <- which(ids == target)
    accver <- unique(features$accver[rows])
    strand <- unique(features$strand[rows])
    if (length(accver) != 1L || length(strand) != 1L) ggchord_stop(caller,
      ": ambiguous feature_id '", target,
      "'; its segments must have one sequence and strand")
    if (!strand %in% c("+", "-") && upstream != downstream) ggchord_stop(caller,
      ": target '", target,
      "' has no single known strand; use equal upstream and downstream")
    first <- min(features$start[rows])
    last <- max(features$end[rows])
    sequence_row <- match(accver, seq_data$accver)
    left <- if (strand == "-") downstream else upstream
    right <- if (strand == "-") upstream else downstream
    requested_start <- first - left
    requested_end <- last + right
    if (!is.finite(requested_start) || !is.finite(requested_end)) ggchord_stop(
      caller, ": expanded coordinates must remain finite")
    out <- data.frame(target_id = target, locus_id = target, accver = accver,
      start = max(1, requested_start),
      end = min(seq_data$length[sequence_row], requested_end), strand = strand,
      target_start = first, target_end = last,
      upstream = upstream, downstream = downstream,
      clipped_left = requested_start < 1,
      clipped_right = requested_end > seq_data$length[sequence_row],
      stringsAsFactors = FALSE)
    out$target_source_rows <- I(list(rows))
    for (column in intersect(c("genome", "genome_id"), names(seq_data))) {
      out[[column]] <- seq_data[[column]][sequence_row]
    }
    out
  })
  out <- do.call(rbind, windows)
  rownames(out) <- NULL
  out
}
