# Explicit target selection produces independent genomic windows for focus.

#' Locate windows around known feature IDs
#'
#' Extensions follow biological strand. A target without a single
#' known strand can be extended symmetrically. The result preserves
#' individual overlapping windows and clips them to source boundaries.
#'
#' @param sequences Complete reference dictionary.
#' @param features Feature table with accver/start/end and feature_id.
#'   Repeated IDs must identify segments of one feature on one reference and
#'   strand.
#' @param targets Distinct, exact feature IDs.
#' @param upstream,downstream Integral upstream and downstream lengths in
#'   source coordinates.
#' @return A window table with locus_id, target_id, coordinates, source
#'   member rows, and clipping flags.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' loci <- find_ggchord_loci(seqs, genes, targets="g1", upstream=50, downstream=100)
#' loci
find_ggchord_loci <- function(sequences, features, targets, upstream = 0, downstream = 0) {
  seq_data <- sequences
  feature_data <- features
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
