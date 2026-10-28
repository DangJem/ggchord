# Explicit selection retains the row identity of the full resolved layer input.
ggchord_select_ribbons <- function(data, seqs, params, collinear = FALSE) {
  if (is.null(data)) return(list(data = NULL, report = data.frame()))
  n <- nrow(data)
  qi <- match(as.character(data$qaccver), seqs)
  si <- match(as.character(data$saccver), seqs)
  adjacent <- abs(qi - si) == 1L
  if (!collinear && length(seqs) > 1L) {
    adjacent <- adjacent | abs(qi - si) == length(seqs) - 1L
  }
  adjacent[is.na(adjacent)] <- FALSE
  valid <- !is.na(qi) & !is.na(si) & qi != si
  explicit_pairs <- params$pairs
  only_adjacent <- params$adjacent_only %||%
    (isTRUE(collinear) && is.null(explicit_pairs))
  keep <- valid
  reason <- ifelse(valid, NA_character_, "outside_sequence_set_or_self")
  if (only_adjacent) {
    keep <- keep & adjacent
    reason[valid & !adjacent] <- "non_adjacent"
  }
  if (!is.null(explicit_pairs)) {
    q <- as.character(explicit_pairs$qaccver)
    s <- as.character(explicit_pairs$saccver)
    if (anyNA(q) || anyNA(s) || any(!q %in% seqs | !s %in% seqs) || any(q == s)) {
      ggchord_stop("geom_link_ribbon(): pairs must name distinct displayed sequences")
    }
    if (isTRUE(params$adjacent_only)) {
      d <- abs(match(q, seqs) - match(s, seqs))
      allowed <- d == 1L | (!collinear & d == length(seqs) - 1L)
      if (any(!allowed)) ggchord_stop(
        "geom_link_ribbon(): non-adjacent pairs conflict with adjacent_only = TRUE; use FALSE"
      )
    }
    pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "\r")
    selected <- pair_key(as.character(data$qaccver), as.character(data$saccver)) %in%
      pair_key(q, s)
    reason[keep & !selected] <- "unselected_pair"
    keep <- keep & selected
  }
  report <- data.frame(source_row = seq_len(n), qaccver = data$qaccver,
    saccver = data$saccver, retained = keep, reason = reason,
    stringsAsFactors = FALSE)
  out <- data[keep, , drop = FALSE]
  out$.ggchord_source_row <- which(keep)
  list(data = out, report = report)
}

# The baseline clearance may shrink in a narrow channel; entity clearance may
# not. A fixed request that cannot fit is an actionable conflict.
ggchord_collinear_check_gap <- function(q, s, qgap, sgap, ref, caller) {
  distance <- abs(unname(ref$y[q] - ref$y[s]))
  if (distance == 0) return(invisible(NULL))
  if (any(!is.finite(c(qgap, sgap))) || any(c(qgap, sgap) < 0) ||
      max(qgap) + max(sgap) >= distance) ggchord_stop(
    caller, ": endpoint clearances do not fit between ", q, " and ", s,
    "; increase geom_seq(seq_gap = ...), reduce an explicit gap, or move feature bands"
  )
  invisible(NULL)
}

ggchord_match_feature_rows <- function(labels, features) {
  if (is.null(labels) || is.null(features)) return(NULL)
  if ("feature_id" %in% names(labels) && "feature_id" %in% names(features)) {
    a <- as.character(labels$feature_id); b <- as.character(features$feature_id)
    if (anyNA(a) || anyNA(b) || any(!nzchar(a)) || any(!nzchar(b)) ||
        anyDuplicated(a) || anyDuplicated(b)) ggchord_stop(
      "Linear label inheritance requires unique, non-missing feature_id values"
    )
    index <- match(a, b)
    if (anyNA(index)) return(NULL)
    for (column in c("accver", "start", "end", "strand")) {
      if (!all(as.character(labels[[column]]) == as.character(features[[column]][index]))) {
        ggchord_stop("Matching feature_id has conflicting ", column,
          "; keep label/source intervals consistent or use distinct feature IDs")
      }
    }
    return(index)
  }
  columns <- c("accver", "start", "end", "strand")
  if (nrow(labels) == nrow(features) && all(columns %in% names(labels)) &&
      all(columns %in% names(features)) &&
      identical(labels[columns], features[columns])) return(seq_len(nrow(labels)))
  NULL
}

#' Convert genomic positions to or from a linear ggchord layout
#'
#' Maps positions for ordinary ggplot2 annotation layers without modifying
#' source tables. Inverse mapping returns the chosen relative or source bp.
#' @param plot A ggchord plot using [coord_collinear()].
#' @param accver Sequence identifiers, scalar or one per position.
#' @param position Numeric genomic positions, or layout x when `inverse = TRUE`.
#' @param inverse Convert layout x back to bp.
#' @param source Use original source coordinates after focusing.
#' @return A data frame with accver, position, x and backbone y.
#' @export
#' @examples
#' p <- ggchord(data.frame(accver = "A", length = 1000)) +
#'   geom_seq() + coord_collinear()
#' map_ggchord_position(p, "A", 500)
map_ggchord_position <- function(plot, accver, position, inverse = FALSE, source = FALSE) {
  if (!inherits(plot, "ggchord") || !isTRUE(plot$coordinates$ggchord_collinear)) {
    ggchord_stop("map_ggchord_position(): use a ggchord plot with coord_collinear()")
  }
  if (!is.logical(inverse) || length(inverse) != 1L || is.na(inverse) ||
      !is.logical(source) || length(source) != 1L || is.na(source)) {
    ggchord_stop("map_ggchord_position(): inverse/source must be TRUE or FALSE")
  }
  if (!is.numeric(position) || any(!is.finite(position)) ||
      !length(accver) %in% c(1L, length(position))) {
    ggchord_stop("map_ggchord_position(): supply finite positions and matching accessions")
  }
  layout <- get_chord_layout(plot, build = TRUE)
  ref <- layout$sequence_reference
  ids <- rep_len(as.character(accver), length(position))
  if (anyNA(ids) || any(!ids %in% layout$seqs)) ggchord_stop(
    "map_ggchord_position(): unknown displayed accession")
  origin <- layout$sequence_layout$source_origin[match(ids, layout$sequence_layout$accver)]
  if (!length(origin)) origin <- rep(1, length(ids))
  if (inverse) {
    local <- (position - unname(ref$x_offset[ids])) / ref$x_scale
    local <- ifelse(ref$orientation[ids] == 1, local, ref$lens[ids] - local)
    bp <- local + if (source) origin - 1 else 0
    x <- position
  } else {
    local <- position - (if (source) origin - 1 else 0)
    x <- ggchord_collinear_x(local, ids, ref)
    bp <- position
  }
  data.frame(accver = ids, position = unname(bp), x = unname(x),
    y = unname(ref$y[ids]), stringsAsFactors = FALSE)
}
