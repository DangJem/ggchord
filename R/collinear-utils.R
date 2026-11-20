# Cross-row outline conventions belong to linear rows, not every derived Coord.
ggchord_collinear_row_styles <- function(coord) {
  isTRUE(coord$ggchord_collinear) &&
    !isTRUE(coord$ggchord_radial) && !isTRUE(coord$ggchord_sankey)
}

# Explicit selection retains the row identity of the full resolved layer input.
ggchord_select_ribbons <- function(data, seqs, params, collinear = FALSE, display_rows = NULL, cyclic = FALSE) {
  if (is.null(data)) return(list(data = NULL, report = data.frame()))
  n <- nrow(data)
  qi <- match(as.character(data$qaccver), seqs)
  si <- match(as.character(data$saccver), seqs)
  same_row <- rep(FALSE, n)
  if (!is.null(display_rows)) {
    qi <- unname(display_rows[as.character(data$qaccver)])
    si <- unname(display_rows[as.character(data$saccver)])
    same_row <- !is.na(qi) & !is.na(si) & qi == si
  }
  adjacent <- abs(qi - si) == 1L
  if ((!collinear || cyclic) && length(seqs) > 1L) {
    adjacent <- adjacent | abs(qi - si) == length(seqs) - 1L
  }
  adjacent[is.na(adjacent)] <- FALSE
  valid <- !is.na(qi) & !is.na(si) & qi != si
  explicit_pairs <- params$pairs
  only_adjacent <- params$adjacent_only %||%
    ((isTRUE(collinear) || cyclic) && is.null(explicit_pairs))
  keep <- valid
  reason <- ifelse(valid, NA_character_, "outside_sequence_set_or_self")
  reason[same_row] <- "same_display_row"
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
    if (!is.null(display_rows) && any(display_rows[q] == display_rows[s])) ggchord_stop(
      "geom_link_ribbon(): explicit pairs in one display row are unsupported; use distinct row_id values")
    if (isTRUE(params$adjacent_only)) {
      d <- if (is.null(display_rows)) abs(match(q, seqs) - match(s, seqs)) else
        abs(display_rows[q] - display_rows[s])
      allowed <- d == 1L | ((!collinear || cyclic) & d == length(seqs) - 1L)
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
  if (cyclic) {
    unknown <- ggchord_alignment_unknown(data, "coord_radial_synteny()")
    if (any(unknown)) warning("coord_radial_synteny(): omitted ", sum(unknown),
      " HSPs with unknown direction", call. = FALSE)
    reason[unknown] <- "unknown_direction"
    keep <- keep & !unknown
  }
  report <- data.frame(source_row = seq_len(n), qaccver = data$qaccver,
    saccver = data$saccver, retained = keep, reason = reason,
    stringsAsFactors = FALSE)
  out <- data[keep, , drop = FALSE]
  out$.ggchord_source_row <- which(keep)
  out$.ggchord_adjacent <- adjacent[keep]
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
    "; increase layout_sequences(gap = ...), reduce an explicit gap, or move feature bands"
  )
  invisible(NULL)
}

ggchord_match_feature_rows <- function(labels, features, same_source, strict = FALSE) {
  coordinates <- c("accver", "start", "end", "strand")
  if (is.null(labels) || is.null(features) ||
      !all(coordinates %in% intersect(names(labels), names(features)))) return(NULL)
  if (!nrow(labels)) return(integer())
  if (!nrow(features)) return(NULL)
  index <- NULL
  if ("feature_id" %in% names(labels) && "feature_id" %in% names(features)) {
    key <- function(d) {
      id <- as.character(d$feature_id)
      # Repeated accessions have distinct internal display keys at this stage.
      ifelse(is.na(id) | !nzchar(id), NA_character_,
        paste(nchar(d$accver), d$accver, id, sep = ":"))
    }
    a <- key(labels)
    b <- key(features)
    if (!anyNA(a)) {
      index <- match(a, b)
      if (anyNA(index)) {
        if (strict) ggchord_stop("No matching feature anchor for label feature_id: ",
          paste(unique(labels$feature_id[is.na(index)]), collapse = ", "))
        return(NULL)
      }
      if (anyDuplicated(b[b %in% a])) ggchord_stop(
        "Feature anchor has multiple rows for a requested feature_id in one displayed accession")
    }
  }
  if (is.null(index)) {
    if (!same_source || nrow(labels) != nrow(features)) return(NULL)
    index <- seq_len(nrow(labels))
  }
  fields <- c(coordinates, intersect(c(".source_file", ".source_accver", ".source_input",
    "source_sha256", "sequence_unit", ".feature_biological_strand"),
    intersect(names(features), names(labels))))
  for (field in fields) {
    a <- as.character(labels[[field]])
    b <- as.character(features[[field]][index])
    if (!identical(a, b)) {
      mismatch <- which(is.na(a) != is.na(b) | (!is.na(a) & !is.na(b) & a != b))
      if (strict) ggchord_stop("Feature anchor conflicts with ", field,
        " for label input row(s) ", paste(utils::head(mismatch, 5L), collapse = ", "))
      return(NULL)
    }
  }
  index
}
