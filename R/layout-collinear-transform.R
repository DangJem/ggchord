# Sequence display translations and direction evidence belong to the linear
# coordinate strategy. Biological positions remain in the source tables.

ggchord_validate_collinear_anchor <- function(anchor) {
  if (is.null(anchor)) return(invisible(NULL))
  if (is.numeric(anchor) && !is.object(anchor)) {
    ids <- names(anchor)
    values <- anchor
    if (is.null(ids) || length(ids) != length(anchor)) ggchord_stop(
      "coord_collinear(): numeric anchor positions must be named by accver"
    )
  } else if (is.data.frame(anchor)) {
    if (!"accver" %in% names(anchor)) ggchord_stop(
      "coord_collinear(): anchor data must contain accver"
    )
    point <- "position" %in% names(anchor)
    interval <- all(c("start", "end") %in% names(anchor))
    if ((point && any(c("start", "end") %in% names(anchor))) ||
        (!point && !interval)) {
      ggchord_stop(
        "coord_collinear(): anchor data must contain either position or start/end, not both"
      )
    }
    cols <- if (point) "position" else c("start", "end")
    if (any(!vapply(anchor[cols], is.numeric, logical(1L)))) ggchord_stop(
      "coord_collinear(): anchor coordinates must be numeric"
    )
    ids <- as.character(anchor$accver)
    values <- unlist(anchor[cols], use.names = FALSE)
  } else ggchord_stop(
    "coord_collinear(): anchor must be a named numeric vector or a data frame"
  )
  if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    ggchord_stop("coord_collinear(): each anchor accver must be non-empty and unique")
  }
  if (any(!is.finite(values))) ggchord_stop(
    "coord_collinear(): anchor coordinates must be finite"
  )
  invisible(NULL)
}

ggchord_collinear_anchor_ranges <- function(anchor, seqs, lens) {
  if (is.null(anchor)) return(data.frame(
    accver = seqs, start = 0, end = unname(lens[seqs]), stringsAsFactors = FALSE
  ))
  ggchord_validate_collinear_anchor(anchor)
  if (is.numeric(anchor)) {
    ranges <- data.frame(accver = names(anchor), start = as.numeric(anchor),
      end = as.numeric(anchor), stringsAsFactors = FALSE)
  } else {
    point <- "position" %in% names(anchor)
    ranges <- data.frame(accver = as.character(anchor$accver),
      start = if (point) anchor$position else anchor$start,
      end = if (point) anchor$position else anchor$end, stringsAsFactors = FALSE)
  }
  missing <- setdiff(seqs, ranges$accver)
  if (length(missing)) ggchord_stop(
    "coord_collinear(): every displayed sequence needs one anchor; missing: ",
    paste(missing, collapse = ", ")
  )
  ranges <- ranges[match(seqs, ranges$accver), , drop = FALSE]
  rownames(ranges) <- NULL
  outside <- ranges$start < 1 | ranges$end < 1 |
    ranges$start > unname(lens[seqs]) | ranges$end > unname(lens[seqs])
  if (any(outside)) ggchord_stop(
    "coord_collinear(): anchor coordinates must be within their sequence; invalid: ",
    paste(seqs[outside], collapse = ", ")
  )
  ranges
}

ggchord_collinear_direction_support <- function(seqs, orientation, ribbon_data,
                                                min_support) {
  out <- data.frame(accver = seqs, forward_support = 0, reverse_support = 0,
    confidence = NA_real_, status = "no_evidence", stringsAsFactors = FALSE)
  out$status[1L] <- "reference"
  data <- ribbon_data
  if (!is.null(data) && nrow(data)) {
    ggchord_require_columns(data,
      c("qaccver", "saccver", "qstart", "qend", "sstart", "send"),
      "coord_collinear(orientation = 'auto')")
    numeric <- c("qstart", "qend", "sstart", "send")
    if (any(!vapply(data[numeric], is.numeric, logical(1L))) ||
        any(!is.finite(as.matrix(data[numeric])))) ggchord_stop(
      "coord_collinear(): automatic orientation requires finite numeric link coordinates"
    )
    q <- as.character(data$qaccver); s <- as.character(data$saccver)
    q_span <- data$qend - data$qstart
    s_span <- data$send - data$sstart
    known <- rep(TRUE, nrow(data))
    if ("orientation_known" %in% names(data)) {
      known <- !is.na(data$orientation_known) &
        as.logical(data$orientation_known) %in% TRUE
    }
    for (strand in intersect(c("qstrand", "sstrand"), names(data))) {
      known <- known & as.character(data[[strand]]) %in% c("+", "-")
    }
    keep <- known & q %in% seqs & s %in% seqs & q != s & q_span != 0 & s_span != 0
    data <- data[keep, , drop = FALSE]
    q <- q[keep]; s <- s[keep]
    relative <- sign(q_span[keep] * s_span[keep])
    weight <- (abs(q_span[keep]) + abs(s_span[keep])) / 2 + 1
    if (length(seqs) > 1L) for (i in seq.int(2L, length(seqs))) {
      id <- seqs[i]
      earlier <- seqs[seq_len(i - 1L)]
      query <- q == id & s %in% earlier
      subject <- s == id & q %in% earlier
      hits <- query | subject
      if (!any(hits)) next
      peer <- ifelse(query[hits], s[hits], q[hits])
      votes <- relative[hits] * unname(orientation[peer])
      positive <- sum(weight[hits][votes > 0])
      negative <- sum(weight[hits][votes < 0])
      confidence <- abs(positive - negative) / (positive + negative)
      out$forward_support[i] <- positive
      out$reverse_support[i] <- negative
      out$confidence[i] <- confidence
      if (positive == negative) {
        out$status[i] <- "ambiguous"
      } else if (confidence < min_support) {
        out$status[i] <- "low_support"
      } else {
        orientation[id] <- if (positive > negative) 1 else -1
        out$status[i] <- "supported"
      }
    }
  }
  out$orientation <- unname(orientation[seqs])
  out$opposing_support <- ifelse(out$orientation == 1,
    out$reverse_support, out$forward_support)
  list(orientation = orientation, diagnostics = out)
}

ggchord_collinear_transform <- function(seqs, lens, orientation, ribbon_data, coord) {
  if (!length(seqs) || anyNA(seqs) || anyDuplicated(seqs)) ggchord_stop(
    "coord_collinear(): displayed sequence IDs must be non-missing and unique"
  )
  lens <- lens[seqs]
  orientation <- orientation[seqs]
  if (any(!is.finite(lens)) || any(lens <= 0) ||
      any(!is.finite(orientation)) || any(!orientation %in% c(-1, 1))) {
    ggchord_stop("coord_collinear(): sequence lengths and directions must be valid")
  }
  mode <- coord$collinear_orientation %||% "manual"
  diagnostics <- data.frame(accver = seqs, forward_support = 0,
    reverse_support = 0, confidence = NA_real_, status = "manual",
    orientation = unname(orientation), opposing_support = 0,
    stringsAsFactors = FALSE)
  if (identical(mode, "auto")) {
    directed <- ggchord_collinear_direction_support(seqs, orientation,
      ribbon_data, coord$collinear_min_support %||% .1)
    orientation <- directed$orientation
    diagnostics <- directed$diagnostics
  }
  ranges <- ggchord_collinear_anchor_ranges(coord$collinear_anchor, seqs, lens)
  align <- coord$collinear_align %||% "left"
  fraction <- switch(align, left = 0, center = .5, right = 1)
  # Left/right alignment is defined after reflection into display direction.
  displayed_start <- ifelse(orientation == 1, ranges$start, lens - ranges$start)
  displayed_end <- ifelse(orientation == 1, ranges$end, lens - ranges$end)
  low <- pmin(displayed_start, displayed_end)
  high <- pmax(displayed_start, displayed_end)
  anchor <- low + fraction * (high - low)
  scale <- 2 / max(lens)
  offset <- stats::setNames(-as.numeric(anchor) * scale, seqs)
  table <- data.frame(accver = seqs, orientation = unname(orientation),
    x_offset = unname(offset), units_per_bp = scale,
    xmin = unname(offset), xmax = unname(offset + lens * scale),
    anchor_start = ranges$start, anchor_end = ranges$end,
    anchor_position = ifelse(orientation == 1, anchor, lens - anchor),
    anchor_x = 0, align = align, stringsAsFactors = FALSE)
  list(orientation = orientation, offset = offset, table = table,
    diagnostics = diagnostics)
}
