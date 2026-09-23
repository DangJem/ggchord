# Feature placement is resolved in the sequence-local frame before polygon
# generation.  The ggplot2 Position objects are deliberately identity at the
# Cartesian stage: their fields are consumed by ggchord's layout engine.

PositionGgchordFeature <- ggplot2::ggproto(
  "PositionGgchordFeature", ggplot2::PositionIdentity,
  ggchord_feature_position = TRUE,
  placement = "identity",
  offset = 0
)

#' Place features on opposite sides according to strand
#'
#' Positive signed offsets point away from the plot centre and negative
#' offsets point inward.  A scalar is expanded to `+ = offset` and
#' `- = -offset`; explicitly strand-named values are already signed and are
#' not flipped again.
#'
#' @param offset A finite scalar or a flexible per-sequence/per-strand
#'   specification.
#' @return A ggplot2 Position object.
#' @export
position_strand <- function(offset = 0.1) {
  ggplot2::ggproto(
    NULL, PositionGgchordFeature,
    placement = "strand", offset = offset
  )
}

#' Place both feature strands on one shared band
#'
#' @param offset Signed local-normal offset. Negative values place the shared
#'   band inside a circular sequence and positive values place it outside.
#' @return A ggplot2 Position object.
#' @export
position_plasmid <- function(offset = -0.1) {
  ggplot2::ggproto(
    NULL, PositionGgchordFeature,
    placement = "plasmid", offset = offset
  )
}

ggchord_as_feature_position <- function(position, caller = "feature layer") {
  if (is.character(position)) {
    if (length(position) != 1L || is.na(position) ||
        !position %in% c("identity", "strand", "plasmid")) {
      ggchord_stop(
        caller, ": position must be 'identity', 'strand', 'plasmid', ",
        "position_strand(), position_plasmid(), or position_feature_stack()"
      )
    }
    return(switch(
      position,
      identity = ggplot2::position_identity(),
      strand = position_strand(),
      plasmid = position_plasmid()
    ))
  }
  if (inherits(position, "PositionIdentity") ||
      isTRUE(position$ggchord_feature_position) ||
      isTRUE(position$ggchord_feature_stack)) return(position)
  ggchord_stop(
    caller, ": only identity, strand, plasmid, and feature-stack Positions ",
    "operate in the sequence-local frame"
  )
}

ggchord_is_identity_position <- function(position) {
  inherits(position, "PositionIdentity") &&
    !isTRUE(position$ggchord_feature_position) &&
    !isTRUE(position$ggchord_feature_stack)
}

ggchord_validate_offset_value <- function(x, caller) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    ggchord_stop(caller, ": offsets must be finite numeric scalars")
  }
  as.numeric(x)
}

# Resolve the same sequence/list forms accepted by process_gene_param(), while
# retaining whether a strand value was explicit (and therefore already signed).
ggchord_resolve_position_offset <- function(offset, seqs, placement,
                                             caller = "position") {
  n <- length(seqs)
  expand_scalar <- function(x) {
    x <- ggchord_validate_offset_value(x, caller)
    if (identical(placement, "strand")) c("+" = x, "-" = -x)
    else c("+" = x, "-" = x)
  }
  explicit_strands <- function(x) {
    if (is.null(names(x)) || !length(x) ||
        any(!names(x) %in% c("+", "-")) || anyDuplicated(names(x))) {
      ggchord_stop(caller, ": strand offsets must be uniquely named '+' and/or '-'")
    }
    if (identical(placement, "plasmid")) {
      ggchord_stop(caller, ": position_plasmid() accepts shared offsets, not strand-specific offsets")
    }
    out <- c("+" = 0, "-" = 0)
    for (st in names(x)) out[st] <- ggchord_validate_offset_value(x[[st]], caller)
    out
  }
  resolve_element <- function(x) {
    if (length(x) == 1L && is.null(names(x))) return(expand_scalar(x))
    explicit_strands(x)
  }

  if (!is.list(offset)) {
    nm <- names(offset)
    if (length(offset) == 1L && is.null(nm)) {
      val <- expand_scalar(offset)
      return(stats::setNames(rep(list(val), n), seqs))
    }
    if (!is.null(nm) && all(nm %in% c("+", "-"))) {
      val <- explicit_strands(offset)
      return(stats::setNames(rep(list(val), n), seqs))
    }
    if (!is.null(nm)) {
      unknown <- setdiff(nm, seqs)
      if (length(unknown)) ggchord_stop(caller, ": unknown accver: ", paste(unknown, collapse = ", "))
      out <- stats::setNames(rep(list(expand_scalar(0)), n), seqs)
      for (id in nm) out[[id]] <- expand_scalar(offset[[id]])
      return(out)
    }
    if (length(offset) != n) {
      ggchord_stop(caller, ": unnamed offsets must have length 1 or match the sequence count")
    }
    return(stats::setNames(lapply(offset, expand_scalar), seqs))
  }

  nm <- names(offset)
  if (length(offset) == 1L && (is.null(nm) || !nzchar(nm[1L]))) {
    val <- resolve_element(offset[[1L]])
    return(stats::setNames(rep(list(val), n), seqs))
  }
  if (!is.null(nm) && length(nm) && all(nm %in% seqs)) {
    out <- stats::setNames(rep(list(expand_scalar(0)), n), seqs)
    for (id in nm) out[[id]] <- resolve_element(offset[[id]])
    return(out)
  }
  idx_names <- as.character(seq_along(seqs))
  if (!is.null(nm) && length(nm) && all(nm %in% idx_names)) {
    out <- stats::setNames(rep(list(expand_scalar(0)), n), seqs)
    for (i in seq_along(nm)) out[[seqs[as.integer(nm[i])]]] <- resolve_element(offset[[i]])
    return(out)
  }
  if (is.null(nm) || all(!nzchar(nm))) {
    if (length(offset) != n) ggchord_stop(caller, ": unnamed offset lists must match the sequence count")
    return(stats::setNames(lapply(offset, resolve_element), seqs))
  }
  ggchord_stop(caller, ": invalid flexible offset specification")
}

ggchord_resolve_legacy_offset <- function(offset, seqs, name) {
  old <- process_gene_param(offset, seqs, name, 0, FALSE)
  lapply(old, function(x) c("+" = unname(x["+"]), "-" = -unname(x["-"])))
}

ggchord_feature_display_priority <- function(data, rows) {
  priority <- rep(0, length(rows))
  if ("display_priority" %in% names(data)) {
    value <- data$display_priority[rows]
    if (is.logical(value)) value <- as.numeric(value)
    value <- suppressWarnings(as.numeric(value))
    value[!is.finite(value)] <- 0
    priority <- value
  }
  if ("prioritized_display" %in% names(data)) {
    prioritized <- data$prioritized_display[rows]
    prioritized[is.na(prioritized)] <- FALSE
    priority[as.logical(prioritized)] <- Inf
  }
  priority
}

ggchord_allocate_feature_lanes <- function(data, base, direction, spacing,
                                            circular = FALSE, lengths = NULL) {
  lane <- integer(nrow(data))
  band <- paste(
    data$accver,
    format(base, digits = 15, scientific = FALSE, trim = TRUE),
    direction,
    sep = "\r"
  )
  groups <- split(seq_len(nrow(data)), factor(band, levels = unique(band)))
  for (idx in groups) {
    sid <- as.character(data$accver[idx[1L]])
    len <- if (!is.null(lengths)) unname(lengths[sid]) else NA_real_
    row_intervals <- lapply(idx, function(i) {
      a <- as.numeric(data$start[i]); b <- as.numeric(data$end[i])
      pad <- 0
      shape <- as.character(data$.feature_shape[i] %||%
        data$.feature_shape_raw[i] %||% "arrow")
      if (isTRUE(is.finite(len)) && identical(shape, "marker")) {
        glyph_width <- as.numeric(data$.feature_width[i] %||% 0)
        radius <- max(.1, abs(1 + base[i]))
        display_bp <- glyph_width * 3 / radius / (2 * pi) * len
        pad <- max(0, display_bp - abs(b - a)) / 2
      }
      if (isTRUE(circular) && is.finite(len) && a > b) {
        rbind(c(max(1, a - pad), len), c(1, min(len, b + pad)))
      } else matrix(c(max(1, min(a, b) - pad),
        if (is.finite(len)) min(len, max(a, b) + pad) else max(a, b) + pad),
        nrow = 1L)
    })
    # Only renderer-created rows from one multi-segment biological feature are
    # co-allocated. User-facing feature names/groups never collapse collisions:
    # local overlap itself is the grouping model.
    group_value <- if (".biological_source_row" %in% names(data)) {
      paste0(".feature.", data$.biological_source_row[idx])
    } else paste0(".row.", idx)
    entity_members <- split(seq_along(idx),
      factor(group_value, levels = unique(group_value)))
    intervals <- lapply(entity_members, function(members) {
      do.call(rbind, row_intervals[members])
    })
    entity_source <- vapply(entity_members, function(members) idx[members[1L]],
      integer(1))
    entity_priority <- vapply(entity_members, function(members) {
      max(ggchord_feature_display_priority(data, idx[members]))
    }, numeric(1))
    entity_span <- vapply(intervals, function(x) {
      sum(pmax(0, x[, 2] - x[, 1]))
    }, numeric(1))
    entity_start <- vapply(intervals, function(x) min(x[, 1]), numeric(1))
    # Longest/highest-priority first, then always reuse the nearest available
    # lane. Non-overlapping neighbours may therefore share a lane, but no
    # cassette-continuity preference may carry them past a parent interval.
    ord <- order(-entity_priority, -entity_span, entity_start, entity_source)
    lane_intervals <- list()
    overlaps <- function(a, occupied) {
      any(vapply(occupied, function(b) any(
        outer(seq_len(nrow(a)), seq_len(nrow(b)), Vectorize(function(i, j) {
          a[i, 1] <= b[j, 2] && b[j, 1] <= a[i, 2]
        }))
      ), logical(1)))
    }
    for (local in ord) {
      members <- entity_members[[local]]
      candidates <- seq_len(nrow(data) + 1L)
      chosen <- candidates[which(vapply(candidates, function(candidate) {
        candidate > length(lane_intervals) ||
          !overlaps(intervals[[local]], lane_intervals[[candidate]])
      }, logical(1)))[1L]]
      if (chosen > length(lane_intervals)) lane_intervals[[chosen]] <- list()
      lane_intervals[[chosen]][[length(lane_intervals[[chosen]]) + 1L]] <- intervals[[local]]
      lane[idx[members]] <- chosen - 1L
    }
  }
  lane
}

ggchord_feature_lane_offsets <- function(data, lane, base, direction,
                                          spacing, circular = FALSE,
                                          lengths = NULL, label_size = 2.5) {
  offsets <- numeric(nrow(data))
  label_text <- if ("anno" %in% names(data)) as.character(data$anno) else
    if ("label" %in% names(data)) as.character(data$label) else
      rep(NA_character_, nrow(data))
  label_metrics <- if (isTRUE(circular) && any(!is.na(label_text) &
      nzchar(label_text))) ggchord_text_boxes(data.frame(
        text = label_text, text_x = 0, text_y = 0,
        size = label_size
      ), units_per_inch = .35) else NULL
  band <- paste(
    data$accver,
    format(base, digits = 15, scientific = FALSE, trim = TRUE),
    direction, sep = "\r"
  )
  groups <- split(seq_len(nrow(data)), factor(band, levels = unique(band)))
  for (idx in groups) {
    local_lanes <- lane[idx]
    widths <- if (".feature_width" %in% names(data)) {
      as.numeric(data$.feature_width[idx])
    } else rep(0, length(idx))
    widths[!is.finite(widths) | widths < 0] <- 0
    lane_width <- vapply(seq.int(0L, max(local_lanes)), function(value) {
      used <- widths[local_lanes == value]
      if (length(used)) max(used) else 0
    }, numeric(1))
    centres <- numeric(length(lane_width))
    if (length(centres) > 1L) {
      sid <- as.character(data$accver[idx[1L]])
      sequence_length <- if (!is.null(lengths)) unname(lengths[sid]) else
        NA_real_
      lane_demand <- rep(1L, length(lane_width))
      if (isTRUE(circular) && !is.null(label_metrics) &&
          length(sequence_length) == 1L && is.finite(sequence_length) &&
          sequence_length > 0) {
        span <- ifelse(data$start[idx] <= data$end[idx],
          data$end[idx] - data$start[idx] + 1,
          sequence_length - data$start[idx] + data$end[idx] + 1)
        midpoint <- ((data$start[idx] - 1 + span / 2) %%
          sequence_length) / sequence_length
        radius <- pmax(.35, abs(1 + base[idx]))
        arc_width <- 2 * pi * radius * span / sequence_length
        needs_gutter <- is.finite(label_metrics$w[idx]) &
          nzchar(label_text[idx]) &
          label_metrics$w[idx] > arc_width * .85
        needs_gutter[is.na(needs_gutter)] <- FALSE
        for (lane_value in seq_len(length(lane_width) - 1L) - 1L) {
          members <- which(local_lanes == lane_value & needs_gutter)
          if (!length(members)) next
          half_span <- pmin(.25, label_metrics$w[idx[members]] /
            (4 * pi * radius[members]))
          overlap <- vapply(seq_along(members), function(j) {
            delta <- abs(midpoint[members] - midpoint[members[j]])
            delta <- pmin(delta, 1 - delta)
            sum(delta < half_span[j] + half_span + .004)
          }, integer(1L))
          lane_demand[lane_value + 1L] <- max(overlap)
        }
      }
      for (i in 2:length(centres)) {
        # Reserve only the corridor demanded by labels that cannot fit their
        # own arrows.  Fixed wide gutters made mostly empty outer bands push
        # the innermost annotations into the centre or even outside the map.
        label_gutter <- if (isTRUE(circular)) {
          base_gutter <- if (i == 2L) max(.125, spacing * 1.20) else
            max(.035, spacing * .40)
          step <- if (!is.null(label_metrics)) max(.05,
            max(label_metrics$h[idx], na.rm = TRUE) * 1.12) else .05
          min(if (i == 2L) .26 else .22,
            base_gutter + (lane_demand[i - 1L] - 1L) * step)
        } else if (i == 2L) {
          max(.125, spacing * 1.20)
        } else {
          max(.080, spacing * .85)
        }
        glyph_clearance <- lane_width[i - 1L] / 2 +
          lane_width[i] / 2 + label_gutter
        centres[i] <- centres[i - 1L] + max(spacing, glyph_clearance)
      }
    }
    offsets[idx] <- direction[idx] * centres[local_lanes + 1L]
  }
  offsets
}

ggchord_apply_feature_position <- function(data, position, seqs, lengths,
                                            circular = FALSE,
                                            label_size = 2.5,
                                            legacy_offset = NULL,
                                            legacy_name = "gene_offset") {
  if (is.null(data) || !nrow(data)) return(data)
  ggchord_require_columns(
    data, c("accver", "start", "end", "strand"), legacy_name
  )
  unknown <- setdiff(unique(as.character(data$accver)), seqs)
  if (length(unknown)) {
    ggchord_stop(legacy_name, ": unknown accver: ", paste(unknown, collapse = ", "))
  }
  if (anyNA(data$strand) || any(!as.character(data$strand) %in% c("+", "-"))) {
    ggchord_stop(legacy_name, ": strand must contain only '+' or '-'")
  }
  if (!is.numeric(data$start) || !is.numeric(data$end) ||
      any(!is.finite(data$start)) || any(!is.finite(data$end))) {
    ggchord_stop(legacy_name, ": start and end must be finite numeric values")
  }
  position <- position %||% ggplot2::position_identity()
  if (!is.null(legacy_offset)) {
    resolved <- ggchord_resolve_legacy_offset(legacy_offset, seqs, legacy_name)
    position_name <- "legacy_strand"
  } else if (isTRUE(position$ggchord_feature_stack)) {
    base_position <- position$base_position
    if (is.null(base_position)) {
      side <- position$side %||% "strand"
      resolved <- switch(
        side,
        strand = ggchord_resolve_position_offset(0.1, seqs, "strand"),
        outside = ggchord_resolve_position_offset(-0.1, seqs, "plasmid"),
        inside = ggchord_resolve_position_offset(0.1, seqs, "plasmid")
      )
      position_name <- "feature_stack_legacy"
    } else {
      resolved <- ggchord_resolve_position_offset(
        base_position$offset %||% 0, seqs,
        base_position$placement %||% "plasmid", "position_feature_stack(base_position)"
      )
      position_name <- paste0("feature_stack_", base_position$placement %||% "identity")
    }
  } else if (isTRUE(position$ggchord_feature_position)) {
    resolved <- ggchord_resolve_position_offset(
      position$offset, seqs, position$placement,
      paste0("position_", position$placement, "()")
    )
    position_name <- position$placement
  } else {
    resolved <- stats::setNames(rep(list(c("+" = 0, "-" = 0)), length(seqs)), seqs)
    position_name <- "identity"
  }

  out <- as.data.frame(data, stringsAsFactors = FALSE)
  base <- vapply(seq_len(nrow(out)), function(i) {
    unname(resolved[[as.character(out$accver[i])]][as.character(out$strand[i])])
  }, numeric(1))
  lane <- integer(nrow(out)); lane_offset <- numeric(nrow(out))
  if (isTRUE(position$ggchord_feature_stack)) {
    direction <- ifelse(base < 0, -1, 1)
    lane <- ggchord_allocate_feature_lanes(
      out, base, direction, position$spacing %||% 0.10,
      circular = circular, lengths = lengths
    )
    lane_offset <- ggchord_feature_lane_offsets(
      out, lane, base, direction, position$spacing %||% 0.10,
      circular = circular, lengths = lengths, label_size = label_size
    )
  }
  out$.position_name <- position_name
  out$.position_base_offset <- base
  out$.feature_stack_lane <- lane
  out$.position_lane_offset <- lane_offset
  out$.normal_offset <- base + lane_offset
  out
}
