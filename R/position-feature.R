# Feature placement is resolved in the sequence-local frame before polygon
# generation.  The ggplot2 Position objects are deliberately identity at the
# Cartesian stage: their fields are consumed by ggchord's layout engine.

PositionGgchordFeature <- ggproto(
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
#' In [coord_collinear()], the positive local normal points upward and the
#' negative local normal downward, regardless of sequence orientation.
#'
#' @param offset A finite scalar or a flexible per-sequence/per-strand
#'   specification.
#' @return A ggplot2 Position object.
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
#' ggchord(data = seqs) + geom_seq() +
#'   geom_feature(data = genes, shape = "arrow",
#'     position = position_strand(list(A = c("+" = 2, "-" = -3), B = 1))) +
#'   coord_collinear()
position_strand <- function(offset = 1) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggproto(
    NULL, PositionGgchordFeature,
    placement = "strand", offset = ggchord_physical_mm(offset, "position_strand()")
  )
}

#' Place both feature strands on one shared band
#'
#' @param offset Signed local-normal offset. Negative values place the shared
#'   band inside a circular sequence and positive values place it outside. Use
#'   one value for all sequences, a named vector keyed by `accver`, or an
#'   existing equivalent list of scalar offsets. An unnamed vector must have
#'   one entry per sequence; omitted names receive zero.
#' @return A ggplot2 Position object.
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
#' ggchord(data = seqs) + geom_seq() +
#'   geom_feature(data = genes, position = position_offset(c(A = 2, B = -1))) +
#'   coord_collinear()
position_offset <- function(offset = 0) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggproto(
    NULL, PositionGgchordFeature,
    placement = "offset", offset = ggchord_physical_mm(offset, "position_offset()")
  )
}

ggchord_local_position <- function(position) {
  inherits(position, "Position") && (isTRUE(position$ggchord_feature_position) ||
    isTRUE(position$ggchord_feature_stack))
}

ggchord_is_data_identity <- function(position) {
  identical(position, "identity") ||
    (inherits(position, "PositionIdentity") && !ggchord_local_position(position))
}

ggchord_as_feature_position <- function(position, caller = "feature layer") {
  if (is.character(position)) {
    if (identical(position, "plasmid")) ggchord_stop(
      caller, ": replace position = 'plasmid' with position_offset(-0.1)")
    if (length(position) != 1L || is.na(position) ||
        !position %in% c("identity", "strand", "feature_stack")) {
      ggchord_stop(
        caller, ": use 'identity', 'strand', 'feature_stack', or a ",
        "sequence-local Position such as position_offset()"
      )
    }
    return(switch(
      position,
      identity = position_identity(),
      strand = position_strand(),
      feature_stack = position_feature_stack()
    ))
  }
  if (inherits(position, "Position")) {
    if (inherits(position, "PositionIdentity") || isTRUE(position$ggchord_feature_position) ||
        isTRUE(position$ggchord_feature_stack)) return(position)
  }
  ggchord_stop(
    caller, ": only identity, strand, offset and feature-stack Positions ",
    "operate in the sequence-local frame"
  )
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
  keys <- names(offset)
  if (anyDuplicated(keys) || (!is.null(keys) && anyNA(keys)))
    ggchord_stop(caller, ": repeated or invalid offset keys")
  expand_scalar <- function(x) {
    x <- ggchord_validate_offset_value(x, caller)
    if (identical(placement, "strand")) {
      if (x < 0) ggchord_stop(caller, ": symmetric strand offsets must be nonnegative")
      c("+" = x, "-" = -x, "." = 0)
    } else c("+" = x, "-" = x, "." = x)
  }
  explicit_strands <- function(x) {
    if (is.null(names(x)) || !length(x) ||
        any(!names(x) %in% c("+", "-", ".")) || anyDuplicated(names(x))) {
      ggchord_stop(caller, ": strand offsets must be uniquely named '+', '-' or '.'")
    }
    if (identical(placement, "offset")) {
      ggchord_stop(caller, ": for strand-specific offsets, use position_strand()")
    }
    out <- c("+" = 1, "-" = -1, "." = 0)
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
      return(setNames(rep(list(val), n), seqs))
    }
    if (!is.null(nm) && all(nm %in% c("+", "-", "."))) {
      val <- explicit_strands(offset)
      return(setNames(rep(list(val), n), seqs))
    }
    if (!is.null(nm)) {
      unknown <- setdiff(nm, seqs)
      if (length(unknown)) ggchord_stop(caller, ": unknown accver: ", paste(unknown, collapse = ", "))
      out <- setNames(rep(list(expand_scalar(if (identical(placement, "strand")) 1 else 0)), n), seqs)
      for (id in nm) out[[id]] <- expand_scalar(offset[[id]])
      return(out)
    }
    if (length(offset) != n) {
      ggchord_stop(caller, ": unnamed offsets must have length 1 or match the sequence count")
    }
    return(setNames(lapply(offset, expand_scalar), seqs))
  }

  nm <- names(offset)
  if (length(offset) == 1L && (is.null(nm) || !nzchar(nm[1L]))) {
    val <- resolve_element(offset[[1L]])
    return(setNames(rep(list(val), n), seqs))
  }
  if (!is.null(nm) && length(nm) && all(nm %in% seqs)) {
    out <- setNames(rep(list(expand_scalar(if (identical(placement, "strand")) 1 else 0)), n), seqs)
    for (id in nm) out[[id]] <- resolve_element(offset[[id]])
    return(out)
  }
  idx_names <- as.character(seq_along(seqs))
  if (!is.null(nm) && length(nm) && all(nm %in% idx_names)) {
    out <- setNames(rep(list(expand_scalar(if (identical(placement, "strand")) 1 else 0)), n), seqs)
    for (i in seq_along(nm)) out[[seqs[as.integer(nm[i])]]] <- resolve_element(offset[[i]])
    return(out)
  }
  if (is.null(nm) || all(!nzchar(nm))) {
    if (length(offset) != n) ggchord_stop(caller, ": unnamed offset lists must match the sequence count")
    return(setNames(lapply(offset, resolve_element), seqs))
  }
  ggchord_stop(caller, ": invalid flexible offset specification")
}

ggchord_lane_intervals <- function(data, row, length, circular, pad = 0) {
  parts <- data$.feature_runs[[row]]
  if (!NROW(parts)) parts <- data[row, c("start", "end"), drop = FALSE]
  intervals <- lapply(seq_len(nrow(parts)), function(i) {
    a <- as.numeric(parts$start[i])
    b <- as.numeric(parts$end[i])
    if (isTRUE(circular) && is.finite(length) && a > b) {
      rbind(c(max(1, a - pad), length), c(1, min(length, b + pad)))
    } else cbind(max(1, min(a, b) - pad),
      if (is.finite(length)) min(length, max(a, b) + pad) else max(a, b) + pad)
  })
  do.call(rbind, intervals)
}

ggchord_allocate_feature_lanes <- function(data, base, direction, spacing,
                                            circular = FALSE, lengths = NULL,
                                            lane_order = "source") {
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
      if (isTRUE(is.finite(len)) && identical(shape, "diamond")) {
        glyph_width <- as.numeric(data$.feature_width[i] %||% 0)
        radius <- max(.1, abs(1 + base[i]))
        display_bp <- glyph_width * 3 / radius / (2 * pi) * len
        pad <- max(0, display_bp - abs(b - a)) / 2
      }
      ggchord_lane_intervals(data, i, len, circular, pad)
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
    entity_start <- vapply(intervals, function(x) min(x[, 1]), numeric(1))
    entity_end <- vapply(intervals, function(x) max(x[, 2]), numeric(1))
    ids <- if ("feature_id" %in% names(data)) as.character(data$feature_id[entity_source]) else
      rep(NA_character_, length(entity_source))
    key <- data.frame(start = entity_start, end = entity_end, id = ids)
    if (anyDuplicated(key)) ggchord_stop(
      "position_feature_stack(): indistinguishable intervals need distinct feature_id values")
    ord <- order(entity_start, entity_end, ids, na.last = TRUE)
    if (lane_order == "priority") {
      priority <- suppressWarnings(as.numeric(
        as.character(data$display_priority %||% rep(0, nrow(data)))))
      priority[!is.finite(priority)] <- 0
      entity_priority <- vapply(entity_members, function(members)
        max(priority[idx[members]]), numeric(1))
      occupied_width <- vapply(intervals, function(parts)
        sum(parts[, 2L] - parts[, 1L] + 1), numeric(1))
      ord <- order(-entity_priority, -occupied_width, entity_start, entity_end,
        ids, na.last = TRUE)
    }
    lane_intervals <- list()
    overlaps <- function(parts, occupied) {
      for (i in seq_len(nrow(parts))) {
        if (any(parts[i, 1L] <= occupied[, 2L] & occupied[, 1L] <= parts[i, 2L])) return(TRUE)
      }
      FALSE
    }
    for (local in ord) {
      members <- entity_members[[local]]
      parts <- intervals[[local]]
      chosen <- 1L
      while (chosen <= length(lane_intervals) && overlaps(parts, lane_intervals[[chosen]]))
        chosen <- chosen + 1L
      lane_intervals[[chosen]] <- if (chosen > length(lane_intervals)) parts else
        rbind(lane_intervals[[chosen]], parts)
      lane[idx[members]] <- chosen - 1L
    }
  }
  lane
}

# Measure the inward extent of the displayed axis before placing features.
# Movable parallel labels use the same adjacent-to-backbone channel as the
# final axis resolver. Explicit positions keep their measured extent.
ggchord_circular_axis_inner_limit <- function(limit, ticks, plot, units_per_mm,
    rot_rad, rotation) {
  if (!nrow(ticks)) return(limit)
  edge <- limit
  tick_radius <- sqrt(ticks$x1^2 + ticks$y1^2)
  for (major in c(TRUE, FALSE)) {
    element <- ggchord_theme_element(plot,
      if (major) "ggchord.axis.ticks" else "ggchord.axis.minor.ticks")
    if (!inherits(element, "element_line")) next
    at <- ticks$is_major == major
    reach <- tapply(tick_radius[at] - (element@linewidth %||% 0) * units_per_mm / 2,
      ticks$accver[at], min)
    limit[names(reach)] <- pmin(limit[names(reach)], reach)
  }
  if (!inherits(ggchord_theme_element(plot, "ggchord.axis.text"), "element_text"))
    return(limit)
  at <- which(!is.na(ticks$label) & nzchar(ticks$label))
  if (!length(at)) return(limit)
  ticks <- ggchord_transform_axis_ticks(ticks[at, , drop = FALSE], TRUE, rot_rad, rotation)
  text <- transform(ticks, text = label, text_x = label_x, text_y = label_y,
    text_angle = label_angle, hjust = label_hjust, vjust = label_vjust)
  curved <- text$label_along_axis & !(text$.axis_plotmath %in% TRUE)
  units <- units_per_mm * 25.4
  # The distance to a rotated rectangle is the length of its shortest vector
  # to the origin. Taking the nearest corner would underestimate its extent.
  nearest <- function(boxes) {
    ca <- cos(boxes$angle); sa <- sin(boxes$angle)
    x <- abs(boxes$cx * ca + boxes$cy * sa)
    y <- abs(boxes$cx * sa - boxes$cy * ca)
    sqrt(pmax(x - boxes$ow / 2, 0)^2 + pmax(y - boxes$oh / 2, 0)^2)
  }
  depth <- numeric(nrow(text))
  if (any(!curved)) depth[!curved] <- nearest(ggchord_text_boxes(
    text[!curved, , drop = FALSE], units_per_inch = units))
  if (any(curved)) {
    arc <- ggchord_arc_text_layout(text[curved, , drop = FALSE], units_per_inch = units,
      respect_justification = TRUE, boxes_only = TRUE)
    rows <- which(curved)[arc$owners]
    inner <- tapply(nearest(arc$boxes), rows, min)
    outer <- tapply(sqrt(rowSums(ggchord_boxes_corners(arc$boxes)^2)),
      rep(rows, each = 4L), max)
    idx <- as.integer(names(inner))
    movable <- text$.axis_automatic_breaks[idx] %in% TRUE &
      text$.axis_allow_clearance[idx] %in% TRUE & !(text$.axis_offset_explicit[idx] %in% TRUE)
    depth[idx] <- inner
    depth[idx[movable]] <- pmax(inner[movable], edge[text$accver[idx[movable]]] -
      .85 * units_per_mm - (outer[movable] - inner[movable]))
  }
  reach <- tapply(depth[is.finite(depth) & depth > 0],
    text$accver[is.finite(depth) & depth > 0], min)
  limit[names(reach)] <- pmin(limit[names(reach)], reach)
  limit
}

# Reserve physical rows for nearby labels that cannot fit their own feature.
# Compare label footprints on the same track across the circular seam, and
# share any constrained radial budget before laying out the individual names.
ggchord_inner_track_gutters <- function(data, lanes, widths, spacing, length,
                                        radius, units_per_mm, labels) {
  gaps <- rep(spacing, length(widths) - 1L)
  if (!length(gaps) || !NROW(labels) || !is.finite(length) || length <= 0) return(gaps)
  labels <- labels[!is.na(labels$text) & nzchar(labels$text), , drop = FALSE]
  if (!nrow(labels)) return(gaps)
  intervals <- lapply(labels$source_row, function(i)
    ggchord_lane_intervals(data, i, length, circular = TRUE))
  spans <- vapply(intervals, function(parts) sum(parts[, 2L] - parts[, 1L] + 1), numeric(1))
  source <- data[labels$source_row, , drop = FALSE]
  center <- ggchord_feature_label_position(source, rep(length, nrow(source)),
    "middle", source$strand) / length
  boxes <- ggchord_text_boxes(transform(labels, text_x = 0, text_y = 0),
    units_per_inch = units_per_mm * 25.4)
  step <- (head(widths, -1L) + tail(widths, -1L)) / 2
  radii <- pmax(radius - c(0, cumsum(step + gaps)), .2 * radius)
  demand <- vapply(seq_along(gaps), function(i) {
    circumference <- 2 * pi * radii[i]
    rows <- which(lanes[labels$source_row] == i - 1L &
      boxes$w > .85 * spans / length * circumference)
    if (!length(rows)) return(0)
    # Count competing rows at each actual label anchor, scanning their
    # measured physical footprints and allowing a short tangential move.
    half_width <- pmin(.25, (boxes$w[rows] + labels$padding[rows] * units_per_mm) /
      (2 * circumference))
    starts <- center[rows] - half_width
    ends <- center[rows] + half_width
    all_starts <- sort(c(starts - 1, starts, starts + 1))
    all_ends <- sort(c(ends - 1, ends, ends + 1))
    competing <- findInterval(ends, all_starts, left.open = TRUE) -
      findInterval(starts, all_ends)
    pitch <- 1.15 * max(boxes$h[rows]) + max(.35, max(labels$padding[rows]) * .55) * units_per_mm
    max(competing) * pitch + .65 * units_per_mm
  }, numeric(1L))
  minimum_radius <- max(.25 * radius, tail(widths, 1L) / 2 + 5 * units_per_mm)
  available <- max(0, radius - minimum_radius - sum(step + gaps))
  extra <- pmax(0, pmin(.18 * radius, demand) - gaps)
  gaps + extra * if (sum(extra) > available) available / sum(extra) else 1
}

ggchord_feature_lane_offsets <- function(data, lane, base, direction, spacing,
    circular = FALSE, lengths = NULL, radii = NULL, units_per_mm = 1, labels = NULL) {
  offsets <- numeric(nrow(data))
  width <- data$.feature_clearance_width %||% data$.feature_width
  bands <- interaction(data$accver, base, direction, drop = TRUE)
  for (idx in split(seq_len(nrow(data)), bands)) {
    lanes <- lane[idx]
    widths <- vapply(seq_len(max(lanes) + 1L) - 1L, function(at) {
      max(width[idx[lanes == at]])
    }, numeric(1))
    neighbours <- (head(widths, -1L) + tail(widths, -1L)) / 2
    gaps <- spacing[idx[1L]]
    if (circular && direction[idx[1L]] < 0 && length(neighbours) && !is.null(radii)) {
      id <- as.character(data$accver[idx[1L]])
      local_labels <- if (NROW(labels)) {
        matched <- match(labels$source_row, idx)
        transform(labels[!is.na(matched), , drop = FALSE], source_row = matched[!is.na(matched)])
      }
      gaps <- ggchord_inner_track_gutters(data[idx, , drop = FALSE], lanes,
        widths, gaps, lengths[[id]], radii[[id]] + base[idx[1L]], units_per_mm,
        local_labels)
    }
    centres <- c(0, cumsum(neighbours + gaps))
    offsets[idx] <- direction[idx] * centres[lanes + 1L]
  }
  offsets
}

ggchord_apply_feature_position <- function(data, position, seqs, lengths,
                                            circular = FALSE,
                                            units_per_mm = ggchord_units_per_mm(),
                                            radii = NULL, inner_limit = NULL,
                                            label_size = NULL, label_requests = NULL,
                                            caller = "feature position") {
  if (is.null(data) || !nrow(data)) return(data)
  ggchord_require_columns(
    data, c("accver", "start", "end", "strand"), caller
  )
  unknown <- setdiff(unique(as.character(data$accver)), seqs)
  if (length(unknown)) {
    ggchord_stop(caller, ": unknown accver: ", paste(unknown, collapse = ", "))
  }
  if (anyNA(data$strand) || any(!as.character(data$strand) %in% c("+", "-"))) {
    ggchord_stop(caller, ": strand must contain only '+' or '-'")
  }
  if (!is.numeric(data$start) || !is.numeric(data$end) ||
      any(!is.finite(data$start)) || any(!is.finite(data$end))) {
    ggchord_stop(caller, ": start and end must be finite numeric values")
  }
  position <- position %||% position_identity()
  if (isTRUE(position$ggchord_feature_stack)) {
    base_position <- position$base_position
    placement <- base_position$placement %||% "identity"
    resolved <- ggchord_resolve_position_offset(
      base_position$offset %||% 0, seqs, placement, "position_feature_stack(base_position)"
    )
    position_name <- paste0("feature_stack_", placement)
  } else if (isTRUE(position$ggchord_feature_position)) {
    resolved <- ggchord_resolve_position_offset(
      position$offset, seqs, position$placement,
      paste0("position_", position$placement, "()")
    )
    position_name <- position$placement
  } else {
    resolved <- setNames(rep(list(c("+" = 0, "-" = 0, "." = 0)), length(seqs)), seqs)
    position_name <- "identity"
  }

  out <- as.data.frame(data, stringsAsFactors = FALSE)
  strand <- out$.feature_biological_strand %||% out$strand
  strand[strand == "+/-"] <- "."
  base <- units_per_mm * vapply(seq_len(nrow(out)), function(i) {
    unname(resolved[[as.character(out$accver[i])]][strand[i]])
  }, numeric(1))
  lane <- integer(nrow(out)); lane_offset <- numeric(nrow(out))
  if (isTRUE(position$ggchord_feature_stack)) {
    direction <- switch(position$direction,
      positive = rep(1, length(base)), negative = rep(-1, length(base)),
      strand = ifelse(strand == "-", -1, 1))
    gap <- ggchord_resolve_position_offset(position$spacing, seqs, "offset",
      "position_feature_stack(spacing)")
    spacing <- units_per_mm * vapply(out$accver, function(id) gap[[id]][[1L]], numeric(1))
    lane <- ggchord_allocate_feature_lanes(out, base, direction, spacing,
      circular = circular, lengths = lengths, lane_order = position$lane_order %||% "source")
    if (circular && isTRUE(position$automatic_base) && !is.null(inner_limit)) {
      inward <- direction < 0
      first <- inward & lane == 0L
      width <- out$.feature_clearance_width %||% out$.feature_width
      by_reference <- tapply(width[first], out$accver[first], max)
      ids <- as.character(out$accver[inward])
      base[inward] <- inner_limit[ids] - radii[ids] -
        by_reference[ids] / 2 - .6 * units_per_mm
    }
    if (circular && !is.null(label_size)) {
      label <- as.character(out$label %||% out$anno)
      if (length(label)) label_requests <- ggchord_rbind_fill(list(label_requests,
        data.frame(source_row = seq_len(nrow(out)), text = label, size = label_size,
          family = "", fontface = 1, lineheight = 1.2, padding = 1)))
    }
    lane_offset <- ggchord_feature_lane_offsets(out, lane, base, direction, spacing,
      circular = circular, lengths = lengths, radii = radii, units_per_mm = units_per_mm,
      labels = label_requests)
  }
  out$.position_name <- position_name
  out$.position_base_offset <- base
  out$.feature_stack_lane <- lane
  out$.position_lane_offset <- lane_offset
  out$.normal_offset <- base + lane_offset
  out
}
