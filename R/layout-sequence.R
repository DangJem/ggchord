#' Arrange reference paths and display windows
#'
#' @param order Complete display order of references or instances.
#' @param orientation Display orientation (1 or -1) by reference or instance.
#' @param length_mode Proportional or equal display lengths.
#' @param angular_gap Angular spacing in degrees for arc references.
#' @param row_gap,segment_gap,column_gap Physical distances in mm. Automatic
#'   Collinear row and Sankey row/column spacing uses the panel size and the
#'   space occupied by features, labels and leaders on facing and exterior
#'   sides of the paths. The automatic layout also leaves room beyond the
#'   sequence ends. Set `row_gap` or `column_gap` to specify that distance.
#' @param radius Relative radius for arc references.
#' @param inner_radius Relative inner radius for radial spokes. With `NULL`,
#'   retain a compact centre and enlarge it if nearby drawn bands need space.
#' @param curvature Signed backbone bend for Chord, Collinear, Sankey and
#'   Radial synteny layouts. Set a scalar, a vector named by source accession,
#'   or `list(instance = ...)` with a vector named by display instance.
#'   Zero draws straight paths and negative
#'   values reverse the bend. The default is 1 (the standard arc) in Chord
#'   and 0 in the other three layouts. In those layouts, a magnitude of 1
#'   displaces the midpoint by one tenth of the along-axis span in the local
#'   coordinate frame.
#' @param align Display alignment of linear paths: start, center or end.
#' @param anchor Source anchor positions or intervals for linear alignment.
#' @param windows An optional table with instance_id, accver, start and end.
#' @return A layout declaration to add to a ggchord plot.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' windows <- data.frame(accver = "A", start = 1, end = 500, window_id = "w1")
#'
#' ggchord(data = seqs) +
#'   geom_seq() + geom_gene(data = genes) +
#'   layout_sequences(
#'     order = c("B", "A"), orientation = c(A = 1, B = -1),
#'     row_gap = 8, align = "center"
#'   ) +
#'   coord_collinear()
#'
#' w <- data.frame(instance_id = c("A_left", "A_right"),
#'                 accver = "A", start = c(1, 501), end = c(500, 1000),
#'                 row_id = "A_split")
#' ggchord(data = one_seq) + geom_seq() +
#'   layout_sequences(windows = w, segment_gap = 3) + coord_collinear()
#'
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_gene(data = genes) +
#'   layout_sequences(curvature = c(A = 0.6, B = -0.4)) +
#'   coord_collinear()
layout_sequences <- function(order = NULL, orientation = 1, length_mode = "proportional",
    angular_gap = NULL, row_gap = NULL, segment_gap = NULL, column_gap = NULL,
    radius = NULL, inner_radius = NULL, curvature = NULL, align = NULL,
    anchor = NULL, windows = NULL) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!is.null(order)) ggchord_id_vector(order, "order")
  length_mode <- ggchord_match_choice(length_mode, c("proportional", "equal"), "length_mode")
  if (!is.null(align)) align <- ggchord_match_choice(align, c("start", "center", "end"), "align")
  if (!is.null(windows)) {
    ggchord_require_columns(windows, c("instance_id", "accver", "start", "end"), "layout_sequences(windows)")
    windows <- as.data.frame(windows)
    ggchord_id_vector(windows$instance_id, "windows instance_id")
    if (anyNA(windows$accver) || any(!nzchar(windows$accver)))
      ggchord_stop("windows must identify a source reference for every instance")
    for (key in c("start", "end")) {
      x <- windows[[key]]
      if (!is.numeric(x) || any(!is.finite(x) | x < 1 | x != trunc(x)))
        ggchord_stop("windows ", key, " must contain positive source positions")
    }
    if (any(windows$start > windows$end)) ggchord_stop("windows start must not exceed end")
    if (!is.null(windows$ring_id) && (!is.character(windows$ring_id) ||
        anyNA(windows$ring_id) || any(!nzchar(trimws(windows$ring_id)))))
      ggchord_stop("windows ring_id must contain nonempty ring identifiers")
  }
  structure(list(order = order, orientation = orientation, length_mode = length_mode,
    angular_gap = angular_gap, row_gap = row_gap, segment_gap = segment_gap,
    column_gap = column_gap, radius = radius, inner_radius = inner_radius,
    curvature = curvature, align = align, anchor = anchor, windows = windows),
    class = "ggchord_sequence_layout")
}

ggchord_id_vector <- function(x, field) {
  if (!is.character(x) || !length(x) || anyNA(x) ||
      any(!nzchar(trimws(x))) || anyDuplicated(x))
    ggchord_stop(field, " must contain distinct nonempty identifiers")
  invisible(x)
}

#' @export
ggplot_add.ggchord_sequence_layout <- function(object, plot, ...) {
  ggchord_local_console_errors()
  if (!inherits(plot, "ggchord")) ggchord_stop("Add a reference layout to a ggchord plot")
  if (!is.null(plot$ggchord$sequence_layout) || !is.null(plot$ggchord$dotplot_layout))
    message("A reference layout is already present and will be replaced.")
  plot$ggchord$sequence_layout <- plot$ggchord$dotplot_layout <- NULL
  field <- if (inherits(object, "ggchord_dotplot_layout")) "dotplot_layout" else "sequence_layout"
  plot$ggchord[[field]] <- object
  plot
}

#' Arrange the query and subject axes of a Dotplot
#'
#' @param x_order,y_order Unique reference IDs in query and subject order.
#' @param x_orientation,y_orientation Display orientations (1 or -1).
#' @param x_gap,y_gap Gaps in the corresponding source coordinate unit.
#' @return A Dotplot layout declaration.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' ggchord(data = seqs) +
#'   geom_alignment_segment(data = links) +
#'   layout_dotplot(x_order = "A", y_order = "B", y_orientation = -1) +
#'   coord_dotplot()
layout_dotplot <- function(x_order = NULL, y_order = NULL, x_orientation = 1,
    y_orientation = 1, x_gap = 0, y_gap = 0) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  if (!is.null(x_order)) ggchord_id_vector(x_order, "x_order")
  if (!is.null(y_order)) ggchord_id_vector(y_order, "y_order")
  for (gap in list(x_gap, y_gap)) if (!is.numeric(gap) || any(!is.finite(gap) | gap < 0)) ggchord_stop("Dotplot gaps must be nonnegative finite numbers")
  structure(list(x_order = x_order, y_order = y_order, x_orientation = x_orientation,
    y_orientation = y_orientation, x_gap = x_gap, y_gap = y_gap),
    class = c("ggchord_dotplot_layout", "ggchord_sequence_layout"))
}

# Resolve the supported flexible input forms before converting physical units.
ggchord_layout_values <- function(value, source, instance = source, field, default) {
  expand <- function(value, ids, keys, named_only = FALSE) {
    if (is.null(value)) return(rep(NA, length(ids)))
    if (anyNA(value) || (is.list(value) && !is.unit(value)))
      ggchord_stop(field, " must contain finite scalar settings")
    if (is.null(names(value))) {
      if (named_only || !length(value) %in% c(1L, length(ids)))
        ggchord_stop(field, " needs an unnamed scalar, exact-length vector or named settings")
      return(rep_len(value, length(ids)))
    }
    if (anyDuplicated(names(value)) || anyNA(names(value)) || any(!names(value) %in% keys))
      ggchord_stop(field, " has unknown or repeated configuration keys")
    value[match(ids, names(value))]
  }
  if (is.list(value) && !is.unit(value)) {
    if (is.null(names(value)) || anyDuplicated(names(value)) ||
        any(!names(value) %in% c("source", "instance")))
      ggchord_stop(field, " accepts list(source = ..., instance = ...)")
    out <- expand(value$source, source, unique(source), TRUE)
    local <- expand(value$instance, instance, instance, TRUE)
    out[!is.na(local)] <- local[!is.na(local)]
  } else out <- expand(value, source, unique(source))
  out[is.na(out)] <- default
  setNames(out, instance)
}

ggchord_layout_mm <- function(value, field) {
  if (is.unit(value)) {
    value <- vapply(seq_along(value), function(i) ggchord_unit_inches(value[i]) * 25.4,
      numeric(1))
  } else if (is.list(value)) {
    value <- lapply(value, ggchord_layout_mm, field = field)
  }
  if (!is.null(value) && !is.list(value) &&
      (!is.numeric(value) || any(!is.finite(value) | value < 0)))
    ggchord_stop(field, " must use nonnegative mm or absolute grid units")
  value
}

ggchord_configure_layout <- function(plot) {
  ggchord_validate_coord_contract(plot)
  if (isTRUE(plot$ggchord$layout_configured)) return(plot)
  coord <- plot$coordinates
  dotplot <- isTRUE(coord$ggchord_dotplot)
  layout <- plot$ggchord$sequence_layout
  if (dotplot) {
    if (!is.null(layout)) ggchord_stop("Use layout_dotplot() with coord_dotplot()")
    return(plot)
  }
  if (!is.null(plot$ggchord$dotplot_layout)) ggchord_stop("layout_dotplot() requires coord_dotplot()")
  layout <- layout %||% layout_sequences()
  source <- plot$ggchord$data$seq_data
  windows <- layout$windows
  mode <- if (isTRUE(coord$ggchord_radial)) "radial" else
    if (isTRUE(coord$ggchord_sankey)) "sankey" else
    if (isTRUE(coord$ggchord_collinear)) "linear" else
    if (isTRUE(coord$ggchord_circular)) "circular" else "chord"
  supported <- switch(mode, radial = c("inner_radius", "curvature"),
    sankey = c("row_gap", "column_gap", "curvature"),
    linear = c("row_gap", "segment_gap", "align", "anchor", "curvature"),
    circular = c("angular_gap", "radius"), chord = c("angular_gap", "radius", "curvature"))
  optional <- c("angular_gap", "row_gap", "segment_gap", "column_gap", "radius",
    "inner_radius", "curvature", "align", "anchor")
  bad <- setdiff(optional[!vapply(layout[optional], is.null, logical(1))], supported)
  if (length(bad)) ggchord_stop("Unsupported layout parameter for ", mode, ": ", bad[1L])
  if (!is.null(layout$segment_gap) && (is.null(windows$row_id) || !anyDuplicated(windows$row_id)))
    ggchord_stop("segment_gap applies to windows sharing a Collinear display row")
  if (mode != "linear" && anyDuplicated(windows$row_id))
    ggchord_stop("Multiple windows in one display row require coord_collinear()")
  if (!mode %in% c("circular", "chord") && !is.null(windows$ring_id))
    ggchord_stop("windows ring_id requires coord_circular() or coord_chord()")
  if (!is.null(windows)) {
    at <- match(windows$accver, source$accver)
    if (anyNA(at) || any(windows$end > source$length[at]))
      ggchord_stop("Layout windows must lie within a known complete reference")
  }
  keys <- windows$instance_id %||% source$accver
  order <- layout$order %||% sort(keys, method = "radix")
  if (!setequal(order, keys) || length(order) != length(keys))
    ggchord_stop("layout_sequences(order) must list every reference or display instance exactly once")
  if (!is.null(windows)) windows <- windows[match(order, keys), , drop = FALSE]
  ids <- windows$accver %||% order
  layout$order <- order
  values <- function(name, default) ggchord_layout_values(layout[[name]], ids, order, name, default)
  layout$orientation <- values("orientation", 1)
  if (!is.null(windows$orientation)) {
    known <- !is.na(windows$orientation)
    layout$orientation[known] <- windows$orientation[known]
  }
  if (any(!layout$orientation %in% c(-1, 1))) ggchord_stop("orientation must be 1 or -1")
  if (mode %in% c("chord", "circular")) {
    layout$gap <- values("angular_gap", if (mode == "chord") 2 else 0) / 360
    layout$radius <- values("radius", 1)
  }
  layout$curvature <- values("curvature", if (mode %in% c("chord", "circular")) 1 else 0)
  if (mode %in% c("linear", "sankey")) {
    layout$auto_row_gap <- is.null(layout$row_gap)
    layout$row_gap <- ggchord_layout_mm(layout$row_gap, "row_gap")
    layout$row_gap <- values("row_gap", 6)
  }
  updates <- list()
  if (mode == "radial") {
    inner <- layout$inner_radius %||% .12
    if (!is.numeric(inner) || length(inner) != 1L || !is.finite(inner) || inner < 0 || inner >= 1)
      ggchord_stop("inner_radius must be in [0, 1)")
    updates$radial_inner_radius <- 2 * inner / (1 - inner)
    updates$radial_inner_automatic <- is.null(layout$inner_radius)
  }
  if (mode == "linear") {
    updates$collinear_align <- c(start = "left", center = "center", end = "right")[[layout$align %||% "start"]]
    if (!is.null(layout$anchor)) {
      anchor <- layout$anchor
      if (!is.data.frame(anchor)) ggchord_stop("anchor must be a data frame")
      if ("pos" %in% names(anchor)) anchor$position <- anchor$pos
      if ("instance_id" %in% names(anchor)) anchor$accver <- anchor$instance_id
      ggchord_validate_collinear_anchor(anchor)
      updates$collinear_anchor <- anchor
    }
    if (!is.null(windows$row_id)) {
      updates$collinear_segments <- data.frame(accver = order,
        row_id = windows$row_id, orientation = unname(layout$orientation))
      updates$collinear_segment_gap_mm <- ggchord_layout_mm(layout$segment_gap %||% 2, "segment_gap")
    }
  }
  if (mode == "sankey" && !is.null(layout$column_gap)) {
    gap <- ggchord_layout_mm(layout$column_gap, "column_gap")
    if (length(gap) != 1L || gap <= 0) ggchord_stop("column_gap must be one positive distance")
    updates$sankey_column_gap_mm <- gap
  }
  layout$windows <- windows
  plot$coordinates <- ggproto(NULL, coord, !!!updates)
  plot$ggchord$sequence_layout <- layout
  plot$ggchord$layout_configured <- TRUE
  plot
}

ggchord_layout_units_per_mm <- function(axis) {
  measured <- ggchord_render_context()$units_per_mm
  if (length(measured) == 2L && all(is.finite(measured))) return(unname(measured[axis]))
  ggchord_device_units_per_inch(c(0, 2), c(0, 2)) / 25.4
}

ggchord_path_normal_axis <- function(coord) {
  if (!isTRUE(coord$ggchord_collinear) || isTRUE(coord$ggchord_radial)) return(NULL)
  vertical <- if (isTRUE(coord$ggchord_sankey)) identical(coord$sankey_axis, "y") else
    identical(coord$collinear_seq_axis, "vertical")
  if (vertical) "x" else "y"
}

# Use the measured panel after guides and margins. For a fixed aspect ratio,
# measuring the fitted panel would feed the chosen gap back into its own budget.
ggchord_spacing_extent_mm <- function(plot, axis) {
  size <- ggchord_render_context()$panel_mm[axis]
  if (is.null(plot$coordinates$ratio) && length(size) && is.finite(size) && size > 0)
    return(unname(size))
  size <- ggchord_render_dimensions()[axis]
  if (!is.finite(size)) size <- 7 * 25.4
  facets <- plot$ggchord$source_facet_setup$panels
  if (NROW(facets)) size <- size / max(if (axis == 1L) facets$COL else facets$ROW)
  unname(size)
}

# Divide the free span after allowing for the drawn geometry on facing and
# exterior sides. Each gap gets the same pitch so references stay aligned.
ggchord_auto_spacing_mm <- function(plot, axis, gaps, share, maximum, clearance = 0) {
  panel <- ggchord_spacing_extent_mm(plot, axis)
  footprint <- ggchord_render_context()$reference_spacing[[axis]]
  inside <- max(clearance, footprint$inside %||% 0)
  outside <- footprint$outside %||% 0
  gaps <- max(1L, gaps)
  free <- max(0, panel - outside - gaps * inside)
  min(maximum * panel / gaps, inside + share * free / gaps)
}

ggchord_auto_row_gap <- function(plot) {
  ids <- plot$ggchord$sequence_layout$order
  if (isTRUE(plot$coordinates$ggchord_sankey)) {
    axis <- if (identical(plot$coordinates$sankey_axis, "x")) 1L else 2L
    rows <- sum(plot$ggchord$display_instances$side == "left")
    # Most of this axis is needed for the lengths of the stacked sequences.
    gap <- ggchord_auto_spacing_mm(plot, axis, rows - 1L, .065, .30, clearance = 4)
  } else {
    segments <- plot$coordinates$collinear_segments
    rows <- if (is.null(segments)) length(ids) else
      length(unique(segments$row_id[segments$accver %in% ids]))
    axis <- if (identical(plot$coordinates$collinear_seq_axis, "vertical")) 1L else 2L
    share <- .32 + .06 * min(4L, max(0L, rows - 2L))
    gap <- ggchord_auto_spacing_mm(plot, axis, rows - 1L, share, .72)
  }
  setNames(rep(gap, length(ids)), ids)
}

# Measure how far a reference's own polygons, text boxes and leader ends
# extend beyond each end or side of its backbone, in final projected axes.
ggchord_reference_spacing <- function(layout) {
  paths <- layout$seq_arcs
  metric <- ggchord_annotation_metric(layout)
  pieces <- unlist(layout$layer_geometry, recursive = FALSE)
  pieces <- Filter(function(d) NROW(d) && all(c("accver", "x", "y") %in% names(d)), pieces)
  lapply(c("x", "y"), function(axis) {
    path <- vapply(paths, function(d) range(d[[axis]]), numeric(2))
    low <- high <- setNames(numeric(length(paths)), names(paths))
    for (d in pieces) {
      columns <- intersect(c(axis, paste0(axis, c("0", "1", "end"))), names(d))
      lower <- c(columns, intersect(paste0("bbox_", axis, "min"), names(d)))
      upper <- c(columns, intersect(paste0("bbox_", axis, "max"), names(d)))
      lo <- do.call(pmin, c(unname(as.list(d[lower])), list(na.rm = TRUE)))
      hi <- do.call(pmax, c(unname(as.list(d[upper])), list(na.rm = TRUE)))
      by_path <- split(seq_len(nrow(d)), as.character(d$accver))
      for (id in intersect(names(by_path), names(paths))) {
        at <- by_path[[id]]
        low[id] <- max(low[id], path[1L, id] - lo[at], na.rm = TRUE)
        high[id] <- max(high[id], hi[at] - path[2L, id], na.rm = TRUE)
      }
    }
    # Mirrored paths and multiple windows on one row share the same slot.
    group <- match(path[1L, ], sort(unique(path[1L, ])))
    negative <- as.numeric(tapply(low, group, max))
    positive <- as.numeric(tapply(high, group, max))
    inside <- max(c(0, head(positive, -1L) + tail(negative, -1L)))
    list(inside = unname(inside / metric[axis]),
      outside = unname((negative[1L] + tail(positive, 1L)) / metric[axis]),
      limits = range(path[1L, ] - low, path[2L, ] + high))
  })
}

ggchord_row_gap <- function(plot, layout, seqs) {
  coord <- plot$coordinates
  if (is.null(layout$row_gap)) return(rep(0, length(seqs)))
  axis <- if (identical(coord$collinear_seq_axis, "vertical")) 1L else 2L
  if (isTRUE(coord$ggchord_sankey)) axis <- if (identical(coord$sankey_axis, "x")) 1L else 2L
  gap <- layout$row_gap[seqs]
  # Parallel Collinear references reserve space for both neighbouring sides.
  if (!isTRUE(coord$ggchord_sankey)) {
    tracks <- Filter(function(layer) identical(layer$ggchord_type, "track"), plot$layers)
    ids <- unique(vapply(tracks, function(layer) layer$ggchord_params$track_id, character(1)))
    reserved <- 0
    for (id in ids) {
      band <- plot$ggchord$track_bands[[id]]
      if (is.null(band)) next
      dimensions <- ggchord_track_band_dimensions(band, seqs, plot)
      reserved <- reserved + dimensions$width$mm + dimensions$gap$mm
    }
    gap <- if (isTRUE(layout$auto_row_gap)) pmax(gap, reserved) else gap + reserved
  }
  gap * ggchord_layout_units_per_mm(axis) / 2
}

# Sequence allocation, references and track geometry.
ggchord_layout_sequence_step <- function(
    axisGap, axisMajLen, axisMinLen, circular, clockwise, geometry_cache, labelOffset, lens, orientation,
    rotation, seq_curvature, seq_gap, seqRadius, seqs, length_mode, ring) {
  n <- length(seqs)

  # ====================================================================
  # Step 1: compute angle allocation
  # ====================================================================
  total_circ <- 2 * pi

  weights <- if (identical(length_mode, "equal")) rep(1, length(lens)) else lens
  groups <- if (is.null(ring)) list(seq_len(n)) else split(seq_len(n), ring)
  starts <- ends <- setNames(numeric(n), seqs)
  for (members in groups) {
    gap <- seq_gap[members]
    available <- 1 - sum(gap)
    if (available <= 0)
      ggchord_stop("Reference gaps must leave drawable space in each ring")
    weight <- weights[members]
    span <- weight / sum(weight) * total_circ * available
    offsets <- c(0, head(cumsum(span + total_circ * gap), -1L))
    starts[members] <- offsets
    ends[members] <- offsets + span
  }
  if (clockwise) {
    original <- starts
    starts <- -ends
    ends <- -original
  }

  # ====================================================================
  # Step 2: convert to radians
  # ====================================================================
  rot_rad <- rotation * pi / 180

  # ====================================================================
  # Step 3: generate reference paths and map_to_curve
  # ====================================================================
  nSeg <- 500
  nRef <- 2000

  # High-resolution reference path for each sequence
  reference_key <- list(
    seqs = seqs, starts = starts, ends = ends,
    radius = seqRadius, curvature = seq_curvature, circular = circular
  )
  if (!is.null(geometry_cache) &&
      identical(geometry_cache$reference_key, reference_key)) {
    seq_refs <- geometry_cache$seq_refs
  } else {
    seq_refs <- lapply(seqs, function(id) {
      path <- generate_curvature_path(
        starts[id], ends[id], seqRadius[id], seq_curvature[id], n_points = nRef
      )
      angles <- seq(starts[id], ends[id], length.out = nRef)
      list(path = path, angles = angles, r0 = seqRadius[id],
        circle = isTRUE(circular) && seq_curvature[id] == 1)
    })
    names(seq_refs) <- seqs
    if (!is.null(geometry_cache)) {
      geometry_cache$reference_key <- reference_key
      geometry_cache$seq_refs <- seq_refs
    }
  }

  gene_track_radius <- function(gene, sid, strand, angle) {
    signed_outward <- if (".normal_offset" %in% names(gene)) {
      as.numeric(gene[[".normal_offset"]])
    } else 0
    if (!is.finite(signed_outward) || signed_outward == 0) {
      return(unname(seqRadius[sid]))
    }
    ref <- seq_refs[[sid]]
    base <- map_to_curve(angle, seqRadius[sid], ref)
    plus <- map_to_curve(angle, seqRadius[sid] + 1e-5, ref)
    left_normal <- plus - base
    # Convert the stable signed-outward Position contract to the curve
    # mapper's left-normal radius convention at this exact genomic location.
    # This remains valid for signed and non-unit sequence curvature.
    radius_sign <- if (sum(left_normal * base) >= 0) 1 else -1
    unname(seqRadius[sid]) + radius_sign * signed_outward
  }

  # ====================================================================
  # Step 4: generate sequence arcs (outer layer)
  # ====================================================================
  seq_arcs <- setNames(lapply(seqs, function(id) {
    path_data <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id], seq_curvature[id], nSeg
    )
    path_data$accver <- id
    if (orientation[id] == -1) {
      path_data <- path_data[nrow(path_data):1, ]
    }
    attr(path_data, "ggchord_path_direction") <- unname(orientation[id])
    path_data
  }), seqs)

  axis_reference_x <- unlist(lapply(seq_arcs, `[[`, "x"), use.names = FALSE)
  axis_reference_y <- unlist(lapply(seq_arcs, `[[`, "y"), use.names = FALSE)
  axis_units_per_inch <- ggchord_device_units_per_inch(
    axis_reference_x, axis_reference_y
  )
  axisGap <- axisGap * axis_units_per_inch
  axisMajLen <- axisMajLen * axis_units_per_inch
  axisMinLen <- axisMinLen * axis_units_per_inch
  labelOffset <- labelOffset * axis_units_per_inch

  # ====================================================================
  list(n = n, starts = starts, ends = ends, rot_rad = rot_rad, nSeg = nSeg, seq_refs = seq_refs,
    gene_track_radius = gene_track_radius, seq_arcs = seq_arcs,
    axis_units_per_inch = axis_units_per_inch, axisGap = axisGap, axisMajLen = axisMajLen,
    axisMinLen = axisMinLen, labelOffset = labelOffset)
}

#' Suggest reference order and orientation from link evidence
#'
#' @param sequences A complete reference table.
#' @param links An interval relationship table supplying direction evidence.
#' @param order Optional complete display order.
#' @param reference Optional source accession to keep in the positive orientation.
#' @param weight Use aligned length, observation count, or a numeric input column.
#' @param min_support Minimum relative imbalance required to choose a direction.
#' @return A list of layout arguments, evidence and diagnostics.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' suggestion <- suggest_sequence_layout(seqs, links, weight = "count")
#' suggestion$evidence
#' ggchord(data = seqs) + geom_seq() + geom_link_ribbon(data = links) +
#'   do.call(layout_sequences, suggestion$layout) + coord_collinear()
suggest_sequence_layout <- function(sequences, links, order = NULL, reference = NULL,
    weight = "aligned_length", min_support = 0.1) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  caller <- "suggest_sequence_layout()"
  sequences <- ggchord_reference_table(sequences, caller)
  links <- as_ggchord_link(links)
  ids <- sort(sequences$accver, method = "radix")
  if (!is.null(order)) {
    ggchord_id_vector(order, "order")
    if (!setequal(order, ids)) ggchord_stop(caller, ": order must list all references")
  }
  if (!is.null(reference) && (length(reference) != 1L || !reference %in% ids))
    ggchord_stop(caller, ": reference must select one known accession")
  if (!is.numeric(min_support) || length(min_support) != 1L ||
      !is.finite(min_support) || min_support < 0 || min_support > 1)
    ggchord_stop(caller, ": min_support must be in [0, 1]")
  if (!is.character(weight) || length(weight) != 1L || is.na(weight))
    ggchord_stop(caller, ": weight must identify one weighting rule or column")
  ggchord_check_reference_ids(links, c("qaccver", "saccver"), ids, caller)
  unit <- setNames(sequences$sequence_unit, sequences$accver)
  q <- links$qaccver
  s <- links$saccver
  aligned <- pmin(links$qend - links$qstart, links$send - links$sstart) + 1
  if (weight == "aligned_length" && any(unit[q] != unit[s]))
    ggchord_stop(caller, ": aligned length requires comparable endpoint units")
  values <- switch(weight, count = rep(1, nrow(links)), aligned_length = aligned, links[[weight]])
  if (!is.numeric(values) || any(!is.finite(values) | values < 0))
    ggchord_stop(caller, ": missing or invalid nonnegative link weights")
  direction <- ggchord_reference_direction(links, caller)
  adjacency <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  for (i in which(q != s)) {
    adjacency[q[i], s[i]] <- adjacency[q[i], s[i]] + values[i]
    adjacency[s[i], q[i]] <- adjacency[s[i], q[i]] + values[i]
  }
  if (is.null(order)) {
    order <- reference %||% ids[which.max(rowSums(adjacency))]
    while (length(order) < length(ids)) {
      remaining <- setdiff(ids, order)
      order <- c(order, remaining[which.max(rowSums(adjacency[remaining, order, drop = FALSE]))])
    }
  }
  reference <- reference %||% order[1L]
  orientation <- setNames(rep(1, length(ids)), ids)
  evidence <- data.frame(accver = order, forward_support = 0, reverse_support = 0,
    confidence = NA_real_, orientation = 1, status = "unresolved")
  evidence$status[order == reference] <- "reference"
  evidence$source_rows <- rep(list(integer()), length(ids))
  visited <- reference
  for (id in setdiff(order, reference)) {
    at <- match(id, order)
    query <- q == id & s %in% visited
    subject <- s == id & q %in% visited
    rows <- which((query | subject) & !is.na(direction) & values > 0)
    peer <- ifelse(query[rows], s[rows], q[rows])
    vote <- direction[rows] * orientation[peer]
    forward <- sum(values[rows][vote > 0])
    reverse <- sum(values[rows][vote < 0])
    total <- forward + reverse
    confidence <- if (total > 0) abs(forward - reverse) / total else NA_real_
    evidence$source_rows[[at]] <- rows
    evidence[at, c("forward_support", "reverse_support", "confidence")] <-
      list(forward, reverse, confidence)
    if (forward != reverse && confidence >= min_support) {
      orientation[id] <- if (forward > reverse) 1 else -1
      evidence$status[at] <- "supported"
    }
    visited <- c(visited, id)
  }
  evidence$orientation <- unname(orientation[order])
  unresolved <- evidence[evidence$status == "unresolved", c("accver", "confidence")]
  unresolved$reason <- ifelse(rowSums(adjacency)[unresolved$accver] == 0, "isolated",
    ifelse(is.na(unresolved$confidence), "missing_direction", "insufficient_support"))
  list(layout = list(order = order, orientation = orientation[order]),
    evidence = evidence, diagnostics = unresolved)
}
