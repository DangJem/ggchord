#' Radial spokes for multi-sequence synteny
#'
#' Displays linear sequences as equally spaced outward spokes. Adjacent
#' spokes, including the last and first, may be joined by existing interval
#' alignments. No missing relationship is created. Two spokes display each
#' input HSP once; one spoke is not a multi-sequence synteny view.
#'
#' Uses the same sequence, feature, track, label and link layers as
#' [coord_collinear()]. `geom_seq(seq_order = ...)` orders spokes and declares
#' named repeated display instances; ambiguous relationships require `between`.
#' Lengths are proportional by default; `geom_seq(length_mode = "equal")`
#' explicitly normalizes display lengths while retaining native positions.
#' Positive feature/track offsets point counterclockwise from each outward
#' spoke, independently of biological direction and display order. Sequence
#' names are placed just beyond the outer ends to keep the central hole clear.
#'
#' Ribbons default to cyclic neighbors. Use explicit `pairs` or
#' `adjacent_only = FALSE` for nonadjacent ribbons; point links to nonneighbors
#' require explicit `between`. Default curved connections stay in an annular
#' channel, retaining the complete interval fronts and inversion twists.
#' Explicit control points use final layout x/y units. Shared link branches,
#' circular-only annotations, shared straight scale bars and joined multi-window
#' rows are unsupported. Use the native position axis on each spoke.
#'
#' Interval fronts use closed one-based interval boundaries
#' `[min(start,end)-1, max(start,end)]` in input direction. A one-position
#' interval retains nonzero width. Axis ticks and position mapping use
#' one-based positions; this new view does not change existing view conventions.
#'
#' @param start_angle Direction of the first outward spoke in degrees,
#'   counterclockwise from the positive x axis. Default 90 starts at noon.
#' @param clockwise Logical, order subsequent spokes clockwise, default TRUE.
#' @param inner_radius Positive inner radius in layout units, default 0.6.
#'   The longest proportional spoke has length 2. The central hole is display
#'   spacing and does not remove any source positions.
#' @inheritParams coord_collinear
#' @return A ggchord coordinate object.
#' @export
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' ggchord(seq_data_example, ribbon_data_example) +
#'   geom_seq() +
#'   geom_link_ribbon() +
#'   geom_seq_label() +
#'   coord_radial_synteny(start_angle = 90, inner_radius = 0.7)
coord_radial_synteny <- function(start_angle = 90, clockwise = TRUE,
    inner_radius = 0.6, ratio = 1, xlim = NULL, ylim = NULL,
    expand = FALSE, clip = "off", fit = c("labels", "geometry", "manual")) {
  for (name in c("start_angle", "inner_radius")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value))
      ggchord_stop("coord_radial_synteny(): ", name, " must be a finite number")
  }
  if (inner_radius <= 0) ggchord_stop("coord_radial_synteny(): inner_radius must be positive")
  if (!is.logical(clockwise) || length(clockwise) != 1L || is.na(clockwise))
    ggchord_stop("coord_radial_synteny(): clockwise must be TRUE or FALSE")
  base <- coord_collinear(ratio = ratio, xlim = xlim, ylim = ylim,
    expand = expand, clip = clip, fit = fit)
  ggplot2::ggproto("CoordRadialSynteny", base, ggchord_radial = TRUE,
    radial_start_angle = start_angle, radial_clockwise = clockwise,
    radial_inner_radius = inner_radius)
}

# An isometry per source frame: along-axis distances and physical band widths
# retain their units. Projection never rotates/rewrites the source tables.
ggchord_radial_points <- function(x, y, ids, frames) {
  at <- match(ids, frames$key)
  r <- frames$inner[at] + x - frames$axis_min[at]
  v <- y - frames$source_y[at]
  cbind(x = r * frames$tx[at] + v * frames$nx[at],
    y = r * frames$ty[at] + v * frames$ny[at])
}

ggchord_project_radial_data <- function(data, frames, id = NULL, fixed_angle = FALSE) {
  if (!is.data.frame(data) || !nrow(data)) return(data)
  ids <- data$accver %||% rep(id, nrow(data))
  if (!length(ids)) return(data)
  for (pair in list(c("x", "y"), c("x0", "y0"), c("x1", "y1"),
      c("xend", "yend"), c("text_x", "text_y"), c("label_x", "label_y"), c("anchor_x", "anchor_y"))) {
    if (!all(pair %in% names(data))) next
    xy <- ggchord_radial_points(data[[pair[1]]], data[[pair[2]]], ids, frames)
    data[[pair[1]]] <- xy[, 1]; data[[pair[2]]] <- xy[, 2]
  }
  if (fixed_angle) return(data)
  at <- match(ids, frames$key)
  for (field in intersect(c("angle", "text_angle", "label_angle"), names(data))) {
    angle <- data[[field]] * pi / 180 + frames$theta[at]
    value <- atan2(sin(angle) * frames$ratio[at], cos(angle)) * 180 / pi
    flipped <- which(value > 90 | value < -90)
    above <- which(value > 90)
    value[above] <- value[above] - 180
    below <- which(value < -90)
    value[below] <- value[below] + 180
    data[[field]] <- value
    # Flipping the baseline requires reflecting both justifications as well.
    prefix <- if (field == "label_angle") "label_" else ""
    for (just in if (field != "angle" || !"text_angle" %in% names(data)) paste0(prefix, c("hjust", "vjust")) else character())
      if (just %in% names(data) && is.numeric(data[[just]]))
        data[[just]][flipped] <- 1 - data[[just]][flipped]
  }
  data
}

# Polar interpolation keeps the default connection out of the central hole.
# Endpoint r is interpolated smoothly, so the curve meets each front along
# its local normal, including the two-spoke semicircular channel.
ggchord_radial_edge <- function(a, b, delta, controls = NULL, n = 50L) {
  if (!is.null(controls)) return(bezier_pts(a, b, controls[[1]], controls[[2]], n = n))
  t <- seq(0, 1, length.out = n)
  ra <- sqrt(sum(a^2)); rb <- sqrt(sum(b^2))
  aa <- atan2(a[2], a[1]); ab <- atan2(b[2], b[1])
  change <- (ab - aa + pi) %% (2 * pi) - pi
  # Choose the declared channel, not the opposite semicircle in N=2.
  if (delta > 0 && change < 0) change <- change + 2 * pi
  if (delta < 0 && change > 0) change <- change - 2 * pi
  r <- ra + (rb - ra) * (3 * t^2 - 2 * t^3)
  angle <- aa + change * t
  cbind(x = r * cos(angle), y = r * sin(angle))
}

ggchord_radial_connections <- function(d, input, params, component, layout, frames) {
  if (!is.data.frame(d) || !nrow(d) || is.null(input)) return(d)
  if (!identical(params$link_branch %||% "none", "none"))
    ggchord_stop("coord_radial_synteny(): shared branches are unsupported; use link_branch='none'")
  ref <- layout$sequence_reference
  configured <- params$ribbon_gap %||% params$link_gap
  gaps <- process_sequence_param(configured %||% .035, layout$seqs, "link_gap", .035)
  pieces <- list(); endpoints <- list()
  for (i in unique(d$source_row)) {
    template <- d[d$source_row == i, , drop = FALSE]
    q <- as.character(input$qaccver[i]); s <- as.character(input$saccver[i])
    qi <- match(q, frames$key); si <- match(s, frames$key); count <- nrow(frames)
    forward <- (si - qi) %% count
    adjacent <- forward %in% c(1L, count - 1L)
    if (component == "link" && !adjacent && is.null(params$between))
      ggchord_stop("coord_radial_synteny(): nonadjacent point links require explicit between")
    # For two spokes, anchor the one physical channel to declared-order pair,
    # independent of input query/subject order. Reciprocal source rows survive.
    step <- if (forward <= count / 2) forward else forward - count
    if (count == 2L) step <- if (qi < si) 1 else -1
    delta <- step * (if (isTRUE(layout$radial_clockwise)) -1 else 1) * 2 * pi / count
    normal_q <- sign(delta); normal_s <- -sign(delta)
    front <- function(id, start, end, sign) {
      if (!is.numeric(start) || !is.numeric(end) || any(!is.finite(c(start, end))) ||
          any(c(start, end) < 1 | c(start, end) > ref$lens[id]))
        ggchord_stop("coord_radial_synteny(): endpoints must be finite and within sequence length")
      bounds <- if (component != "ribbon") c(start, end) else
        if (start <= end) c(start - 1, end) else c(start, end - 1)
      pos <- seq(bounds[1], bounds[2], length.out = if (component == "ribbon") 16L else 1L)
      gap <- ggchord_collinear_gap_profile(pos, id, sign, layout$obstacles,
        avoid = if (is.null(configured)) params$link_avoid %||% "none" else "none",
        close_gap = unname(gaps[id]))
      r <- frames$inner[match(id, frames$key)] + ggchord_collinear_x(pos, id, ref) - frames$axis_min[match(id, frames$key)]
      # Oversized fixed/automatic bands would invert the front or consume its
      # angular channel. Do not silently clamp the user's physical geometry.
      if (any(!is.finite(gap)) || any(gap < 0) || any(atan2(gap, r) >= abs(delta) * .45))
        ggchord_stop("coord_radial_synteny(): endpoint clearance does not fit its angular channel; increase inner_radius or reduce bands/gaps")
      ggchord_radial_points(ggchord_collinear_x(pos, id, ref), ref$y[id] + sign * gap, id, frames)
    }
    if (component == "ribbon") {
      qxy <- front(q, input$qstart[i], input$qend[i], normal_q)
      sxy <- front(s, input$sstart[i], input$send[i], normal_s)
      edge <- function(a, b, which) {
        cp <- if (is.null(params$ribbon_ctrl_point)) NULL else
          rep(ggchord_link_controls(params$ribbon_ctrl_point, i, (a + b) / 2)[which], 2L)
        ggchord_radial_edge(a, b, delta, cp)
      }
      a <- edge(qxy[1, ], sxy[1, ], 1L)
      b <- edge(qxy[nrow(qxy), ], sxy[nrow(sxy), ], 2L)
      xy <- rbind(qxy, b, sxy[nrow(sxy):1L, , drop = FALSE], a[nrow(a):1L, , drop = FALSE])
    } else {
      qxy <- front(q, input$qpos[i], input$qpos[i], normal_q)
      sxy <- front(s, input$spos[i], input$spos[i], normal_s)
      cp <- if (is.null(params$link_ctrl_point)) NULL else
        ggchord_link_controls(params$link_ctrl_point, i, (qxy[1, ] + sxy[1, ]) / 2)
      xy <- if (identical(params$link_type, "straight")) rbind(qxy, sxy) else
        ggchord_radial_edge(qxy[1, ], sxy[1, ], delta, cp)
    }
    piece <- template[rep(1L, nrow(xy)), , drop = FALSE]
    piece$x <- xy[, 1]; piece$y <- xy[, 2]; piece$group <- i
    piece$.ggchord_adjacent <- adjacent
    if (component == "ribbon") {
      piece$.q_n <- nrow(qxy); piece$.s_n <- nrow(sxy)
      endpoints[[length(endpoints) + 1L]] <- data.frame(source_row = i,
        qaccver = q, saccver = s, qstart = input$qstart[i], qend = input$qend[i],
        sstart = input$sstart[i], send = input$send[i],
        qboundary_start = if (input$qstart[i] <= input$qend[i]) input$qstart[i] - 1 else input$qstart[i],
        qboundary_end = if (input$qstart[i] <= input$qend[i]) input$qend[i] else input$qend[i] - 1,
        sboundary_start = if (input$sstart[i] <= input$send[i]) input$sstart[i] - 1 else input$sstart[i],
        sboundary_end = if (input$sstart[i] <= input$send[i]) input$send[i] else input$send[i] - 1,
        qx_start = qxy[1, 1], qy_start = qxy[1, 2], qx_end = tail(qxy[, 1], 1), qy_end = tail(qxy[, 2], 1),
        sx_start = sxy[1, 1], sy_start = sxy[1, 2], sx_end = tail(sxy[, 1], 1), sy_end = tail(sxy[, 2], 1),
        adjacent = adjacent, channel_angle = delta)
    }
    pieces[[length(pieces) + 1L]] <- piece
  }
  list(data = ggchord_rbind_fill(pieces), endpoints = ggchord_rbind_fill(endpoints))
}

ggchord_project_radial_layout <- function(layout, plot) {
  if (isTRUE(layout$radial_projected)) return(layout)
  n <- length(layout$seqs)
  if (n < 2L) ggchord_stop("coord_radial_synteny(): at least two displayed sequences are required")
  theta <- (plot$coordinates$radial_start_angle +
    (if (plot$coordinates$radial_clockwise) -1 else 1) * 360 * (seq_len(n) - 1) / n) * pi / 180
  tab <- layout$sequence_layout
  frames <- data.frame(key = layout$seqs, theta = theta, tx = cos(theta), ty = sin(theta),
    nx = -sin(theta), ny = cos(theta), inner = plot$coordinates$radial_inner_radius,
    source_y = unname(layout$sequence_reference$y[layout$seqs]), axis_min = tab$xmin,
    ratio = plot$coordinates$ratio %||% 1)
  layout$radial_clockwise <- plot$coordinates$radial_clockwise
  endpoints <- list()
  for (owner in names(layout$layer_geometry)) for (component in names(layout$layer_geometry[[owner]])) {
    d <- layout$layer_geometry[[owner]][[component]]
    if (!component %in% c("ribbon", "link")) {
      layer <- plot$layers[[which(vapply(plot$layers, function(l) identical(l$ggchord_layer_id, owner), logical(1)))[1]]]
      params <- layer$ggchord_params
      if (component == "seq_label" && nrow(d)) {
        at <- match(d$accver, frames$key)
        d$text_x <- tab$xmax[at] + .08
        d$text_y <- frames$source_y[at]
        horizontal <- identical(params$seq_label_orientation, "horizontal")
        upright <- abs(frames$tx[at]) < .2
        if (is.null(params$seq_label_hjust)) d$hjust <- if (horizontal)
          ifelse(upright, .5, ifelse(frames$tx[at] > 0, 0, 1)) else 0
        if (is.null(params$seq_label_vjust)) d$vjust <- if (horizontal)
          ifelse(upright, ifelse(frames$ty[at] > 0, 0, 1), .5) else .5
      }
      fixed <- identical(params$seq_label_orientation, "horizontal") ||
        identical(params$gene_label_orientation, "horizontal") || !is.null(params$gene_label_rotation)
      layout$layer_geometry[[owner]][[component]] <- ggchord_project_radial_data(d, frames, fixed_angle = fixed)
      next
    }
    if (!is.data.frame(d) || !nrow(d)) next
    input <- layout$layer_inputs[[owner]][[component]]
    layer <- plot$layers[[which(vapply(plot$layers, function(l) identical(l$ggchord_layer_id, owner), logical(1)))[1]]]
    result <- ggchord_radial_connections(d, input, layer$ggchord_params, component, layout, frames)
    layout$layer_geometry[[owner]][[component]] <- result$data
    if (!is.null(result$endpoints)) {
      result$endpoints$layer_id <- owner
      result$endpoints$source_hsp_id <- if ("source_hsp_id" %in% names(input))
        as.character(input$source_hsp_id[result$endpoints$source_row]) else NA_character_
      endpoints[[owner]] <- result$endpoints
    }
  }
  layout$seq_arcs <- lapply(layout$seq_arcs, ggchord_project_radial_data, frames = frames)
  axis_angle <- ggchord_plot_settings(plot)$axis$text.orientation
  fixed_axis <- is.numeric(axis_angle) || identical(axis_angle, "horizontal")
  for (field in c("axis_lines", "axis_ticks", "gene_polys", "gene_labels", "gene_label_segments",
      "region_polys", "seq_labels_df", "sequence_scales", "track_axes", "seq_breaks"))
    layout[[field]] <- ggchord_project_radial_data(layout[[field]], frames,
      fixed_angle = field == "axis_ticks" && fixed_axis)
  for (pair in list(c("ribbon_polys", "ribbon"), c("link_lines", "link"))) {
    combined <- ggchord_rbind_fill(lapply(layout$layer_geometry, function(owner) owner[[pair[2]]]))
    layout[[pair[1]]] <- if (is.data.frame(combined) && nrow(combined)) combined else NULL
  }
  layout$seq_labels_df <- ggchord_rbind_fill(lapply(layout$layer_geometry, function(owner) owner$seq_label)) %||% data.frame()
  # Expose both final bounds and the unchanged local axis used for inversion.
  tab$axis_min <- tab$xmin; tab$axis_max <- tab$xmax
  a <- ggchord_radial_points(tab$axis_min, frames$source_y, frames$key, frames)
  b <- ggchord_radial_points(tab$axis_max, frames$source_y, frames$key, frames)
  tab$xmin <- pmin(a[,1], b[,1]); tab$xmax <- pmax(a[,1], b[,1])
  tab$ymin <- pmin(a[,2], b[,2]); tab$ymax <- pmax(a[,2], b[,2])
  tab$y <- (tab$ymin + tab$ymax) / 2
  tab$theta <- theta; tab$inner_radius <- frames$inner
  layout$sequence_layout <- tab
  layout$radial_frames <- frames
  layout$radial_endpoints <- ggchord_rbind_fill(endpoints) %||% data.frame()
  geometry <- c(unlist(layout$layer_geometry, recursive = FALSE), list(layout$axis_lines, layout$axis_ticks))
  xx <- yy <- numeric()
  for (d in geometry) if (is.data.frame(d)) for (pair in list(c("x", "y"), c("x0", "y0"), c("x1", "y1"), c("text_x", "text_y"), c("label_x", "label_y")))
    if (all(pair %in% names(d))) { xx <- c(xx, d[[pair[1]]]); yy <- c(yy, d[[pair[2]]]) }
  layout$extremes <- list(x_min = min(xx, na.rm = TRUE), x_max = max(xx, na.rm = TRUE),
    y_min = min(yy, na.rm = TRUE), y_max = max(yy, na.rm = TRUE))
  layout$seq_axis <- "radial"; layout$radial_projected <- TRUE
  layout
}

ggchord_map_radial <- function(plot, accver, position, inverse, source) {
  if (!is.logical(inverse) || length(inverse) != 1L || is.na(inverse) ||
      !is.logical(source) || length(source) != 1L || is.na(source))
    ggchord_stop("map_ggchord_position(): inverse/source must be TRUE or FALSE")
  if (inverse) {
    if (!(is.matrix(position) || is.data.frame(position)) || ncol(position) != 2L ||
        !all(c("x", "y") %in% colnames(position)))
      ggchord_stop("Radial inverse mapping requires a matrix/data frame with x and y columns")
    xy <- as.matrix(position[, c("x", "y"), drop = FALSE])
    if (!is.numeric(xy) || any(!is.finite(xy))) ggchord_stop("Radial inverse mapping requires finite x/y")
    count <- nrow(xy)
  } else {
    if (!is.numeric(position) || any(!is.finite(position))) ggchord_stop("Supply finite sequence positions")
    count <- length(position)
  }
  if (anyNA(accver) || !length(accver) %in% c(1L, count)) ggchord_stop("Supply matching accessions")
  layout <- get_chord_layout(plot)
  frames <- layout$radial_frames; ref <- layout$sequence_reference
  ids <- rep_len(as.character(accver), count)
  if (!is.null(layout$display_instances)) {
    instances <- layout$display_instances
    ids <- vapply(ids, function(id) {
      if (id %in% instances$instance_id) return(instances$key[match(id, instances$instance_id)])
      at <- which(instances$accver == id)
      if (length(at) != 1L) ggchord_stop("Ambiguous display instance; specify a named instance ID")
      instances$key[at]
    }, character(1))
  }
  at <- match(ids, frames$key)
  if (anyNA(at)) ggchord_stop("Unknown displayed Radial accession")
  origin <- layout$sequence_layout$source_origin[match(ids, layout$sequence_layout$accver)]
  if (inverse) {
    r <- xy[,1] * frames$tx[at] + xy[,2] * frames$ty[at]
    normal <- xy[,1] * frames$nx[at] + xy[,2] * frames$ny[at]
    if (any(abs(normal) > 1e-7)) ggchord_stop("Radial inverse position must lie on the specified backbone")
    local <- (r - frames$inner[at] + frames$axis_min[at] - unname(ref$x_offset[ids])) / ggchord_collinear_scale(ids, ref)
    pos <- ifelse(ref$orientation[ids] == 1, local, ref$lens[ids] - local)
  } else {
    pos <- position - if (source) origin - 1 else 0
    xy <- ggchord_radial_points(ggchord_collinear_x(pos, ids, ref), ref$y[ids], ids, frames)
  }
  if (any(pos < 1 - 1e-7 | pos > ref$lens[ids] + 1e-7))
    ggchord_stop("Radial position lies outside the specified sequence/window")
  out <- data.frame(accver = ids, position = unname(pos + if (source) origin - 1 else 0),
    x = xy[,1], y = xy[,2], sequence_unit = layout$sequence_unit, stringsAsFactors = FALSE)
  ggchord_restore_instance_data(out, layout$display_instances)
}
