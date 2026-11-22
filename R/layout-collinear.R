# Linear sequence-local geometry. Positions and normal offsets are kept
# separate: no circular radius or angular reference is created for this coord.

ggchord_collinear_x <- function(pos, id, ref) {
  sequence_length <- unname(ref$lens[id])
  orientation <- unname(ref$orientation[id])
  ifelse(rep_len(orientation, length(pos)) == 1,
    pos, rep_len(sequence_length, length(pos)) - pos) * ggchord_collinear_scale(id, ref) +
    rep_len(unname(ref$x_offset[id] %||% 0), length(pos))
}

ggchord_collinear_scale <- function(id, ref) {
  if (length(ref$x_scale) == 1L) unname(ref$x_scale) else unname(ref$x_scale[id])
}

# Obstacles record signed distances from a horizontal backbone. Only entities
# on the side facing the destination can move a connection endpoint.
ggchord_collinear_gap_profile <- function(
    pos, id, direction, obstacles, avoid = "none", close_gap = .035) {
  out <- rep(close_gap, length(pos))
  if (identical(avoid, "none") || is.null(obstacles) || !nrow(obstacles) ||
      !length(pos)) return(out)
  obs <- obstacles[obstacles$accver == id, , drop = FALSE]
  if (!nrow(obs)) return(out)
  clearance <- if (direction > 0) obs$normal_max else -obs$normal_min
  keep <- is.finite(clearance) & clearance > 0
  obs <- obs[keep, , drop = FALSE]
  clearance <- clearance[keep]
  if (!nrow(obs)) return(out)
  if (identical(avoid, "uniform")) {
    hit <- obs$start <= max(pos) & obs$end >= min(pos)
    if (any(hit)) out[] <- max(close_gap, clearance[hit] + .025)
    return(out)
  }
  transition <- max(diff(range(pos)) * .08, 1)
  for (i in seq_len(nrow(obs))) {
    distance <- pmax(obs$start[i] - pos, pos - obs$end[i], 0)
    weight <- ifelse(distance < transition,
      (1 + cos(pi * pmin(distance / transition, 1))) / 2, 0)
    out <- pmax(out,
      close_gap + pmax(0, clearance[i] + .025 - close_gap) * weight)
  }
  out
}

ggchord_collinear_controls <- function(from, to, control = NULL,
                                       row = 1L, edge = NULL) {
  if (!is.null(control)) {
    controls <- ggchord_link_controls(control, row, (from + to) / 2)
    if (!is.null(edge)) controls <- rep(controls[edge], 2L)
    return(controls)
  }
  middle_y <- (from[2L] + to[2L]) / 2
  list(c(from[1L], middle_y), c(to[1L], middle_y))
}

# The normal derivative of a cubic whose controls swap the two endpoints
# vanishes at the midpoint. Keep enough separation there to preserve the
# apparent width of a ribbon, even at maximum requested curvature.
ggchord_connection_bend <- function(route) {
  strength <- route$curvature %||% .5
  strength - max(0, strength - .5) / 2
}

# A positive bulge requests an outside route for non-neighbouring rows.
# Solve the cubic at each intervening backbone, using the inside ribbon edge.
# This also works for flipped references and separately aligned display windows.
ggchord_outside_controls <- function(route, frames, endpoints, from, to) {
  if (!inherits(route, "ggchord_link_route") || !isTRUE(route$explicit_bulge) || route$bulge <= 0 ||
      !identical(route$route, "curve") || from[2L] == to[2L]) return(NULL)
  fraction <- (frames$y - from[2L]) / (to[2L] - from[2L])
  between <- which(fraction > 0 & fraction < 1 & !frames$accver %in% endpoints)
  if (!length(between)) return(NULL)
  bend <- ggchord_connection_bend(route)
  y_curve <- function(t) 3 * (1 - t)^2 * t * bend +
    3 * (1 - t) * t^2 * (1 - bend) + t^3
  t <- vapply(fraction[between], function(f)
    stats::uniroot(function(t) y_curve(t) - f, c(0, 1), tol = 1e-8)$root, numeric(1))
  clearance <- route$bulge * max(frames$xmax - frames$xmin)
  control_x <- max((frames$xmax[between] + clearance -
    (1 - t)^3 * from[1L] - t^3 * to[1L]) / (3 * t * (1 - t)))
  dy <- to[2L] - from[2L]
  list(c(control_x, from[2L] + bend * dy), c(control_x, to[2L] - bend * dy))
}

ggchord_collinear_ribbons <- function(p, layout) {
  input <- p$ribbon_data
  if (is.null(input) || !nrow(input)) return(NULL)
  ref <- layout$sequence_reference
  valid <- which(input$qaccver %in% layout$seqs & input$saccver %in% layout$seqs)
  if (!length(valid)) {
    warning("No valid alignment data available for plotting", call. = FALSE)
    return(NULL)
  }
  scheme <- p$ribbon_color_scheme %||% "pident"
  alpha <- rep(p$ribbon_alpha %||% .42, nrow(input))
  if (!is.null(p$ribbon_alpha_by)) {
    av <- as.numeric(input[[p$ribbon_alpha_by]][valid])
    alpha_range <- p$ribbon_alpha_range %||% c(.15, .9)
    span <- range(av, na.rm = TRUE)
    normalized <- if (length(av) > 1L && diff(span) > 0) {
      (av - span[1L]) / diff(span)
    } else rep(.5, length(av))
    alpha[valid] <- alpha[valid] *
      (alpha_range[1L] + diff(alpha_range) * normalized)
  }
  discrete_values <- function(values, supplied, defaults) {
    levels <- unique(values)
    if (is.null(supplied)) return(setNames(defaults(length(levels)), levels))
    if (is.null(names(supplied))) {
      if (!length(supplied) %in% c(1L, length(levels))) {
        ggchord_stop("Ribbon style values must have length 1 or match their categories")
      }
      return(setNames(rep_len(supplied, length(levels)), levels))
    }
    unknown <- setdiff(names(supplied), levels)
    if (length(unknown)) ggchord_stop("Ribbon style values contain unknown categories")
    out <- setNames(defaults(length(levels)), levels)
    out[names(supplied)] <- supplied
    out
  }
  outline_map <- if (!is.null(p$ribbon_outline_by)) discrete_values(
    as.character(input[[p$ribbon_outline_by]][valid]),
    p$ribbon_outline_colors, chord_default_palette
  ) else NULL
  linetype_map <- if (!is.null(p$ribbon_linetype_by)) discrete_values(
    as.character(input[[p$ribbon_linetype_by]][valid]), p$ribbon_linetypes,
    function(n) rep_len(c("solid", "dashed", "dotted", "dotdash",
      "longdash", "twodash"), n)
  ) else NULL
  obstacle_index <- ggchord_obstacle_index(p$ribbon_obstacles, p$link_avoid)
  pieces <- lapply(seq_along(valid), function(j) {
    i <- valid[j]
    q <- as.character(input$qaccver[i]); s <- as.character(input$saccver[i])
    direction_q <- sign(ref$y[s] - ref$y[q])
    channel <- abs(unname(ref$y[q] - ref$y[s]))
    same_row <- channel == 0
    direction_s <- if (same_row) 1 else -direction_q
    if (same_row) direction_q <- 1
    make_front <- function(id, start, end, direction) {
      bounds <- ggchord_interval_front(start, end)
      normal_mm <- p$ribbon_normal_mm
      front <- if (isTRUE(p$ribbon_gap_auto)) {
        gap_mm <- p$ribbon_gap_mm[[id]]
        base <- ggchord_ribbon_base_gap(p$ribbon_clearance, id, direction, normal_mm, gap_mm)
        ggchord_connection_front(bounds, obstacle_index[[id]], direction, p$link_avoid,
          base, normal_mm, ggchord_collinear_scale(id, ref) * normal_mm /
            ggchord_units_per_mm(if (isTRUE(p$collinear_vertical)) "y" else "x"),
          gap_mm = gap_mm, strength = p$ribbon_avoid_strength)
      } else list(position = bounds, gap = rep(unname(p$ribbonGap[[id]]), 2L))
      list(xy = cbind(x = ggchord_collinear_x(front$position, id, ref),
        y = unname(ref$y[id]) + direction * front$gap + unname(p$ribbon_offset[id])),
        gap = max(front$gap))
    }
    q_front <- make_front(q, input$qstart[i], input$qend[i], direction_q)
    s_front <- make_front(s, input$sstart[i], input$send[i], direction_s)
    q_xy <- q_front$xy; s_xy <- s_front$xy
    ggchord_collinear_check_gap(q, s, q_front$gap, s_front$gap, ref,
      "geom_link_ribbon()")
    outside <- ggchord_outside_controls(p$ribbon_ctrl_point, layout$sequence_layout, c(q, s),
      c(min(q_xy[, 1]), q_xy[1, 2]), c(min(s_xy[, 1]), s_xy[1, 2]))
    edge <- function(index, control_edge) {
      from <- q_xy[index[1L], ]; to <- s_xy[index[2L], ]
      route <- p$ribbon_ctrl_point
      if (same_row) {
        if (identical(route$route, "straight"))
          ggchord_stop("geom_link_ribbon(): a same-row ribbon requires a curved route")
        rise <- abs(to[1L] - from[1L]) * (route$curvature %||% .5)
        return(bezier_pts(from, to, from + c(0, rise), to + c(0, rise), n = 50L))
      }
      if (inherits(route, "ggchord_link_route")) {
        if (route$route == "straight") return(rbind(from, to))
        bend <- ggchord_connection_bend(route) * (to[2L] - from[2L])
        cp <- outside %||% list(c(from[1L], from[2L] + bend),
          c(to[1L], to[2L] - bend))
      } else cp <- ggchord_collinear_controls(from, to, route, i, edge = control_edge)
      bezier_pts(from, to, cp[[1L]], cp[[2L]], n = 50L)
    }
    first_edge <- edge(c(1L, 1L), 1L)
    last_edge <- edge(c(nrow(q_xy), nrow(s_xy)), 2L)
    xy <- rbind(q_xy, last_edge, s_xy[nrow(s_xy):1L, , drop = FALSE],
      first_edge[nrow(first_edge):1L, , drop = FALSE])
    biological_direction <- input$.ggchord_relative_direction[i]
    opacity <- alpha[i]
    if (identical(p$ribbon_direction, "alpha") && biological_direction != "unknown") {
      opacity <- opacity * p$ribbon_direction_alpha[biological_direction]
    }
    outline <- if (!is.null(outline_map)) {
      outline_map[as.character(input[[p$ribbon_outline_by]][i])]
    } else if (identical(p$ribbon_direction, "outline")) {
      p$ribbon_direction_colors[biological_direction]
    } else "black"
    linetype <- if (!is.null(linetype_map)) {
      linetype_map[as.character(input[[p$ribbon_linetype_by]][i])]
    } else if (identical(p$ribbon_direction, "linetype")) {
      p$ribbon_direction_linetypes[biological_direction]
    } else if (same_row || !isTRUE(p$collinear_crossrow_styles) ||
        isTRUE(input$.ggchord_adjacent[i])) "solid" else "dashed"
    out <- data.frame(x = xy[, 1L], y = xy[, 2L], group = j,
      source_row = input$.ggchord_source_row[i] %||% i,
      alpha = unname(opacity), outline_col = unname(outline),
      linetype_val = unname(linetype), direction = biological_direction,
      .ggchord_adjacent = isTRUE(input$.ggchord_adjacent[i]),
      .q_n = nrow(q_xy), .s_n = nrow(s_xy),
      q_gap = q_front$gap, s_gap = s_front$gap, stringsAsFactors = FALSE)
    if (identical(scheme, "pident")) {
      out$pident <- input$pident[i] %||% NA_real_
    } else if (identical(scheme, "value")) {
      out$value <- as.numeric(input[[p$ribbon_color_by]][i])
    } else {
      out$fill <- unname(switch(scheme,
        single = p$ribbon_colors[[1L]],
        query = p$ribbon_colors[q], subject = p$ribbon_colors[s]))
    }
    out
  })
  do.call(rbind, pieces)
}

ggchord_collinear_link_geometry <- function(data, params, layout, rows) {
  ref <- layout$sequence_reference
  gaps <- process_sequence_param(params$link_gap %||% .035,
    layout$seqs, "link_gap", .035)
  offsets <- ggchord_connection_offsets(params, layout$seqs,
    if (identical(layout$seq_axis, "vertical")) "x" else "y")
  if (!is.numeric(gaps) || any(!is.finite(gaps))) {
    ggchord_stop("link_gap must contain finite numbers")
  }
  pieces <- lapply(rows, function(i) {
    q_id <- as.character(data$qaccver[i]); s_id <- as.character(data$saccver[i])
    endpoint <- function(id, pos, destination) {
      direction <- sign(ref$y[destination] - ref$y[id])
      gap <- unname(gaps[[id]])
      if (is.null(params$link_gap)) gap <- min(gap,
        abs(unname(ref$y[destination] - ref$y[id])) * .2)
      if (is.null(params$link_gap)) gap <- ggchord_collinear_gap_profile(
        pos, id, direction, layout$obstacles, params$link_avoid %||% "none", gap)
      c(ggchord_collinear_x(pos, id, ref), unname(ref$y[id]) + direction * gap + offsets[[id]])
    }
    q <- endpoint(q_id, data$qpos[i], s_id)
    s <- endpoint(s_id, data$spos[i], q_id)
    curvature <- ggchord_connection_bend(params$route)
    outside <- ggchord_outside_controls(params$route, layout$sequence_layout, c(q_id, s_id), q, s)
    controls <- if (ref$y[q_id] == ref$y[s_id]) {
      rise <- abs(q[1L] - s[1L]) * curvature
      list(q + c(0, rise), s + c(0, rise))
    } else outside %||% list(c(q[1L], q[2L] + (s[2L] - q[2L]) * curvature),
      c(s[1L], s[2L] + (q[2L] - s[2L]) * curvature))
    xy <- if (identical(params$link_type, "straight")) rbind(q, s) else
      bezier_pts(q, s, controls[[1L]], controls[[2L]], n = 60L)
    data.frame(x = xy[, 1L], y = xy[, 2L], group = i, source_row = i,
      .component = "branch", stringsAsFactors = FALSE)
  })
  if (length(pieces)) do.call(rbind, pieces) else
    data.frame(x = numeric(), y = numeric(), group = integer(),
      source_row = integer(), .component = character())
}

ggchord_collinear_entity_obstacles <- function(layer, layout, input) {
  polygons <- layout$gene_polys
  if (is.null(input) || !nrow(input) || !nrow(polygons)) {
    return(ggchord_empty_obstacles())
  }
  ref <- layout$sequence_reference
  stroke_mm <- ggchord_outline_half_mm(layer)
  pieces <- lapply(split(polygons, polygons$group), function(d) {
    row <- input[d$source_row[1L], , drop = FALSE]
    id <- as.character(row$accver)
    local_x <- (d$x - unname(ref$x_offset[id] %||% 0)) / ggchord_collinear_scale(id, ref)
    positions <- if (ref$orientation[id] == 1) local_x else
      unname(ref$lens[id]) - local_x
    offsets <- d$y - unname(ref$y[id])
    low <- min(offsets); high <- max(offsets)
    data.frame(accver = id, start = max(0, min(positions)),
      end = min(unname(ref$lens[id]), max(positions)),
      side = if (low >= 0) "positive" else if (high <= 0) "negative" else "both",
      normal_min = low, normal_max = high, priority = 1,
      source_layer = layer$ggchord_layer_id, stroke_mm = stroke_mm,
      profile = I(list(ggchord_obstacle_profile(positions, offsets))),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}

compute_collinear_layout <- function(...) {
  p <- list(...)
  seqs <- p$seqs; lens <- p$lens; n <- length(seqs)
  y <- setNames(if (n > 1L) c(0, -2 * cumsum(p$seq_gap[-n])) else 0,
    seqs)
  segment_table <- p$collinear_transform$table
  if ("row_index" %in% names(segment_table)) {
    first <- which(!duplicated(segment_table$row_index))
    rows <- length(first)
    gap <- p$seq_gap[first]
    for (i in unique(segment_table$row_index)) {
      values <- p$seq_gap[segment_table$row_index == i]
      if (length(unique(values)) > 1L) ggchord_stop(
        "coord_collinear(): windows in one display row require the same layout_sequences(gap)")
    }
    row_y <- if (rows > 1L) c(0, -2 * cumsum(gap[-rows])) else 0
    y <- setNames(row_y[segment_table$row_index], seqs)
  }
  if (isTRUE(p$collinear_right_normal)) y <- -y
  ref <- list(mode = "collinear", lens = lens, orientation = p$orientation,
    x_scale = if (identical(p$length_mode, "equal")) setNames(
      p$collinear_transform$table$units_per_bp, seqs) else p$collinear_transform$table$units_per_bp[1L], y = y,
    x_offset = p$collinear_transform$offset %||% setNames(rep(0, n), seqs))
  seq_arcs <- setNames(lapply(seqs, function(id) {
    data.frame(x = ggchord_collinear_x(c(0, lens[[id]]), id, ref),
      y = unname(y[id]), accver = id, stringsAsFactors = FALSE)
  }), seqs)
  empty <- data.frame()
  layout <- list(
    sequence_reference = ref,
    sequence_layout = cbind(p$collinear_transform$table,
      order = seq_len(n), y = unname(y)),
    orientation_diagnostics = p$collinear_transform$diagnostics,
    seq_arcs = seq_arcs, ribbon_polys = NULL, link_lines = NULL,
    region_polys = empty, gene_polys = empty, gene_labels = empty,
    gene_label_segments = empty, gene_label_clip_units = NA_real_,
    seq_labels_df = empty, axis_lines = empty, axis_ticks = empty,
    text_units_per_inch = ggchord_device_units_per_inch(
      c(0, 2), range(y)),
    seq_colors = p$seq_colors, seq_labels = p$seq_labels, seqs = seqs,
    ribbon_color_scheme = p$ribbon_color_scheme,
    ribbon_colors = p$ribbon_colors, ribbon_alpha = p$ribbon_alpha,
    ribbon_color_by = p$ribbon_color_by,
    ribbon_color_limits = p$ribbon_color_limits,
    ribbon_color_breaks = p$ribbon_color_breaks,
    ribbon_color_name = p$ribbon_color_name,
    ribbon_use_outline = !is.null(p$ribbon_outline_by) ||
      identical(p$ribbon_direction, "outline"),
    ribbon_use_linetype = !is.null(p$ribbon_linetype_by) ||
      identical(p$ribbon_direction, "linetype"),
    feature_shape_pal = p$feature_shape_pal,
    feature_shape_order = p$feature_shape_order,
    gene_pal = NULL, gene_color_scheme = p$gene_color_scheme,
    final_gene_order = NULL, gene_label_show = p$gene_label_show,
    gene_label_size = p$gene_label_size, gene_label_layout = p$gene_label_layout,
    show_axis = p$show_axis, axisLabelOrientation = p$axisLabelOrientation,
    rotation = 0, n_sequences = n, circular = FALSE, collinear = TRUE
  )
  layout$ribbon_polys <- ggchord_collinear_ribbons(p, layout)
  layout <- ggchord_collinear_annotations(p, layout)
  layout$extremes <- get_plot_extremes(
    allRibbon = layout$ribbon_polys, seqArcs = layout$seq_arcs,
    axisLines = layout$axis_lines, axisTicks = layout$axis_ticks,
    gene_polys = layout$gene_polys, gene_arrows = layout$gene_labels,
    seq_labels = layout$seq_labels_df, show_axis = layout$show_axis)
  class(layout) <- "chord_layout"
  layout
}
