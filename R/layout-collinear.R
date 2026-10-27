# Linear sequence-local geometry. Positions and normal offsets are kept
# separate: no circular radius or angular reference is created for this coord.

ggchord_collinear_x <- function(pos, id, ref) {
  sequence_length <- unname(ref$lens[id])
  orientation <- unname(ref$orientation[id])
  ifelse(rep_len(orientation, length(pos)) == 1,
    pos, rep_len(sequence_length, length(pos)) - pos) * ref$x_scale
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

ggchord_collinear_ribbons <- function(p, layout) {
  input <- p$ribbon_data
  if (is.null(input) || !nrow(input)) return(NULL)
  ref <- layout$sequence_reference
  valid <- which(input$qaccver != input$saccver &
    input$qaccver %in% layout$seqs & input$saccver %in% layout$seqs)
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
    if (is.null(supplied)) return(stats::setNames(defaults(length(levels)), levels))
    if (is.null(names(supplied))) {
      if (!length(supplied) %in% c(1L, length(levels))) {
        ggchord_stop("Ribbon style values must have length 1 or match their categories")
      }
      return(stats::setNames(rep_len(supplied, length(levels)), levels))
    }
    unknown <- setdiff(names(supplied), levels)
    if (length(unknown)) ggchord_stop("Ribbon style values contain unknown categories")
    out <- stats::setNames(defaults(length(levels)), levels)
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
  pieces <- lapply(seq_along(valid), function(j) {
    i <- valid[j]
    q <- as.character(input$qaccver[i]); s <- as.character(input$saccver[i])
    direction_q <- sign(ref$y[s] - ref$y[q])
    direction_s <- -direction_q
    make_front <- function(id, start, end, direction) {
      positions <- seq(start, end, length.out = 50L)
      configured <- unname(p$ribbonGap[[id]])
      automatic <- isTRUE(p$ribbon_gap_auto)
      avoid <- p$link_avoid %||% "none"
      obstacles <- p$ribbon_obstacles
      if (automatic && identical(avoid, "smooth") &&
          !is.null(obstacles) && nrow(obstacles)) {
        obs <- obstacles[obstacles$accver == id, , drop = FALSE]
        margin <- max(abs(end - start) * .08, 1)
        knots <- c(obs$start - margin, obs$start, obs$end, obs$end + margin)
        positions <- sort(unique(c(positions,
          knots[knots > min(start, end) & knots < max(start, end)])),
          decreasing = end < start)
      }
      gap <- if (automatic) ggchord_collinear_gap_profile(
        positions, id, direction, obstacles, avoid,
        close_gap = min(configured, .035)
      ) else rep(configured, length(positions))
      list(xy = cbind(x = ggchord_collinear_x(positions, id, ref),
        y = unname(ref$y[id]) + direction * gap), gap = max(gap))
    }
    q_front <- make_front(q, input$qstart[i], input$qend[i], direction_q)
    s_front <- make_front(s, input$sstart[i], input$send[i], direction_s)
    q_xy <- q_front$xy; s_xy <- s_front$xy
    edge <- function(index, control_edge) {
      from <- q_xy[index[1L], ]; to <- s_xy[index[2L], ]
      cp <- ggchord_collinear_controls(from, to, p$ribbon_ctrl_point,
        i, edge = control_edge)
      bezier_pts(from, to, cp[[1L]], cp[[2L]], n = 50L)
    }
    first_edge <- edge(c(1L, 1L), 1L)
    last_edge <- edge(c(nrow(q_xy), nrow(s_xy)), 2L)
    xy <- rbind(q_xy, last_edge, s_xy[nrow(s_xy):1L, , drop = FALSE],
      first_edge[nrow(first_edge):1L, , drop = FALSE])
    biological_direction <- if ("direction" %in% names(input) &&
        !is.na(input$direction[i]) && input$direction[i] %in% c("same", "reverse")) {
      as.character(input$direction[i])
    } else if ((input$qend[i] - input$qstart[i]) *
        (input$send[i] - input$sstart[i]) >= 0) "same" else "reverse"
    opacity <- alpha[i]
    if (identical(p$ribbon_direction, "alpha")) {
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
    } else "solid"
    out <- data.frame(x = xy[, 1L], y = xy[, 2L], group = j,
      source_row = i, alpha = unname(opacity), outline_col = unname(outline),
      linetype_val = unname(linetype), direction = biological_direction,
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

ggchord_collinear_link_geometry <- function(data, params, layout) {
  if (is.null(data)) ggchord_stop(
    "geom_link_line(): supply data with qaccver/saccver/qpos/spos")
  ggchord_require_columns(data, c("qaccver", "saccver", "qpos", "spos"),
    "geom_link_line()")
  ref <- layout$sequence_reference
  for (pair in list(c("qaccver", "qpos"), c("saccver", "spos"))) {
    ids <- as.character(data[[pair[1L]]]); pos <- data[[pair[2L]]]
    if (anyNA(ids) || any(!ids %in% names(ref$lens))) {
      ggchord_stop("geom_link_line(): unknown sequence accession")
    }
    if (!is.numeric(pos) || any(!is.finite(pos)) ||
        any(pos < 1 | pos > ref$lens[ids])) {
      ggchord_stop("geom_link_line(): positions must be finite and within sequence length")
    }
  }
  gaps <- process_sequence_param(params$link_gap %||% .035,
    layout$seqs, "link_gap", .035)
  if (!is.numeric(gaps) || any(!is.finite(gaps))) {
    ggchord_stop("link_gap must contain finite numbers")
  }
  pieces <- lapply(seq_len(nrow(data)), function(i) {
    q_id <- as.character(data$qaccver[i]); s_id <- as.character(data$saccver[i])
    endpoint <- function(id, pos, destination) {
      direction <- sign(ref$y[destination] - ref$y[id])
      gap <- unname(gaps[[id]])
      if (is.null(params$link_gap)) gap <- ggchord_collinear_gap_profile(
        pos, id, direction, layout$obstacles, params$link_avoid %||% "none", gap)
      c(ggchord_collinear_x(pos, id, ref), unname(ref$y[id]) + direction * gap)
    }
    q <- endpoint(q_id, data$qpos[i], s_id)
    s <- endpoint(s_id, data$spos[i], q_id)
    controls <- ggchord_collinear_controls(q, s, params$link_ctrl_point, i)
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
  pieces <- lapply(split(polygons, polygons$group), function(d) {
    row <- input[d$source_row[1L], , drop = FALSE]
    id <- as.character(row$accver)
    positions <- if (ref$orientation[id] == 1) d$x / ref$x_scale else
      unname(ref$lens[id]) - d$x / ref$x_scale
    offsets <- d$y - unname(ref$y[id])
    low <- min(offsets); high <- max(offsets)
    data.frame(accver = id, start = max(1, min(positions)),
      end = min(unname(ref$lens[id]), max(positions)),
      side = if (low >= 0) "positive" else if (high <= 0) "negative" else "both",
      normal_min = low, normal_max = high, priority = 1,
      source_layer = layer$ggchord_layer_id, stringsAsFactors = FALSE)
  })
  do.call(rbind, pieces)
}

compute_collinear_layout <- function(...) {
  p <- list(...)
  seqs <- p$seqs; lens <- p$lens; n <- length(seqs)
  y <- stats::setNames(if (n > 1L) c(0, -2 * cumsum(p$seq_gap[-n])) else 0,
    seqs)
  ref <- list(mode = "collinear", lens = lens, orientation = p$orientation,
    x_scale = 2 / max(lens), y = y)
  seq_arcs <- stats::setNames(lapply(seqs, function(id) {
    data.frame(x = ggchord_collinear_x(c(0, lens[[id]]), id, ref),
      y = unname(y[id]), accver = id, stringsAsFactors = FALSE)
  }), seqs)
  empty <- data.frame()
  layout <- list(
    sequence_reference = ref,
    sequence_layout = data.frame(accver = seqs, order = seq_len(n), xmin = 0,
      xmax = unname(lens[seqs]) * ref$x_scale, y = unname(y),
      orientation = unname(p$orientation[seqs]), units_per_bp = ref$x_scale,
      stringsAsFactors = FALSE),
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
