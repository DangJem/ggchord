# Sequence-local annotations for the Cartesian collinear layout. Positive
# normal offsets always mean up, including on reverse-oriented sequences.

ggchord_collinear_feature_parts <- function(shape, x0, x1, center, width,
                                            backbone, p, draw_head = TRUE,
                                            draw_start_head = FALSE,
                                            bidirectional = FALSE) {
  if (identical(shape, "lollipop")) {
    mid <- (x0 + x1) / 2
    stem_half <- max(.001, min(abs(x1 - x0) * .08, width * .1))
    theta <- seq(0, 2 * pi, length.out = 48L)
    return(list(
      cbind(x = mid + c(-1, 1, 1, -1) * stem_half,
        y = c(backbone, backbone, center, center)),
      cbind(x = mid + width * .58 * cos(theta),
        y = center + width * .58 * sin(theta))
    ))
  }
  # These shape templates use two local coordinates. A unit reference makes
  # their along-interval dimensions linear; no circular reference, radius or
  # physical track is stored in the resulting layout or exported geometry.
  short <- p$short_feature %||% "auto"
  if (identical(short, "auto")) short <- "arrow"
  parts <- ggchord_feature_geometry(
    shape, x0, x1, 1, width, 1, NULL,
    arrow_head_length = p$arrow_head_length %||% .04,
    arrow_head_width = p$arrow_head_width %||% 1,
    arrow_head_style = p$arrow_head_style %||% "shouldered",
    short_feature = short,
    draw_head = draw_head, draw_start_head = draw_start_head,
    bidirectional = bidirectional
  )
  lapply(parts, function(part) cbind(x = part$angle,
    y = center + part$radius - 1))
}

ggchord_collinear_position_fields <- function(row) {
  list(
    position_name = as.character(row$.position_name %||% "identity"),
    base_offset = as.numeric(row$.position_base_offset %||% 0),
    lane = as.integer(row$.feature_stack_lane %||% 0L),
    lane_offset = as.numeric(row$.position_lane_offset %||% 0),
    normal_offset = as.numeric(row$.normal_offset %||% 0)
  )
}

ggchord_collinear_gene_geometry <- function(p, layout) {
  input <- p$gene_data
  if (is.null(input) || !nrow(input)) return(data.frame())
  input$.source_row <- seq_len(nrow(input))
  input <- input[input$accver %in% p$seqs, , drop = FALSE]
  genes <- ggchord_expand_feature_render_runs(input)
  if (!isTRUE(p$draw_gene_geometry)) return(data.frame())
  ref <- layout$sequence_reference
  out <- list()
  for (i in seq_len(nrow(genes))) {
    row <- genes[i, , drop = FALSE]
    id <- as.character(row$accver)
    strand <- as.character(row$strand)
    width <- as.numeric(p$geneWidth[[id]][strand]) *
      as.numeric(row$.feature_width_factor %||% 1)
    if (!is.finite(width) || width <= 0) width <- .05
    fields <- ggchord_collinear_position_fields(row)
    backbone <- unname(ref$y[id])
    center <- backbone + fields$normal_offset
    shape <- as.character(row$.feature_shape %||% "arrow")
    biological <- as.character(row$.feature_biological_strand %||% strand)
    biological_x <- range(ggchord_collinear_x(c(row$start, row$end), id, ref))
    cut_left <- isTRUE(row$.truncated_left)
    cut_right <- isTRUE(row$.truncated_right)
    head_cut <- if (identical(strand, "-")) cut_left else cut_right
    start_head_cut <- if (identical(strand, "-")) cut_right else cut_left
    part_index <- 0L
    add_part <- function(xy, component = "polygon", boundary = NA_character_) {
      part_index <<- part_index + 1L
      out[[length(out) + 1L]] <<- data.frame(
        x = xy[, 1L], y = xy[, 2L], group = i * 100L + part_index,
        anno = as.character(row$anno), strand = strand,
        feature_shape = if (identical(component, "point")) "point" else shape,
        biological_strand = biological,
        feature_fill_explicit = as.character(row$feature_color %||% NA_character_),
        accver = id, sequence_y = backbone,
        position_name = fields$position_name, base_offset = fields$base_offset,
        lane = fields$lane, lane_offset = fields$lane_offset,
        normal_offset = fields$normal_offset,
        display_xmin = min(xy[, 1L]), display_xmax = max(xy[, 1L]),
        display_extended = min(xy[, 1L]) < biological_x[1L] - 1e-10 ||
          max(xy[, 1L]) > biological_x[2L] + 1e-10,
        display_extension_reason = if (shape %in% c("lollipop", "capped_line"))
          "symbol_geometry" else NA_character_,
        head_truncated = head_cut, start_head_truncated = start_head_cut,
        source_row = row$.source_row, ord = seq_len(nrow(xy)),
        boundary_linetype = boundary,
        boundary_draw_linetype = if (boundary %in% c("dotted", "dashed"))
          "solid" else boundary,
        .component = component, stringsAsFactors = FALSE
      )
    }
    if (row$start == row$end) {
      x <- ggchord_collinear_x(row$start, id, ref)
      half <- max(.014, min(.028, width * .36))
      add_part(cbind(rep(x, 2L), center + c(-half, half)), "point")
      next
    }
    pieces <- ggchord_feature_intervals(row$start, row$end, ref$lens[id], strand)
    for (piece in pieces) {
      x <- ggchord_collinear_x(c(piece$start, piece$end), id, ref)
      parts <- ggchord_collinear_feature_parts(
        shape, x[1L], x[2L], center, width, backbone, p,
        draw_head = piece$draw_head && isTRUE(row$.feature_draw_head %||% TRUE) && !head_cut,
        draw_start_head = isTRUE(row$.feature_draw_start_head %||%
          identical(biological, "+/-")) && !start_head_cut,
        bidirectional = identical(biological, "+/-")
      )
      for (xy in parts) add_part(xy)
    }
    boundaries <- if (".feature_boundaries" %in% names(row))
      row$.feature_boundaries[[1L]] else numeric()
    styles <- if (".feature_boundary_styles" %in% names(row))
      as.character(row$.feature_boundary_styles[[1L]]) else
        rep(NA_character_, length(boundaries))
    keep <- boundaries > min(row$start, row$end) &
      boundaries < max(row$start, row$end)
    for (j in which(keep)) {
      pattern <- if (identical(styles[j], "dotted")) {
        matrix(c(0, .13, .29, .42, .58, .71, .87, 1), ncol = 2L, byrow = TRUE)
      } else if (identical(styles[j], "dashed")) {
        matrix(c(0, .38, .62, 1), ncol = 2L, byrow = TRUE)
      } else matrix(c(0, 1), ncol = 2L)
      for (k in seq_len(nrow(pattern))) {
        x <- ggchord_collinear_x(boundaries[j], id, ref)
        add_part(cbind(rep(x, 2L), center - width / 2 +
          pattern[k, ] * width), "boundary", styles[j])
      }
    }
  }
  if (length(out)) ggchord_rbind_fill(out) else data.frame()
}

ggchord_collinear_region_geometry <- function(p, layout) {
  d <- p$region_data
  if (is.null(d) || !nrow(d)) return(data.frame())
  ggchord_require_columns(d, c("accver", "start", "end"), "geom_seq_region()")
  if (all(c("colour", "color") %in% names(d))) {
    ggchord_stop("geom_seq_region(): data may contain only one of colour and color")
  }
  ref <- layout$sequence_reference
  out <- lapply(seq_len(nrow(d)), function(i) {
    row <- d[i, , drop = FALSE]
    id <- as.character(row$accver)
    if (!id %in% p$seqs) return(NULL)
    positions <- as.numeric(c(row$start, row$end))
    if (any(!is.finite(positions))) return(NULL)
    positions <- pmax(1, pmin(ref$lens[id], positions))
    x <- sort(ggchord_collinear_x(positions, id, ref))
    center <- unname(ref$y[id]) +
      if (identical(p$region_side, "outside")) p$region_offset else -p$region_offset
    colour_col <- intersect(c("colour", "color"), names(d))
    fill <- if (length(colour_col) && !is.na(row[[colour_col]]))
      as.character(row[[colour_col]]) else p$region_fill
    data.frame(x = x[c(1, 2, 2, 1)],
      y = center + c(-1, -1, 1, 1) * p$region_width / 2,
      group = i, source_row = i, accver = id,
      zregionfill = fill, colour = p$region_color, alpha = p$region_alpha,
      label = as.character(row$label %||% NA_character_),
      category = as.character(row$category %||% NA_character_),
      stringsAsFactors = FALSE)
  })
  out <- Filter(Negate(is.null), out)
  if (length(out)) ggchord_rbind_fill(out) else data.frame()
}

ggchord_collinear_axis <- function(p, layout) {
  layout$axis_lines <- layout$axis_ticks <- data.frame()
  if (!isTRUE(p$show_axis)) return(layout)
  ref <- layout$sequence_reference
  units <- layout$text_units_per_inch
  lines <- ticks <- list()
  for (id in p$seqs) {
    gap <- unname(p$axisGap[id]) * units
    y <- unname(ref$y[id]) - gap
    lines[[id]] <- data.frame(
      x = ggchord_collinear_x(c(0, ref$lens[id]), id, ref),
      y = rep(y, 2L), accver = id
    )
    majors <- p$axis_breaks[[id]] %||% breakPointsFunc(ref$lens[id], p$axisMaj[id])
    major_labels <- p$axis_labels[[id]] %||% as.character(majors)
    minors <- p$axis_minor_breaks[[id]]
    if (is.null(minors)) {
      minors <- if (length(majors) > 1L) unlist(lapply(
        seq_len(length(majors) - 1L), function(j) {
          k <- as.integer(p$axisMin[id])
          if (!k) return(numeric())
          seq(majors[j], majors[j + 1L], length.out = k + 2L)[-c(1L, k + 2L)]
        })) else numeric()
    }
    pos <- c(majors, minors)
    major <- c(rep(TRUE, length(majors)), rep(FALSE, length(minors)))
    labels <- c(as.character(major_labels), rep(NA_character_, length(minors)))
    keep <- is.finite(pos) & pos >= 0 & pos <= ref$lens[id]
    pos <- pos[keep]; major <- major[keep]; labels <- labels[keep]
    if (!length(pos)) next
    len <- ifelse(major, p$axisMajLen[id], p$axisMinLen[id]) * units
    direction <- if (gap >= 0) -1 else 1
    angle <- p$axisLabelOrientation[[id]]
    angle <- if (is.character(angle)) {
      if (identical(tolower(angle), "perpendicular")) 90 else 0
    } else as.numeric(angle)
    label_y <- y + direction * (len + p$labelOffset[id] * units)
    ticks[[id]] <- data.frame(
      x0 = ggchord_collinear_x(pos, id, ref), y0 = y,
      x1 = ggchord_collinear_x(pos, id, ref), y1 = y + direction * len,
      label = labels, label_x = ggchord_collinear_x(pos, id, ref),
      label_y = label_y, size = unname(p$labelSize[id]),
      label_angle = angle, label_angle_relative = FALSE,
      label_along_axis = FALSE, label_gap = unname(p$labelOffset[id]) * units,
      label_hjust = .5, label_vjust = if (direction < 0) 1 else 0,
      is_major = major, is_origin = pos == 0, accver = id,
      stringsAsFactors = FALSE
    )
  }
  layout$axis_lines <- if (length(lines)) do.call(rbind, lines) else data.frame()
  layout$axis_ticks <- if (length(ticks)) do.call(rbind, ticks) else data.frame()
  if (isTRUE(p$axis_label_hide_overlaps) && nrow(layout$axis_ticks)) {
    layout$axis_ticks <- ggchord_hide_text_overlaps(layout$axis_ticks,
      ggchord_repel_points(layout$seq_arcs, layout$gene_polys,
        data.frame(x = numeric(), y = numeric()), data.frame(), show_axis = FALSE),
      units_per_inch = units)
  }
  layout
}

ggchord_collinear_seq_labels <- function(p, layout) {
  if (is.null(p$seq_label_text)) return(data.frame())
  ref <- layout$sequence_reference
  do.call(rbind, lapply(p$seqs, function(id) data.frame(
    text_x = min(ggchord_collinear_x(c(0, ref$lens[id]), id, ref)) - .04,
    text_y = unname(ref$y[id]),
    label = as.character(p$seq_label_text[id]),
    text_angle = if (identical(p$seq_label_orientation, "horizontal")) 0 else
      unname(p$seq_label_rotation[id]),
    size = unname(p$seq_label_size[id]),
    hjust = if (is.null(p$seq_label_hjust)) 1 else unname(p$seq_label_hjust[id]),
    vjust = if (is.null(p$seq_label_vjust)) .5 else unname(p$seq_label_vjust[id]),
    accver = id, stringsAsFactors = FALSE
  )))
}

ggchord_collinear_gene_labels <- function(p, layout) {
  d <- p$gene_data
  layout$gene_labels <- layout$gene_label_segments <- data.frame()
  layout$gene_label_layout <- "collinear"
  layout$gene_label_clip_units <- layout$text_units_per_inch
  if (!isTRUE(p$gene_label_show) || is.null(d) || !nrow(d)) return(layout)
  ref <- layout$sequence_reference
  labels <- lapply(seq_len(nrow(d)), function(i) {
    row <- d[i, , drop = FALSE]
    id <- as.character(row$accver); strand <- as.character(row$strand)
    if (!id %in% p$seqs) return(NULL)
    text <- as.character(if (is.null(p$gene_label_text)) row$anno else
      p$gene_label_text[i])
    if (is.na(text) || !nzchar(text)) return(NULL)
    fields <- ggchord_collinear_position_fields(row)
    width <- as.numeric(p$geneWidth[[id]][strand])
    center <- unname(ref$y[id]) + fields$normal_offset
    head_anchor <- identical(as.character(row$.feature_label_anchor %||% "midpoint"), "head")
    anchor_pos <- if (head_anchor) {
      if (strand == "-") row$start else row$end
    } else mean(c(row$start, row$end))
    circum <- p$geneLabelCircumOffset[[id]][strand]
    if (isTRUE(p$geneLabelCircumLimit[[id]][strand])) {
      limit <- abs(row$end - row$start) / ref$lens[id] / 2
      circum <- max(-limit, min(limit, circum))
    }
    pos <- max(0, min(ref$lens[id], anchor_pos + circum * ref$lens[id]))
    side <- if (identical(p$gene_label_side, "auto")) {
      if (fields$normal_offset != 0) sign(fields$normal_offset) else
        if (strand == "-") -1 else 1
    } else if (identical(p$gene_label_side, "inside")) -1 else 1
    feature_mode <- identical(p$gene_label_orientation, "feature") ||
      (isTRUE(p$gene_label_repel_layer) && identical(p$gene_label_layout, "feature"))
    angle <- if (isTRUE(p$gene_label_repel_layer)) 0 else
      (if (identical(p$gene_label_orientation, "radial")) 90 else 0) +
        p$geneLabelRotation[[id]][strand]
    x <- ggchord_collinear_x(pos, id, ref)
    offset <- p$geneLabelRadialOffset[[id]][strand]
    data.frame(text = text, text_x = x,
      text_y = center + side * (width / 2 + .025 + offset),
      text_angle = as.numeric(angle), hjust = .5, vjust = .5,
      size = p$gene_label_size, family = p$gene_label_family,
      fontface = p$gene_label_fontface, lineheight = p$gene_label_lineheight,
      accver = id, sequence_y = unname(ref$y[id]), group = i, source_row = i,
      anchor_position = anchor_pos, sequence_length = unname(ref$lens[id]),
      position_name = fields$position_name, base_offset = fields$base_offset,
      lane = fields$lane, lane_offset = fields$lane_offset,
      normal_offset = fields$normal_offset,
      anchor_x = ggchord_collinear_x(anchor_pos, id, ref), anchor_y = center,
      side_flipped = (side < 0) != (strand == "-"),
      feature_label_colour = "#202020",
      feature_label_fill = as.character(row$feature_label_fill %||% "#F2F3F4"),
      feature_label_orientation = if (angle == 90) "vertical" else "horizontal",
      feature_label_inside = FALSE, feature_label_mode = "adjacent",
      label_layout = "collinear", label_lane = 0L, .feature_width = width,
      label_side = side, label_overlap = p$gene_label_overlap,
      label_segment_linetype = rep_len(
        p$gene_label_segment_linetype %||% "auto", nrow(d))[i],
      label_segment_overlap = p$gene_label_segment_overlap,
      label_segment_overlap_alpha = p$gene_label_segment_overlap_alpha,
      .feature_mode = feature_mode, .label_side = side,
      .inside_colour = as.character(row$feature_label_colour %||% "#202020"),
      .feature_span = abs(row$end - row$start) * ref$x_scale,
      stringsAsFactors = FALSE)
  })
  labels <- Filter(Negate(is.null), labels)
  if (!length(labels)) return(layout)
  gl <- do.call(rbind, labels)
  gl$text <- ggchord_label_wrap_text(gl$text, p$gene_label_wrap)
  units <- layout$text_units_per_inch
  obstacles <- ggchord_text_obstacle_boxes(layout$seq_labels_df,
    layout$axis_ticks, p$show_axis, units_per_inch = units, box_padding = .015)
  fits_inside <- function(metric) gl$.feature_mode & metric$w < pmax(0,
    gl$.feature_span - (p$arrow_head_length %||% .04) - .02) &
    metric$h < gl$.feature_width * .9
  initial_inside <- fits_inside(ggchord_text_boxes(gl,
    units_per_inch = units, box_padding = .025))
  gl$text_y[initial_inside] <- gl$anchor_y[initial_inside]
  if (isTRUE(p$gene_label_repel_layer) && is.null(p$gene_label_wrap)) {
    fitted <- ggchord_fit_label_text(gl, fit = p$gene_label_fit,
      max_lines = p$gene_label_max_lines, units_per_inch = units,
      box_padding = .025, repel_boxes = obstacles)
    # A label that already fits the feature should not be wrapped because an
    # unused adjacent position would collide with the backbone's axis.
    gl$text[!initial_inside] <- fitted$text[!initial_inside]
  }
  metric <- ggchord_text_boxes(gl, units_per_inch = units, box_padding = .025)
  inside <- fits_inside(metric)
  gl$text_y[inside] <- gl$anchor_y[inside]
  gl$feature_label_inside[inside] <- TRUE
  gl$feature_label_mode[inside] <- "inside"
  gl$feature_label_colour[inside] <- gl$.inside_colour[inside]
  # Final collision placement is shared across all layer registries below.
  gl <- gl[!is.na(gl$text) & nzchar(gl$text), , drop = FALSE]
  if (nrow(gl) && isTRUE(p$gene_label_repel_layer)) {
    draw <- !gl$feature_label_inside & abs(gl$text_y - gl$anchor_y) > .015
    rows <- which(draw)
    if (length(rows)) {
      segment <- data.frame(x0 = gl$anchor_x[rows], y0 = gl$anchor_y[rows],
        x1 = gl$text_x[rows], y1 = gl$text_y[rows],
        group = gl$group[rows], source_row = gl$source_row[rows],
        accver = gl$accver[rows], alpha = 1, occluded = FALSE,
        linetype = if (identical(p$gene_label_segment_linetype, "auto"))
          ifelse(gl$side_flipped[rows], "dashed", "solid") else
            gl$label_segment_linetype[rows],
        stringsAsFactors = FALSE)
      # Keep group identities usable by the shared composite label renderer.
      group_map <- match(segment$group, gl$group)
      segment$group <- group_map
      gl$group <- seq_len(nrow(gl))
      segment <- ggchord_clip_segments_to_labels(segment, gl,
        units_per_inch = units, include_own = TRUE,
        overlap = p$gene_label_segment_overlap,
        overlap_alpha = p$gene_label_segment_overlap_alpha)
      layout$gene_label_segments <- segment
    }
  }
  gl$.feature_mode <- gl$.label_side <- gl$.feature_span <- gl$.inside_colour <- NULL
  layout$gene_labels <- gl
  layout
}

# Resolve independent label layers against one shared Cartesian obstacle set.
# Fixed labels retain their requested collision policy and become obstacles
# before automatic labels are considered, regardless of layer addition order.
ggchord_linear_line_hits <- function(a, b, x, y) {
  if (length(x) < 2L || max(x) < min(a[1], b[1]) ||
      min(x) > max(a[1], b[1]) || max(y) < min(a[2], b[2]) ||
      min(y) > max(a[2], b[2])) return(FALSE)
  x <- c(x, x[1L]); y <- c(y, y[1L])
  dx <- diff(x); dy <- diff(y); r <- b - a
  denominator <- r[1] * dy - r[2] * dx
  valid <- abs(denominator) > 1e-12
  t <- ((head(x, -1L) - a[1]) * dy - (head(y, -1L) - a[2]) * dx) /
    ifelse(valid, denominator, 1)
  u <- ((head(x, -1L) - a[1]) * r[2] - (head(y, -1L) - a[2]) * r[1]) /
    ifelse(valid, denominator, 1)
  any(valid & t > 1e-7 & t < 1 - 1e-7 & u >= 0 & u <= 1)
}

ggchord_linear_line_hits_box <- function(a, b, box) {
  xy <- ggchord_box_corners(box)
  ggchord_linear_line_hits(a, b, xy[, 1L], xy[, 2L]) ||
    ggchord_point_in_polygon(mean(c(a[1], b[1])), mean(c(a[2], b[2])),
      data.frame(x = xy[, 1L], y = xy[, 2L]))
}

ggchord_collinear_share_labels <- function(registry, layout) {
  units <- layout$text_units_per_inch
  ref <- layout$sequence_reference
  polygons <- list()
  fixed <- automatic <- list()
  diagnostics <- list()
  placed_routes <- list()
  for (id in names(registry)) {
    components <- registry[[id]]
    geometry <- components$gene_poly
    if (is.data.frame(geometry) && nrow(geometry)) {
      keep <- geometry$.component %in% "polygon"
      pieces <- split(geometry[keep, , drop = FALSE], geometry$group[keep])
      polygons <- c(polygons, unname(pieces))
    }
    for (component in intersect(c("gene_text", "gene_text_repel", "gene_label_repel"),
        names(components))) {
      d <- components[[component]]
      if (!is.data.frame(d) || !nrow(d)) next
      if (identical(component, "gene_label_repel")) {
        d <- d[d$.component %in% "text", , drop = FALSE]
      }
      if (!nrow(d)) next
      item <- list(layer_id = id, component = component, labels = d)
      if (identical(component, "gene_text")) fixed[[length(fixed) + 1L]] <- item else
        automatic[[length(automatic) + 1L]] <- item
    }
  }
  if (!length(fixed) && !length(automatic)) return(registry)
  sequence_labels <- Filter(function(d) is.data.frame(d) && nrow(d),
    lapply(registry, `[[`, "seq_label"))
  sequence_labels <- if (length(sequence_labels)) ggchord_rbind_fill(sequence_labels) else
    layout$seq_labels_df
  placed <- ggchord_text_obstacle_boxes(sequence_labels, layout$axis_ticks,
    layout$show_axis, units_per_inch = units, box_padding = .015)
  own_polygons <- function(row) vapply(polygons, function(polygon) {
      if (
          !all(c("accver", "start", "end") %in% names(row)) ||
          !all(c("accver", "start", "end") %in% names(polygon))) return(FALSE)
      first <- polygon[1L, , drop = FALSE]
      identical(as.character(first$accver), as.character(row$accver)) &&
        first$start == row$start && first$end == row$end &&
        abs((first$normal_offset %||% 0) - (row$normal_offset %||% 0)) < 1e-9
    }, logical(1L))
  geometry_hits <- function(box, row) {
    own <- own_polygons(row) & isTRUE(row$feature_label_inside)
    if (ggchord_boxes_hit_features(box, polygons[!own])) return(TRUE)
    any(vapply(layout$seqs, function(id) {
      if (isTRUE(row$feature_label_inside) && identical(id, as.character(row$accver)))
        return(FALSE)
      limits <- range(ggchord_collinear_x(c(0, ref$lens[id]), id, ref))
      box$xmin < limits[2] && box$xmax > limits[1] &&
        box$ymin < ref$y[id] + .006 && box$ymax > ref$y[id] - .006
    }, logical(1L)))
  }
  route_for <- function(row) {
    if (isTRUE(row$feature_label_inside) || abs(row$text_y - row$anchor_y) < .015) return(NULL)
    side <- sign(row$text_y - row$anchor_y)
    list(a = c(row$anchor_x, row$anchor_y + side * row$.feature_width / 2),
      b = c(row$text_x, row$text_y))
  }
  route_conflicts <- function(row, box, repel) {
    text_hits <- if (length(placed_routes)) sum(vapply(placed_routes, function(route) {
      ggchord_linear_line_hits_box(route$a, route$b, box)
    }, logical(1L))) else 0L
    route <- if (repel) route_for(row) else NULL
    if (is.null(route)) return(c(text_hits = text_hits, crossings = 0L, geometry_hits = 0L))
    if (!is.null(placed) && nrow(placed)) text_hits <- text_hits + sum(vapply(
      seq_len(nrow(placed)), function(j) ggchord_linear_line_hits_box(
        route$a, route$b, placed[j, , drop = FALSE]), logical(1L)))
    geometry <- sum(vapply(polygons[!own_polygons(row)], function(poly) {
      ggchord_linear_line_hits(route$a, route$b, poly$x, poly$y) ||
        ggchord_point_in_polygon(mean(c(route$a[1], route$b[1])),
          mean(c(route$a[2], route$b[2])), poly)
    }, logical(1L)))
    geometry <- geometry + sum(vapply(layout$seqs, function(id) {
      limits <- range(ggchord_collinear_x(c(0, ref$lens[id]), id, ref))
      route$a[1] >= limits[1] && route$a[1] <= limits[2] &&
        ref$y[id] > min(route$a[2], route$b[2]) + 1e-7 &&
        ref$y[id] < max(route$a[2], route$b[2]) - 1e-7
    }, logical(1L)))
    crosses <- if (length(placed_routes)) sum(vapply(placed_routes, function(old) {
      ggchord_segments_cross(route$a[1], route$a[2], route$b[1], route$b[2],
        old$a[1], old$a[2], old$b[1], old$b[2]) ||
        (abs(route$a[1] - route$b[1]) < 1e-9 &&
          abs(old$a[1] - old$b[1]) < 1e-9 && abs(route$a[1] - old$a[1]) < 1e-9 &&
          min(max(route$a[2], route$b[2]), max(old$a[2], old$b[2])) >
          max(min(route$a[2], route$b[2]), min(old$a[2], old$b[2])) + 1e-7)
    }, logical(1L))) else 0L
    c(text_hits = text_hits, crossings = crosses, geometry_hits = geometry)
  }
  append_box <- function(box) {
    placed <<- if (is.null(placed) || !nrow(placed)) box else
      ggchord_rbind_fill(list(placed, box))
  }
  solve <- function(item, repel) {
    gl <- item$labels
    visible <- rep(TRUE, nrow(gl))
    for (i in seq_len(nrow(gl))) {
      row <- gl[i, , drop = FALSE]
      policy <- as.character(row$label_overlap %||% "hide")
      movable <- repel || identical(policy, "nudge")
      side <- as.numeric(row$label_side %||%
        ifelse(row$text_y >= row$anchor_y, 1, -1))
      measured <- ggchord_text_boxes(row, units_per_inch = units, box_padding = .025)
      step <- measured$h + .025
      initial_y <- row$text_y
      found <- FALSE
      for (lane in if (movable) 0:64 else 0L) {
        row$text_y <- initial_y + side * lane * step
        if (lane > 0L) {
          row$feature_label_inside <- FALSE
          row$feature_label_mode <- "adjacent"
          row$feature_label_colour <- "#202020"
        }
        box <- ggchord_text_boxes(row, units_per_inch = units,
          box_padding = if (isTRUE(row$feature_label_inside)) .005 else .025)
        collision <- (!is.null(placed) && nrow(placed) &&
          any(ggchord_oriented_box_overlaps(box, placed))) || geometry_hits(box, row)
        route_hits <- route_conflicts(row, box, repel)
        collision <- collision || any(route_hits > 0)
        if (!collision || (!repel && identical(policy, "allow"))) {
          found <- TRUE
          break
        }
      }
      if (found) {
        row$label_lane <- as.integer(row$label_lane %||% 0L) + lane
        if ("x" %in% names(row)) row$x <- row$text_x
        if ("y" %in% names(row)) row$y <- row$text_y
        append_box(box)
        if (repel) {
          route <- route_for(row)
          if (!is.null(route)) placed_routes[[length(placed_routes) + 1L]] <<- route
        }
        gl[i, ] <- row
      } else visible[i] <- FALSE
      diagnostics[[length(diagnostics) + 1L]] <<- data.frame(
        layer_id = item$layer_id, source_row = row$source_row,
        status = if (found) "placed" else "hidden",
        reason = if (found) {
          if (lane > 0L) "shared_obstacle" else NA_character_
        } else if (any(route_hits > 0)) "leader_conflict" else
          if (movable) "no_collision_free_row" else "fixed_label_overlap",
        label_lane = as.integer(row$label_lane %||% 0L),
        leader_text_hits = unname(route_hits["text_hits"]),
        leader_crossings = unname(route_hits["crossings"]),
        leader_geometry_hits = unname(route_hits["geometry_hits"]),
        stringsAsFactors = FALSE)
    }
    gl <- gl[visible, , drop = FALSE]
    if (!repel) {
      registry[[item$layer_id]][[item$component]] <<- gl
      return(invisible(NULL))
    }
    original <- registry[[item$layer_id]][[item$component]]
    segments <- original[original$.component %in% "segment", , drop = FALSE]
    rebuilt <- list()
    for (i in seq_len(nrow(gl))) {
      if (isTRUE(gl$feature_label_inside[i]) ||
          abs(gl$text_y[i] - gl$anchor_y[i]) < .015) next
      group <- gl$group[i]
      template <- segments[segments$group %in% group, , drop = FALSE]
      base <- if (nrow(template)) template[1L, , drop = FALSE] else gl[i, , drop = FALSE]
      base$x0 <- base$x <- gl$anchor_x[i]
      base$y0 <- base$y <- gl$anchor_y[i] + sign(gl$text_y[i] - gl$anchor_y[i]) *
        gl$.feature_width[i] / 2
      base$x1 <- base$xend <- gl$text_x[i]
      base$y1 <- base$yend <- gl$text_y[i]
      base$group <- i
      base$source_row <- gl$source_row[i]
      base$alpha <- 1
      base$occluded <- FALSE
      requested_linetype <- gl$label_segment_linetype[i] %||% "auto"
      base$linetype <- if (!identical(requested_linetype, "auto")) {
        requested_linetype
      } else if (nrow(template)) template$linetype[1L] else
        if (isTRUE(gl$side_flipped[i])) "dashed" else "solid"
      rebuilt[[length(rebuilt) + 1L]] <- base
    }
    gl$group <- seq_len(nrow(gl))
    segments <- if (length(rebuilt)) ggchord_rbind_fill(rebuilt) else data.frame()
    if (nrow(segments)) {
      clip_text <- gl
      clip_text$feature_label_mode <- "linear"
      segments <- ggchord_clip_segments_to_labels(segments, clip_text,
        units_per_inch = units, include_own = TRUE,
        overlap = gl$label_segment_overlap[1L] %||% "fade",
        overlap_alpha = gl$label_segment_overlap_alpha[1L] %||% .18)
    }
    if (identical(item$component, "gene_label_repel")) {
      registry[[item$layer_id]][[item$component]] <<- ggchord_repel_geometry(list(
        gene_labels = gl, gene_label_segments = segments,
        gene_label_layout = "collinear", text_units_per_inch = units))
    } else {
      registry[[item$layer_id]][[item$component]] <<- gl
    }
    invisible(NULL)
  }
  for (item in fixed) solve(item, FALSE)
  for (item in automatic) solve(item, TRUE)
  attr(registry, "collinear_label_diagnostics") <- if (length(diagnostics))
    do.call(rbind, diagnostics) else data.frame()
  registry
}

ggchord_collinear_annotations <- function(p, layout) {
  layout$text_units_per_inch <- ggchord_device_units_per_inch(
    unlist(lapply(layout$seq_arcs, `[[`, "x"), use.names = FALSE),
    unlist(lapply(layout$seq_arcs, `[[`, "y"), use.names = FALSE))
  layout$gene_polys <- ggchord_collinear_gene_geometry(p, layout)
  layout$region_polys <- ggchord_collinear_region_geometry(p, layout)
  layout <- ggchord_collinear_axis(p, layout)
  layout$seq_labels_df <- ggchord_collinear_seq_labels(p, layout)
  layout <- ggchord_collinear_gene_labels(p, layout)
  genes <- p$gene_data
  layout$gene_pal <- if (!is.null(genes) && nrow(genes)) {
    if (identical(p$gene_color_scheme, "strand")) process_strand_colors(p$gene_colors) else
      process_manual_colors(p$gene_colors, unique(genes$anno), p$gene_order)
  } else character()
  layout$final_gene_order <- if (!is.null(genes))
    unique(c(p$gene_order, genes$anno)) else character()
  layout$gene_color_scheme <- p$gene_color_scheme
  layout$gene_label_show <- p$gene_label_show
  layout$gene_label_size <- p$gene_label_size
  layout$feature_shape_pal <- p$feature_shape_pal
  layout$feature_shape_order <- p$feature_shape_order
  layout$show_axis <- p$show_axis
  layout$axisLabelOrientation <- p$axisLabelOrientation
  layout
}
