# Automatically positioned circular strings may acquire bounded radial
# clearance. Tick anchors, explicit breaks and physical offsets stay put.
ggchord_axis_endpoints <- function(data) {
  fields <- intersect(c("is_origin", "is_endpoint", ".axis_seam_combined"), names(data))
  Reduce(function(a, field) a | data[[field]] %in% TRUE, fields,
    init = rep(FALSE, nrow(data)))
}

# Reserve endpoint text before routing annotations or thinning interior labels.
# Work in device inches so tangents and rotated text have the same metric,
# including on nonsquare panels. Only text with a competing endpoint moves.
ggchord_space_axis_endpoints <- function(layout) {
  ticks <- layout$axis_ticks
  if (!isTRUE(layout$axis_label_hide_overlaps) || !NROW(ticks)) return(ticks)
  rows <- which(ggchord_axis_endpoints(ticks) & !is.na(ticks$label) & nzchar(ticks$label))
  if (isTRUE(layout$circular)) {
    curved <- ticks$label_along_axis %in% TRUE & !(ticks$.axis_plotmath %in% TRUE)
    rows <- rows[!curved[rows]]
  }
  if (length(rows) < 2L) return(ticks)
  metric <- ggchord_annotation_metric(layout) * 25.4
  labels <- ticks[rows, , drop = FALSE]
  labels$label_x <- labels$label_x / metric[1L]
  labels$label_y <- labels$label_y / metric[2L]
  padding <- .3 / 25.4
  boxes <- ggchord_text_boxes(labels, "label_x", "label_y", "label", "label_angle",
    hjust_col = "label_hjust", vjust_col = "label_vjust", units_per_inch = 1,
    box_padding = padding)
  crowded <- which(vapply(seq_len(nrow(boxes)), function(i)
    any(ggchord_oriented_box_overlaps(boxes, boxes, i = i)[-i]), logical(1)))
  if (!length(crowded)) return(ticks)

  paths <- if (NROW(layout$axis_lines)) split(layout$axis_lines, layout$axis_lines$accver)
  reserved <- setdiff(seq_len(nrow(boxes)), crowded)
  occupied <- ggchord_annotation_text_obstacles(layout$layer_geometry, layout,
    include_axis = FALSE, scale = metric)
  polygons <- ggchord_annotation_polygons(layout$layer_geometry, metric)
  bounds <- ggchord_polygon_bounds(polygons)
  for (i in crowded) {
    id <- as.character(labels$accver[i])
    path <- paths[[id]] %||% layout$seq_arcs[[id]]
    if (NROW(path) < 2L) next
    n <- nrow(path)
    edges <- c(1L, n)
    distance <- ((path$x[edges] - labels$x0[i]) / metric[1L])^2 +
      ((path$y[edges] - labels$y0[i]) / metric[2L])^2
    end <- which.min(distance)
    neighbour <- c(2L, n - 1L)[end]
    tangent <- c(path$x[neighbour] - path$x[edges[end]],
      path$y[neighbour] - path$y[edges[end]]) / metric
    magnitude <- sqrt(sum(tangent^2))
    if (!is.finite(magnitude) || magnitude <= 0) next
    tangent <- tangent / magnitude
    box <- boxes[i, ]
    # Place the nearest edge just inside its own sequence. At a crowded
    # row or a very short sequence, try the other side of the same tick.
    corners <- ggchord_boxes_corners(box)
    projection <- (corners[, 1L] - box$x) * tangent[1L] +
      (corners[, 2L] - box$y) * tangent[2L]
    shifts <- c(max(0, -min(projection)), 0, min(0, -max(projection)))
    candidates <- ggchord_project_text_boxes(box$x + tangent[1L] * shifts,
      box$y + tangent[2L] * shifts, box$w, box$h, box$angle,
      labels$label_hjust[i], labels$label_vjust[i], 1, padding)
    clear <- which(vapply(seq_along(shifts), function(j)
      !any(ggchord_oriented_box_overlaps(candidates, boxes[reserved, , drop = FALSE], i = j)) &&
        !any(ggchord_oriented_box_overlaps(candidates, occupied, i = j)) &&
        !ggchord_boxes_hit_features(candidates[j, , drop = FALSE], polygons, bounds = bounds),
      logical(1)))
    if (length(clear)) {
      chosen <- candidates[clear[1L], , drop = FALSE]
      boxes[i, ] <- chosen
      ticks$label_x[rows[i]] <- chosen$x * metric[1L]
      ticks$label_y[rows[i]] <- chosen$y * metric[2L]
    }
    reserved <- c(reserved, i)
  }
  ticks
}

ggchord_clear_projected_axis_text <- function(text, layout) {
  if (!isTRUE(layout$axis_label_hide_overlaps) || !length(layout$layer_geometry))
    return(text)
  scale <- ggchord_annotation_metric(layout) * 25.4
  occupied <- ggchord_annotation_text_obstacles(layout$layer_geometry, layout,
    include_axis = FALSE, scale = scale)
  leaders <- ggchord_rbind_fill(lapply(layout$layer_geometry, function(owner) {
    d <- owner$domain_label
    if (!NROW(d)) return(NULL)
    ggchord_annotation_segments(d[d$.component %in% "leader", , drop = FALSE], scale)
  })) %||% data.frame()
  polygons <- ggchord_annotation_polygons(layout$layer_geometry, scale)
  bounds <- ggchord_polygon_bounds(polygons)
  boxes <- ggchord_text_boxes(transform(text, x = x / scale[1L], y = y / scale[2L],
    text = label), x_col = "x", y_col = "y", angle_col = "angle",
    units_per_inch = 1, box_padding = .3 / 25.4)
  endpoints <- ggchord_axis_endpoints(text)
  accepted <- which(endpoints)
  hidden <- integer()
  for (i in which(!endpoints)) {
    box <- boxes[i, , drop = FALSE]
    if (ggchord_boxes_overlap_any(box, occupied) ||
        ggchord_boxes_overlap_any(box, boxes[accepted, , drop = FALSE]) ||
        ggchord_inner_leader_hits(leaders, box) ||
        ggchord_boxes_hit_features(box, polygons, tol = .1 / 25.4, bounds = bounds)) {
      hidden <- c(hidden, i)
    } else {
      accepted <- c(accepted, i)
    }
  }
  if (length(hidden)) text <- text[-hidden, , drop = FALSE]
  text
}

ggchord_clear_circular_axis_text <- function(text, layout) {
  required <- c('.axis_automatic_breaks', '.axis_offset_explicit', '.axis_allow_clearance')
  if (!all(required %in% names(text)) || !length(layout$layer_geometry)) return(text)
  critical <- ggchord_axis_endpoints(text)
  at <- which(text$.axis_automatic_breaks %in% TRUE & text$.axis_allow_clearance %in% TRUE)
  at <- at[order(!critical[at],at)]
  if (!length(at)) return(text)
  polygons <- list()
  for (owner in layout$layer_geometry) {
    d <- owner$gene_poly
    if (!is.data.frame(d) || !nrow(d)) next
    keep <- if ('.component' %in% names(d)) is.na(d$.component) | d$.component %in% 'polygon' else rep(TRUE, nrow(d))
    if ('feature_shape' %in% names(d)) keep <- keep & !d$feature_shape %in% 'capped_line'
    d <- d[keep & is.finite(d$x) & is.finite(d$y), , drop = FALSE]
    if (!nrow(d)) next
    parts <- split(d, d$group, drop = TRUE)
    polygons <- c(polygons, Filter(function(p) nrow(p) >= 3L, parts))
  }
  units <- layout$text_units_per_inch %||% .35
  # Axis placement precedes leader routing. Depending on the current leader
  # path here would move text again after its obstacle was already routed.
  roots <- ggchord_annotation_leader_roots(layout$layer_geometry)
  # Bare gapped circles can have two colliding endpoints even when there
  # are no annotation obstacles. Check the actual curved glyphs first.
  endpoints <- which(critical)
  metrics <- lapply(seq_len(nrow(text)), function(i) new.env(parent = emptyenv()))
  boxes_for <- function(d, i) ggchord_arc_text_layout(d, units_per_inch = units,
    box_padding = if (length(endpoints) > 1L && critical[i]) .3 / 25.4 else 0,
    respect_justification = TRUE, line_metrics = metrics[[i]], boxes_only = TRUE)$boxes
  if (!length(polygons) && !nrow(roots) && !any(text$.axis_seam_combined %in% TRUE)) {
    if (length(endpoints) < 2L) return(text)
    endpoint_boxes <- lapply(endpoints, function(i) boxes_for(text[i, , drop = FALSE], i))
    if (!any(vapply(seq_along(endpoints)[-1L], function(j) ggchord_boxes_overlap_any(
        endpoint_boxes[[j]], do.call(rbind, endpoint_boxes[seq_len(j - 1L)])), logical(1))))
      return(text)
  }
  active_roots <- roots
  root_hit <- function(boxes) {
    if (!nrow(active_roots) || !nrow(boxes)) return(FALSE)
    boxes$ow <- boxes$ow+2*.3*units/25.4
    boxes$oh <- boxes$oh+2*.3*units/25.4
    ggchord_inner_leader_hits(active_roots, boxes)
  }
  backbone_radii <- unlist(lapply(layout$layer_geometry, function(owner) {
    d <- owner[['seq']]
    if (!is.data.frame(d) || !all(c('x', 'y') %in% names(d))) return(numeric())
    sqrt(d$x^2+d$y^2)
  }), use.names = FALSE)
  backbone_radii <- backbone_radii[is.finite(backbone_radii)]
  band <- if (length(backbone_radii)) range(backbone_radii) + c(-1, 1)*.15*units/25.4 else numeric()
  backbone_hit <- function(boxes) {
    if (!length(band) || !nrow(boxes)) return(FALSE)
    # Clamp the origin in each glyph's local frame for the nearest radius.
    ca <- cos(boxes$angle)
    sa <- sin(boxes$angle)
    x <- -boxes$cx * ca - boxes$cy * sa
    y <- boxes$cx * sa - boxes$cy * ca
    nearest <- sqrt(pmax(abs(x) - boxes$ow / 2, 0)^2 +
      pmax(abs(y) - boxes$oh / 2, 0)^2)
    corners <- matrix(sqrt(rowSums(ggchord_boxes_corners(boxes)^2)), nrow = 4L)
    farthest <- pmax(corners[1, ], corners[2, ], corners[3, ], corners[4, ])
    any(nearest <= band[2L] & farthest >= band[1L])
  }
  polygon_bounds <- ggchord_polygon_bounds(polygons)
  feature_hit <- function(boxes) ggchord_boxes_hit_features(boxes, polygons, tol = 1e-8, bounds = polygon_bounds)
  static <- ggchord_annotation_text_obstacles(layout$layer_geometry, layout, include_axis = FALSE)
  text$.axis_clearance_mm <- 0
  text$.axis_clearance_status <- NA_character_
  text$.axis_label_side <- "following"
  text$.axis_tangent_gap_mm <- text$label_gap*25.4/units
  inward <- function(row) !identical(row$.axis_inward, FALSE)
  near_backbone <- function(boxes) {
    if (!length(band) || !nrow(boxes)) return(TRUE)
    radii <- sqrt(rowSums(ggchord_boxes_corners(boxes)^2))
    gap <- if (inward(original)) band[1L] - max(radii) else
      min(radii) - band[2L]
    is.finite(gap) && gap >= .35 * units / 25.4
  }
  # Glyphs on the other ticks keep their measured boxes until that tick moves.
  tick_boxes <- lapply(seq_len(nrow(text)), function(i) boxes_for(text[i, , drop = FALSE], i))
  angular_span <- function(boxes, center) {
    corners <- ggchord_boxes_corners(boxes)
    angle <- atan2(corners[, 2L], corners[, 1L]) - center
    range(atan2(sin(angle), cos(angle)))
  }
  seam <- which(text$.axis_seam_combined %in% TRUE)
  seam_conflict <- function(boxes, i) {
    if (!isTRUE(layout$axis_label_hide_overlaps) || !nrow(boxes)) return(FALSE)
    for (j in seam[text$accver[seam] == text$accver[i]]) {
      if (!nrow(tick_boxes[[j]])) next
      center <- atan2(text$text_y[j], text$text_x[j])
      a <- angular_span(boxes, center)
      b <- angular_span(tick_boxes[[j]], center)
      # Measured glyphs reserve one angular channel around the origin.
      # Moving an ordinary label deeper does not make it a separate label slot.
      margin <- .45 * units / (25.4 * max(.1, sqrt(text$text_x[j]^2 + text$text_y[j]^2)))
      if (diff(a) < pi && a[1L] < b[2L] + margin && a[2L] > b[1L] - margin)
        return(TRUE)
    }
    FALSE
  }
  for (i in at) {
    original <- text[i, , drop = FALSE]
    initial <- original
    active_roots <- if (critical[i]) roots[
      roots$x0==roots$x1 & roots$y0==roots$y1,,drop=FALSE] else roots
    boxes <- tick_boxes[[i]]
    if (!critical[i] && seam_conflict(boxes, i)) {
      text$text[i] <- NA_character_
      tick_boxes[[i]] <- boxes_for(text[i, , drop = FALSE], i)
      text$.axis_clearance_status[i] <- "hidden_at_seam"
      next
    }
    other_rows <- setdiff(seq_len(nrow(text)),i)
    if (critical[i]) other_rows <- other_rows[critical[other_rows] |
      !(text$.axis_automatic_breaks[other_rows] %in% TRUE)]
    others <- ggchord_rbind_fill(c(list(static), tick_boxes[other_rows])) %||% data.frame()
    text_hit <- function(candidate) ggchord_boxes_overlap_any(candidate, others)
    # Use the narrow channel adjacent to the backbone first. Explicit offsets
    # and breaks keep their contract; glyph corners, rather than only a text
    # centre, determine clearance from the backbone.
    if (length(band) && !isTRUE(original$.axis_offset_explicit) &&
        isTRUE(original$.axis_allow_clearance) && nrow(boxes)) {
      radius <- sqrt(original$text_x^2+original$text_y^2)
      envelope <- sqrt(rowSums(ggchord_boxes_corners(boxes)^2))
      shift <- if (inward(original)) band[1L] - .65 * units / 25.4 - max(envelope) else
        band[2L] + .65 * units / 25.4 - min(envelope)
      candidate <- original
      candidate$text_x <- original$text_x*(radius+shift)/radius
      candidate$text_y <- original$text_y*(radius+shift)/radius
      close_boxes <- boxes_for(candidate, i)
      if (!feature_hit(close_boxes) && !backbone_hit(close_boxes) && !root_hit(close_boxes) &&
          !text_hit(close_boxes)) {
        original <- candidate
        text$text_x[i] <- candidate$text_x; text$text_y[i] <- candidate$text_y
        text$.axis_clearance_mm[i] <- shift*25.4/units
        boxes <- close_boxes
      }
    }
    # A seam string can occupy either readable side of its unchanged tick.
    # Reserve fixed source exit rays here, before the outer queue is routed;
    # using mutable elbows would create a text/leader positioning cycle.
    # Only automatic, unlocked endpoints may acquire another nearby slot.
    exit_hit <- function(candidate) {
      if (!nrow(roots)) return(FALSE)
      candidate$ow <- candidate$ow + .6*units/25.4
      candidate$oh <- candidate$oh + .6*units/25.4
      ggchord_inner_leader_hits(roots, candidate)
    }
    if (critical[i] && length(band) && !isTRUE(initial$.axis_offset_explicit) &&
        (exit_hit(boxes) || feature_hit(boxes) || text_hit(boxes))) {
      direction <- c(cos(initial$text_angle*pi/180), sin(initial$text_angle*pi/180))
      tick <- c(initial$text_x, initial$text_y) - direction*initial$label_gap
      gaps <- unique(c(initial$label_gap*25.4/units, .8, 1.2, 1.6, 2, 2.5, 3, 3.5, 4))
      slots <- expand.grid(side = c(0, 1), gap = gaps, radial = c(.65, 1, 1.5, 2))
      slots$cost <- abs(slots$gap-initial$label_gap*25.4/units) +
        slots$radial-.65 + .2*slots$side
      slots <- slots[order(slots$cost, slots$side, slots$gap), , drop = FALSE]
      for (slot in seq_len(nrow(slots))) {
        option <- slots[slot, ]
        candidate <- initial
        candidate$hjust <- option$side
        point <- tick + direction*(1-2*option$side)*option$gap*units/25.4
        radius <- sqrt(sum(point^2))
        if (!is.finite(radius) || radius <= 0) next
        candidate$text_x <- point[1L]; candidate$text_y <- point[2L]
        moved <- boxes_for(candidate, i)
        envelope <- sqrt(rowSums(ggchord_boxes_corners(moved)^2))
        shift <- if (inward(original)) band[1L] - option$radial * units / 25.4 - max(envelope) else
          band[2L] + option$radial * units / 25.4 - min(envelope)
        candidate$text_x <- point[1L]*(radius+shift)/radius
        candidate$text_y <- point[2L]*(radius+shift)/radius
        moved <- boxes_for(candidate, i)
        if (feature_hit(moved) || backbone_hit(moved) || exit_hit(moved) ||
            text_hit(moved) || !near_backbone(moved)) next
        original <- candidate
        text$text_x[i] <- candidate$text_x; text$text_y[i] <- candidate$text_y
        text$hjust[i] <- candidate$hjust
        text$.axis_clearance_mm[i] <- (sqrt(candidate$text_x^2+candidate$text_y^2)-
          sqrt(initial$text_x^2+initial$text_y^2))*25.4/units
        text$.axis_label_side[i] <- if (option$side == 1) "preceding" else "following"
        text$.axis_tangent_gap_mm[i] <- option$gap
        boxes <- moved
        break
      }
    }
    tick_boxes[[i]] <- boxes
    feature_collision <- feature_hit(boxes)
    anchor_collision <- !feature_collision && root_hit(boxes)
    if (!feature_collision && !anchor_collision && !text_hit(boxes) &&
        (isTRUE(original$.axis_offset_explicit) || !isTRUE(original$.axis_allow_clearance) ||
          near_backbone(boxes))) {
      text$.axis_clearance_status[i] <- 'clear'
      next
    }
    if (isTRUE(original$.axis_offset_explicit) || !isTRUE(original$.axis_allow_clearance)) {
      text$.axis_clearance_status[i] <- if (isTRUE(original$.axis_offset_explicit)) 'explicit_offset' else 'disabled'
      next
    }
    radius <- sqrt(original$text_x^2 + original$text_y^2)
    if (!is.finite(radius) || radius <= 0) next
    height <- max(boxes$oh, na.rm = TRUE)
    step <- max(height/4, .15*units/25.4)
    limit <- max(2*height, .8*units/25.4)
    amounts <- seq(step, limit, by = step)
    amounts <- c(amounts, -seq(step, max(2*height, 2*units/25.4), by = step))
    # Protect endpoints with a deeper lane on the same side.
    if (critical[i] && length(band)) {
      direction <- if (inward(original)) -1 else 1
      amounts <- c(amounts, direction * c(3, 4.5, 6) * units / 25.4)
    }
    amounts <- unique(amounts[order(abs(amounts), amounts < 0)])
    text$.axis_clearance_status[i] <- if (feature_collision) 'feature_unresolved' else 'anchor_unresolved'
    min_radius <- max(2*height, sum(boxes$ow)/(pi/2))
    for (amount in amounts) {
      if (radius+amount < min_radius) next
      candidate <- original
      candidate$text_x <- original$text_x*(radius+amount)/radius
      candidate$text_y <- original$text_y*(radius+amount)/radius
      moved <- boxes_for(candidate, i)
      if (feature_hit(moved) || backbone_hit(moved) || root_hit(moved) ||
          (!critical[i] && seam_conflict(moved, i)) ||
          !near_backbone(moved)) next
      hits_text <- nrow(others) > 0L && any(vapply(seq_len(nrow(moved)), function(k)
        any(ggchord_oriented_box_overlaps(moved[k, , drop = FALSE], others)), logical(1)))
      if (hits_text) next
      tick_boxes[[i]] <- moved
      text$text_x[i] <- candidate$text_x
      text$text_y[i] <- candidate$text_y
      text$.axis_clearance_mm[i] <- text$.axis_clearance_mm[i]+amount*25.4/units
      text$.axis_clearance_status[i] <- if (feature_collision) 'feature_clearance' else 'anchor_clearance'
      break
    }
    unresolved <- grepl("unresolved$",text$.axis_clearance_status[i])
    if (unresolved && critical[i]) {
      # Keep even a crowded endpoint on the requested side of the backbone.
      if (!near_backbone(boxes) && length(band)) {
        candidate <- original
        corner_radius <- sqrt(rowSums(ggchord_boxes_corners(boxes)^2))
        target <- if (inward(original)) band[1L] - max(corner_radius) - height else
          band[2L] - min(corner_radius) + height
        candidate$text_x <- original$text_x * (radius + target) / radius
        candidate$text_y <- original$text_y * (radius + target) / radius
        text$text_x[i] <- candidate$text_x; text$text_y[i] <- candidate$text_y
        tick_boxes[[i]] <- boxes_for(candidate, i)
        text$.axis_clearance_mm[i] <- text$.axis_clearance_mm[i] + target*25.4/units
      }
      text$.axis_clearance_status[i] <- "endpoint_unresolved"
    } else if (unresolved && isTRUE(layout$axis_label_hide_overlaps)) {
      text$text[i] <- NA_character_
      tick_boxes[[i]] <- boxes_for(text[i, , drop = FALSE], i)
      text$.axis_clearance_status[i] <- "hidden_near_backbone"
    }
  }
  text
}
