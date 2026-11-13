# Automatically positioned circular strings may acquire bounded radial
# clearance. Tick anchors, explicit breaks and physical offsets stay put.
ggchord_clear_circular_axis_text <- function(text, layout) {
  required <- c('.axis_automatic_breaks', '.axis_offset_explicit', '.axis_allow_clearance')
  if (!all(required %in% names(text)) || !length(layout$layer_geometry)) return(text)
  critical <- rep(FALSE,nrow(text))
  for(field in intersect(c('is_origin','is_endpoint','.axis_seam_combined'),names(text)))
    critical <- critical | text[[field]] %in% TRUE
  at <- which(text$.axis_automatic_breaks %in% TRUE)
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
  if (!length(polygons) && !nrow(roots)) return(text)
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
  backbone_hit <- function(boxes) length(band) > 0L && any(vapply(seq_len(nrow(boxes)), function(i) {
    b <- boxes[i, , drop = FALSE]
    # A rotated rectangle has a continuous radial range. Its closest point to
    # the centre is obtained by clamping the centre in its local frame.
    x <- -b$cx*cos(b$angle)-b$cy*sin(b$angle)
    y <- b$cx*sin(b$angle)-b$cy*cos(b$angle)
    nearest <- sqrt(max(abs(x)-b$ow/2, 0)^2+max(abs(y)-b$oh/2, 0)^2)
    farthest <- max(sqrt(rowSums(ggchord_box_corners(b)^2)))
    nearest <= band[2L] && farthest >= band[1L]
  }, logical(1)))
  boxes_for <- function(d) ggchord_arc_text_layout(d, units_per_inch = units,
    respect_justification = TRUE)$boxes
  feature_hit <- function(boxes) any(vapply(seq_len(nrow(boxes)), function(i)
    any(vapply(polygons, function(p) ggchord_box_overlaps_feature(
      boxes[i, , drop = FALSE], p, tol = 1e-8), logical(1))), logical(1)))
  static <- ggchord_annotation_text_obstacles(layout$layer_geometry, layout, include_axis = FALSE)
  text$.axis_clearance_mm <- 0
  text$.axis_clearance_status <- NA_character_
  text$.axis_label_side <- "following"
  text$.axis_tangent_gap_mm <- text$label_gap*25.4/units
  near_backbone <- function(boxes) {
    if (!length(band) || !nrow(boxes)) return(TRUE)
    radii <- sqrt(rowSums(ggchord_boxes_corners(boxes)^2))
    gap <- if (max(radii)<band[1L]) band[1L]-max(radii) else
      if (min(radii)>band[2L]) min(radii)-band[2L] else -Inf
    is.finite(gap) && gap>=.5*units/25.4 && gap<=2*units/25.4
  }
  for (i in at) {
    original <- text[i, , drop = FALSE]
    initial <- original
    active_roots <- if (critical[i]) roots[
      roots$x0==roots$x1 & roots$y0==roots$y1,,drop=FALSE] else roots
    boxes <- boxes_for(original)
    other_rows <- setdiff(seq_len(nrow(text)),i)
    if (critical[i]) other_rows <- other_rows[critical[other_rows] |
      !(text$.axis_automatic_breaks[other_rows] %in% TRUE)]
    others <- ggchord_rbind_fill(list(static,
      boxes_for(text[other_rows, , drop = FALSE]))) %||% data.frame()
    text_hit <- function(candidate) nrow(others) > 0L && any(vapply(
      seq_len(nrow(candidate)), function(k) any(ggchord_oriented_box_overlaps(
        candidate[k, , drop = FALSE], others)), logical(1)))
    # Use the narrow channel adjacent to the backbone first. Explicit offsets
    # and breaks keep their contract; glyph corners, rather than only a text
    # centre, determine clearance from the backbone.
    if (length(band) && !isTRUE(original$.axis_offset_explicit) &&
        isTRUE(original$.axis_allow_clearance) && nrow(boxes)) {
      radius <- sqrt(original$text_x^2+original$text_y^2)
      outer <- max(sqrt(rowSums(ggchord_boxes_corners(boxes)^2)))
      shift <- band[1L]-.65*units/25.4-outer
      candidate <- original
      candidate$text_x <- original$text_x*(radius+shift)/radius
      candidate$text_y <- original$text_y*(radius+shift)/radius
      close_boxes <- boxes_for(candidate)
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
        isTRUE(initial$.axis_allow_clearance) && exit_hit(boxes)) {
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
        moved <- boxes_for(candidate)
        shift <- band[1L]-option$radial*units/25.4 -
          max(sqrt(rowSums(ggchord_boxes_corners(moved)^2)))
        candidate$text_x <- point[1L]*(radius+shift)/radius
        candidate$text_y <- point[2L]*(radius+shift)/radius
        moved <- boxes_for(candidate)
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
    feature_collision <- feature_hit(boxes)
    anchor_collision <- !feature_collision && root_hit(boxes)
    if (!feature_collision && !anchor_collision &&
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
    # A narrow inner channel may not fit even at the font floor. One bounded
    # outer lane can clear the backbone without shrinking or deleting text.
    outer <- if (length(band)) band[2L]+height-radius else NA_real_
    if (is.finite(outer) && outer > limit && outer <= max(8*height, 3*units/25.4))
      amounts <- c(amounts, outer)
    amounts <- unique(amounts[order(abs(amounts), amounts < 0)])
    text$.axis_clearance_status[i] <- if (feature_collision) 'feature_unresolved' else 'anchor_unresolved'
    min_radius <- max(2*height, sum(boxes$ow)/(pi/2))
    for (amount in amounts) {
      if (radius+amount < min_radius) next
      candidate <- original
      candidate$text_x <- original$text_x*(radius+amount)/radius
      candidate$text_y <- original$text_y*(radius+amount)/radius
      moved <- boxes_for(candidate)
      if (feature_hit(moved) || backbone_hit(moved) || root_hit(moved) ||
          !near_backbone(moved)) next
      hits_text <- nrow(others) > 0L && any(vapply(seq_len(nrow(moved)), function(k)
        any(ggchord_oriented_box_overlaps(moved[k, , drop = FALSE], others)), logical(1)))
      if (hits_text) next
      text$text_x[i] <- candidate$text_x
      text$text_y[i] <- candidate$text_y
      text$.axis_clearance_mm[i] <- text$.axis_clearance_mm[i]+amount*25.4/units
      text$.axis_clearance_status[i] <- if (feature_collision) 'feature_clearance' else 'anchor_clearance'
      break
    }
    unresolved <- grepl("unresolved$",text$.axis_clearance_status[i])
    if (unresolved && critical[i]) {
      text$.axis_clearance_status[i] <- "endpoint_unresolved"
    } else if (unresolved && isTRUE(layout$axis_label_hide_overlaps)) {
      text$text[i] <- NA_character_
      text$.axis_clearance_status[i] <- "hidden_near_backbone"
    }
  }
  text
}
