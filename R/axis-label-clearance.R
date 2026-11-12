# Automatically positioned circular strings may acquire bounded radial
# clearance. Tick anchors, explicit breaks and physical offsets stay put.
ggchord_clear_circular_axis_text <- function(text, layout) {
  required <- c('.axis_automatic_breaks', '.axis_offset_explicit', '.axis_allow_clearance')
  if (!all(required %in% names(text)) || !length(layout$layer_geometry)) return(text)
  at <- which(text$.axis_automatic_breaks %in% TRUE)
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
  if (!length(polygons)) return(text)
  units <- layout$text_units_per_inch %||% .35
  backbone_radii <- unlist(lapply(layout$layer_geometry, function(owner) {
    d <- owner$seq
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
  for (i in at) {
    original <- text[i, , drop = FALSE]
    boxes <- boxes_for(original)
    if (!feature_hit(boxes)) {
      text$.axis_clearance_status[i] <- 'clear'
      next
    }
    if (isTRUE(original$.axis_offset_explicit) || !isTRUE(original$.axis_allow_clearance)) {
      text$.axis_clearance_status[i] <- if (isTRUE(original$.axis_offset_explicit)) 'explicit_offset' else 'disabled'
      next
    }
    others <- ggchord_rbind_fill(list(static, boxes_for(text[-i, , drop = FALSE]))) %||% data.frame()
    radius <- sqrt(original$text_x^2 + original$text_y^2)
    if (!is.finite(radius) || radius <= 0) next
    height <- max(boxes$oh, na.rm = TRUE)
    step <- max(height/4, .15*units/25.4)
    limit <- max(2*height, .8*units/25.4)
    amounts <- seq(step, limit, by = step)
    # A narrow inner channel may not fit even at the font floor. One bounded
    # outer lane can clear the backbone without shrinking or deleting text.
    outer <- if (length(band)) band[2L]+height-radius else NA_real_
    if (is.finite(outer) && outer > limit && outer <= max(4*height, 3*units/25.4))
      amounts <- c(amounts, outer)
    text$.axis_clearance_status[i] <- 'feature_unresolved'
    for (amount in amounts) {
      candidate <- original
      candidate$text_x <- original$text_x*(radius+amount)/radius
      candidate$text_y <- original$text_y*(radius+amount)/radius
      moved <- boxes_for(candidate)
      if (feature_hit(moved) || backbone_hit(moved)) next
      hits_text <- nrow(others) > 0L && any(vapply(seq_len(nrow(moved)), function(k)
        any(ggchord_oriented_box_overlaps(moved[k, , drop = FALSE], others)), logical(1)))
      if (hits_text) next
      text$text_x[i] <- candidate$text_x
      text$text_y[i] <- candidate$text_y
      text$.axis_clearance_mm[i] <- amount*25.4/units
      text$.axis_clearance_status[i] <- 'feature_clearance'
      break
    }
  }
  text
}
