# Fit the glyph envelopes actually drawn, rather than a whole-word width alone.
ggchord_feature_contains_label <- function(boxes, polygons, source_row) {
  own <- Filter(function(p) identical(p$source_row[1L], source_row), polygons)
  if (!length(own) || !nrow(boxes)) return(FALSE)
  all(vapply(seq_len(nrow(boxes)), function(i) {
    corners <- ggchord_box_corners(boxes[i, , drop = FALSE])
    probes <- rbind(corners, (corners + corners[c(2,3,4,1), ])/2,
      c(boxes$cx[i], boxes$cy[i]))
    any(vapply(own, function(p) all(vapply(seq_len(nrow(probes)), function(j)
      ggchord_point_in_polygon(probes[j,1], probes[j,2], p), logical(1))), logical(1)))
  }, logical(1)))
}

# Intersect the unchanged source-angle ray with its actual drawn contour.
# The side is a display attachment, never a new source position or interval.
ggchord_feature_label_attachment <- function(candidate, polygons, outward, units_per_inch) {
  if (isTRUE(candidate$repeat_count > 1L)) return(candidate)
  radius <- sqrt(candidate$anchor_x^2 + candidate$anchor_y^2)
  if (!is.finite(radius) || radius < 1e-8) return(candidate)
  direction <- c(candidate$anchor_x, candidate$anchor_y)/radius
  cross <- function(x,y) x[1]*y[2]-x[2]*y[1]
  distances <- numeric()
  for (p in polygons) {
    if (!identical(p$source_row[1L], candidate$source_row)) next
    xy <- cbind(p$x, p$y); xy <- rbind(xy, xy[1L, ])
    a <- xy[-nrow(xy), , drop = FALSE]
    edge <- xy[-1L, , drop = FALSE]-a
    den <- direction[1]*edge[,2]-direction[2]*edge[,1]
    at <- which(abs(den) > 1e-12)
    t <- (a[at,1]*edge[at,2]-a[at,2]*edge[at,1])/den[at]
    u <- (a[at,1]*direction[2]-a[at,2]*direction[1])/den[at]
    distances <- c(distances, t[u >= -1e-8 & u <= 1+1e-8 & t > 0])
  }
  if (!length(distances)) return(candidate)
  # A tiny visible-side clearance avoids a numerical re-entry at the boundary.
  edge <- if (outward) max(distances) else min(distances)
  edge <- edge + if (outward) .03*units_per_inch/25.4 else -.03*units_per_inch/25.4
  candidate$anchor_x <- edge*direction[1]
  candidate$anchor_y <- edge*direction[2]
  candidate
}

ggchord_feature_label_shifted <- function(candidate, height) {
  anchor_angle <- atan2(candidate$anchor_y, candidate$anchor_x)
  angle <- atan2(candidate$text_y, candidate$text_x)-anchor_angle
  abs(atan2(sin(angle), cos(angle))) * sqrt(candidate$text_x^2+candidate$text_y^2) > height*.5
}
# Preserve straight source-to-label association. Clip glyph envelopes instead
# of making a visibility-graph detour; transparent feature fills stay honest.
ggchord_simple_feature_leaders <- function(registry, layout) {
  boxes <- ggchord_annotation_text_obstacles(registry, layout)
  units <- layout$text_units_per_inch %||% .35
  for (owner in names(registry)) {
    geometry <- registry[[owner]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry) ||
        !"label_layout" %in% names(geometry) ||
        !any(geometry$label_layout %in% "feature")) next
    segments <- which(geometry$.component %in% "segment")
    if (!length(segments)) next
    parts <- list(geometry[-segments, , drop = FALSE])
    texts <- geometry[geometry$.component %in% "text", , drop = FALSE]
    for (source in unique(geometry$source_row[segments])) {
      row <- texts[texts$source_row %in% source, , drop = FALSE]
      if (nrow(row) != 1L || is.na(row$label) || !nzchar(row$label)) next
      if (row$feature_label_mode %in% "external") {
        # The shared perimeter has already ordered and attached this route.
        # Reconnecting its endpoints can undo that ordering and cross callouts.
        parts[[length(parts)+1L]] <- geometry[segments[
          geometry$source_row[segments] == source], , drop=FALSE]
        next
      }
      own <- boxes$owner == owner & !is.na(boxes$source_row) & boxes$source_row == source
      route <- data.frame(x0 = row$anchor_x, y0 = row$anchor_y,
        x1 = row$x, y1 = row$y)
      own_cuts <- ggchord_inner_leader_intervals(route, boxes[own, , drop = FALSE])
      if (nrow(own_cuts)) {
        stop <- max(0, min(own_cuts$start) - .01*units /
          max(sqrt((route$x1-route$x0)^2+(route$y1-route$y0)^2), 1e-8))
        route$x1 <- route$x0 + stop*(route$x1-route$x0)
        route$y1 <- route$y0 + stop*(route$y1-route$y0)
      }
      policy <- row$.leader_overlap %||% "fade"
      fade <- row$.leader_overlap_alpha %||% .18
      if ((route$x1-route$x0)^2+(route$y1-route$y0)^2 < 1e-16) next
      source_segments <- segments[geometry$source_row[segments] == source]
      base_alpha <- geometry$alpha[source_segments]
      base_alpha <- if (any(is.finite(base_alpha))) max(base_alpha, na.rm=TRUE) else NA_real_
      cuts <- if (identical(policy, "show")) data.frame(start=numeric(),end=numeric()) else
        ggchord_inner_leader_intervals(route, boxes[!own, , drop = FALSE])
      times <- sort(unique(c(0, 1, cuts$start, cuts$end)))
      times <- times[times >= 0 & times <= 1]
      for (i in seq_len(length(times)-1L)) {
        a <- times[i]; b <- times[i+1L]; middle <- (a+b)/2
        covered <- any(middle > cuts$start & middle < cuts$end)
        if (b-a < 1e-10 || (covered && identical(policy, "clip"))) next
        d <- geometry[segments[which(geometry$source_row[segments] == source)[1L]], , drop = FALSE]
        d$x <- route$x0 + a*(route$x1-route$x0)
        d$y <- route$y0 + a*(route$y1-route$y0)
        d$xend <- route$x0 + b*(route$x1-route$x0)
        d$yend <- route$y0 + b*(route$y1-route$y0)
        for (pair in list(c('x0','x'), c('y0','y'), c('x1','xend'), c('y1','yend')))
          if (pair[1L] %in% names(d)) d[[pair[1L]]] <- d[[pair[2L]]]
        d$occluded <- covered
        d$alpha <- if (covered && identical(policy, "fade")) base_alpha * fade else base_alpha
        d$leader_text_status <- 'direct_single_segment'
        parts[[length(parts)+1L]] <- d
      }
    }
    registry[[owner]]$gene_label_repel <- ggchord_rbind_fill(parts)
  }
  registry
}
