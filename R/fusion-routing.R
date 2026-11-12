ggchord_route_inflate <- function(boxes, padding) {
  for (name in intersect(c("xmin", "ymin"), names(boxes))) boxes[[name]] <- boxes[[name]]-padding
  for (name in intersect(c("xmax", "ymax"), names(boxes))) boxes[[name]] <- boxes[[name]]+padding
  boxes
}

# Deterministic physical routes. These adjust drawing, never source intervals.
ggchord_route_segment_hits <- function(a, b, boxes, tolerance = 1e-7) {
  if (!nrow(boxes)) return(logical())
  enter <- rep(0, nrow(boxes)); leave <- rep(1, nrow(boxes))
  for (axis in 1:2) {
    lo <- boxes[[if (axis == 1L) 'xmin' else 'ymin']] + tolerance
    hi <- boxes[[if (axis == 1L) 'xmax' else 'ymax']] - tolerance
    delta <- b[axis] - a[axis]
    if (abs(delta) < 1e-12) {
      leave[a[axis] <= lo | a[axis] >= hi] <- -Inf
    } else {
      first <- (lo - a[axis])/delta; last <- (hi - a[axis])/delta
      enter <- pmax(enter, pmin(first, last)); leave <- pmin(leave, pmax(first, last))
    }
  }
  enter < leave & leave > 0 & enter < 1
}

ggchord_route_polygon_hits <- function(path, boxes) {
  if (!nrow(boxes)) return(logical())
  path <- as.matrix(path)
  hits <- rep(FALSE, nrow(boxes))
  for (i in seq_len(nrow(path))) hits <- hits | ggchord_route_segment_hits(
    path[i, ], path[if (i == nrow(path)) 1L else i + 1L, ], boxes)
  polygon <- data.frame(x = path[, 1], y = path[, 2])
  hits | vapply(seq_len(nrow(boxes)), function(i) ggchord_point_in_polygon(
    mean(c(boxes$xmin[i], boxes$xmax[i])), mean(c(boxes$ymin[i], boxes$ymax[i])), polygon), logical(1))
}

ggchord_route_obstacles <- function(obstacles, endpoints) {
  # Endpoint containers must be entered; their measured text/guide boxes remain obstacles.
  obstacles[!obstacles$kind %in% 'container' | !obstacles$panel_id %in% endpoints, , drop = FALSE]
}

ggchord_route_shortest <- function(start, end, boxes, canvas, clearance = .6) {
  if (any(c(start, end) < 0) || any(c(start[1], end[1]) > canvas[1]) ||
      any(c(start[2], end[2]) > canvas[2])) return(NULL)
  padded <- ggchord_route_inflate(boxes, clearance)
  if (!any(ggchord_route_segment_hits(start, end, padded))) return(rbind(start, end))
  corners <- do.call(rbind, lapply(seq_len(nrow(padded)), function(i)
    as.matrix(expand.grid(x = c(padded$xmin[i], padded$xmax[i]),
      y = c(padded$ymin[i], padded$ymax[i])))))
  corners <- corners[corners[, 1] >= clearance & corners[, 1] <= canvas[1] - clearance &
    corners[, 2] >= clearance & corners[, 2] <= canvas[2] - clearance, , drop = FALSE]
  nodes <- unique(rbind(start, end, corners))
  distance <- rep(Inf, nrow(nodes)); distance[1L] <- 0
  previous <- rep(NA_integer_, nrow(nodes)); visited <- rep(FALSE, nrow(nodes))
  for (step in seq_len(nrow(nodes))) {
    available <- which(!visited)
    at <- available[which.min(distance[available])]
    if (!is.finite(distance[at])) return(NULL)
    if (at == 2L) break
    visited[at] <- TRUE
    for (to in which(!visited)) {
      cost <- sqrt(sum((nodes[to, ] - nodes[at, ])^2))
      if (distance[at] + cost >= distance[to] ||
          any(ggchord_route_segment_hits(nodes[at, ], nodes[to, ], padded))) next
      distance[to] <- distance[at] + cost; previous[to] <- at
    }
  }
  if (!is.finite(distance[2L])) return(NULL)
  route <- 2L
  while (route[1L] != 1L) route <- c(previous[route[1L]], route)
  nodes[route, , drop = FALSE]
}

ggchord_route_ribbon <- function(q, s, direct, obstacles, canvas, route, width, corridor = NULL, tip = "none") {
  visible <- function(path) ggchord_ribbon_visible_polygon(path, nrow(q), nrow(s), tip)
  hits <- ggchord_route_polygon_hits(visible(direct), obstacles)
  outside <- any(direct[, 1] < 0 | direct[, 1] > canvas[1] | direct[, 2] < 0 | direct[, 2] > canvas[2])
  if (route == 'direct' || (route == 'avoid' && !any(hits) && !outside)) return(list(path = direct,
    status = if (outside) 'direct_outside_canvas' else if (any(hits)) 'direct_collision' else 'direct_clear',
    collisions = sum(hits)))
  a <- as.matrix(q[c('x', 'y')]); b <- as.matrix(s[c('x', 'y')])
  qc <- colMeans(a); sc <- colMeans(b)
  normals <- attr(direct, 'endpoint_normals')
  lead <- max(1, width)
  start <- qc + normals['q', ]*lead*2; end <- sc + normals['s', ]*lead*2
  if (sum((start-end)^2) < 1e-12)
    ggchord_stop('Fusion: no distinct routing gates at these endpoint fronts')
  for (padding in c(.6, 1.2, 2.4)) {
    line <- corridor %||% ggchord_route_shortest(start, end, obstacles, canvas, clearance = width/2 + padding)
    if (is.null(line)) next
    if (nrow(line) == 2L) {
      # A clear centreline can still have a wide colliding direct ribbon.
      line <- rbind(start, start + (end-start)*.15, start + (end-start)*.85, end)
    }
    delta <- line[-1L, , drop = FALSE] - line[-nrow(line), , drop = FALSE]
    normal <- cbind(-delta[, 2], delta[, 1])/sqrt(rowSums(delta^2))
    midnormal <- rbind(normal[1L, ], (normal[-nrow(normal), , drop = FALSE] +
      normal[-1L, , drop = FALSE])/2, normal[nrow(normal), ])
    midnormal <- midnormal/pmax(sqrt(rowSums(midnormal^2)), .2)
    if (sum(midnormal[1L, ]*(a[nrow(a), ]-a[1L, ])) < 0) midnormal <- -midnormal
    edge1 <- line - midnormal*width/2; edge2 <- line + midnormal*width/2
    # Short endpoint departures retain the original physical front normals.
    edge1 <- rbind(a[1L, ], a[1L, ]+normals['q', ]*lead, edge1,
      b[1L, ]+normals['s', ]*lead, b[1L, ])
    edge2 <- rbind(a[nrow(a), ], a[nrow(a), ]+normals['q', ]*lead, edge2,
      b[nrow(b), ]+normals['s', ]*lead, b[nrow(b), ])
    path <- rbind(a, edge2, b[nrow(b):1L, , drop = FALSE], edge1[nrow(edge1):1L, , drop = FALSE])
    if (any(path[, 1] < 0 | path[, 1] > canvas[1] | path[, 2] < 0 | path[, 2] > canvas[2])) next
    if (any(ggchord_route_polygon_hits(visible(path), obstacles))) next
    attr(path, 'endpoint_normals') <- attr(direct, 'endpoint_normals')
    return(list(path = path, status = if (route == 'bundle') 'bundled_clear' else 'routed_clear',
      collisions = 0L, corridor = line))
  }
  ggchord_stop('Fusion: no collision-free ribbon corridor preserves the complete interval fronts; move/enlarge panels, reduce the interval, or use route = "direct" and inspect diagnostics')
}

# Device metrics of final native layer text. Rectangles conservatively enclose
# rotated glyphs; arc strings receive a larger envelope, not an exact-glyph claim.
ggchord_fusion_text_obstacles <- function(panel, id) {
  bounds <- panel$context$panel_bounds_mm
  result <- list()
  for (index in seq_along(panel$built$data)) {
    layer <- panel$built$data[[index]]
    if (!all(c('x', 'y', 'label', 'size') %in% names(layer))) next
    if ('.component' %in% names(layer)) layer <- layer[layer$.component %in%
      c('text', 'arc_text', 'track_axis_text', 'track_axis_title', 'label'), , drop = FALSE]
    layer <- layer[!is.na(layer$label) & nzchar(as.character(layer$label)), , drop = FALSE]
    if (!nrow(layer)) next
    xy <- panel$built$layout$coord$transform(layer, panel$built$layout$panel_params[[1L]])
    # Use the same installed ggplot2 resolver as GeomText for inward/outward
    # justification; do not approximate it differently in composition space.
    justify <- utils::getFromNamespace('compute_just', 'ggplot2')
    if ('hjust' %in% names(xy)) layer$hjust <- justify(xy$hjust, xy$x, xy$y, xy$angle %||% 0)
    if ('vjust' %in% names(xy)) layer$vjust <- justify(xy$vjust, xy$y, xy$x, xy$angle %||% 0)
    for (i in seq_len(nrow(layer))) {
      d <- layer[i, , drop = FALSE]
      if ('alpha' %in% names(d) && !is.na(d$alpha) && d$alpha == 0) next
      if ('colour' %in% names(d) && (is.na(d$colour) ||
          grDevices::col2rgb(d$colour, alpha = TRUE)[4L] == 0)) next
      scalar <- function(name, default) if (name %in% names(d)) d[[name]][[1L]] else default
      geom <- panel$built$plot$layers[[index]]
      points_per_unit <- switch(geom$geom_params$size.unit %||% 'mm',
        mm = 72.27/25.4, cm = 72.27/2.54, 'in' = 72.27, pt = 1, pc = 12, 72.27/25.4)
      label <- if (isTRUE(geom$geom_params$parse)) parse(text = as.character(d$label)) else d$label
      glyph <- grid::textGrob(label, gp = grid::gpar(
        fontsize = as.numeric(d$size)*points_per_unit, fontfamily = scalar('family', ''),
        fontface = scalar('fontface', 'plain'), lineheight = scalar('lineheight', 1.2)))
      w <- grid::convertWidth(grid::grobWidth(glyph), 'mm', TRUE)
      h <- grid::convertHeight(grid::grobHeight(glyph), 'mm', TRUE)
      if (inherits(geom$geom, 'GeomResidue')) {
        corners <- d[rep(1L, 2L), , drop = FALSE]
        corners$x <- d$x + c(-.46, .46); corners$y <- as.numeric(d$y) + c(-.4, .4)
        corners <- panel$built$layout$coord$transform(corners, panel$built$layout$panel_params[[1L]])
        if (abs(diff(corners$x))*diff(bounds[1:2]) < w*1.15 ||
            abs(diff(corners$y))*diff(bounds[3:4]) < h*1.15) next
      }
      if (inherits(geom$geom, 'GeomLabel')) {
        padding <- geom$geom_params$label.padding %||% grid::unit(.25, 'lines')
        grid::pushViewport(grid::viewport(gp = grid::gpar(fontsize = as.numeric(d$size)*points_per_unit)))
        pad <- grid::convertUnit(padding, 'mm', valueOnly = TRUE)
        grid::popViewport()
        w <- w + 2*pad; h <- h + 2*pad
      }
      angle <- scalar('angle', 0)*pi/180
      x <- c(-scalar('hjust', .5), 1-scalar('hjust', .5))*w
      y <- c(-scalar('vjust', .5), 1-scalar('vjust', .5))*h
      corners <- as.matrix(expand.grid(x = x, y = y)) %*%
        matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2, byrow = TRUE)
      arc <- isTRUE(scalar('.draw_as_arc', FALSE)) || identical(scalar('.component', ''), 'arc_text')
      if (arc) corners <- rbind(c(-w-h, -w-h), c(w+h, w+h))
      center <- c(bounds[1]+xy$x[i]*diff(bounds[1:2]), bounds[3]+xy$y[i]*diff(bounds[3:4]))
      if (any(!is.finite(c(center, corners)))) next
      result[[length(result)+1L]] <- data.frame(panel_id = id, kind = 'text',
        xmin = center[1]+min(corners[, 1]), xmax = center[1]+max(corners[, 1]),
        ymin = center[2]+min(corners[, 2]), ymax = center[2]+max(corners[, 2]),
        measurement = if (arc) 'conservative_arc_envelope' else 'device_text_envelope')
    }
  }
  ggchord_rbind_fill(result) %||% data.frame()
}

# Resolve guide-cell bounds through the same gtable layout as the final panel.
ggchord_fusion_guide_obstacles <- function(table, id, origin) {
  at <- which(grepl('^(guide-box|title$|subtitle$|caption$|tag$|xlab|ylab|axis-)', table$layout$name) &
    !vapply(table$grobs, inherits, logical(1), 'zeroGrob'))
  if (!length(at)) return(data.frame())
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    length(table$heights), length(table$widths), widths = table$widths,
    heights = table$heights, respect = table$respect)))
  on.exit(grid::popViewport())
  ggchord_rbind_fill(lapply(at, function(i) {
    entry <- table$layout[i, ]
    grid::pushViewport(grid::viewport(layout.pos.row = seq.int(entry$t, entry$b),
      layout.pos.col = seq.int(entry$l, entry$r)))
    on.exit(grid::popViewport())
    lo <- grid::deviceLoc(grid::unit(0, 'npc'), grid::unit(0, 'npc'), valueOnly = TRUE)
    hi <- grid::deviceLoc(grid::unit(1, 'npc'), grid::unit(1, 'npc'), valueOnly = TRUE)
    data.frame(panel_id = id, kind = if (grepl('^guide-box', entry$name)) 'guide' else 'decoration', xmin = lo$x*25.4-origin[1],
      xmax = hi$x*25.4-origin[1], ymin = lo$y*25.4-origin[2], ymax = hi$y*25.4-origin[2],
      measurement = 'resolved_gtable_cell')
  }))
}
