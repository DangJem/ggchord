# Bend along a sequence-local frame before it is rotated or translated.
# A transverse displacement keeps the source coordinate and the endpoint
# positions fixed. Applying it to every offset preserves the track and
# feature clearances used by the connection routing code.
ggchord_backbone_bend <- function(x, id, ref, slope = FALSE) {
  span <- unname(ref$span[id])
  strength <- unname(ref$curvature[id])
  t <- (x - unname(ref$x_offset[id])) / span
  value <- if (slope) .4 * strength * (1 - 2 * t) else
    .4 * strength * span * t * (1 - t)
  value[is.finite(t) & (t < 0 | t > 1)] <- 0
  if (any(is.infinite(value))) ggchord_stop("Curvature produces non-finite coordinates")
  value
}

# Subdivide long edges in the original frame so polygons and paths follow
# the bend between their stored vertices. Keep metadata with its source row.
ggchord_bend_subdivide <- function(d, ref, closed = FALSE) {
  if (nrow(d) < 2L || !all(is.finite(d$x) & is.finite(d$y))) return(d)
  id <- as.character(d$accver[1L])
  if (!isTRUE(unname(ref$curvature[id]) != 0)) return(d)
  n <- nrow(d)
  ends <- if (closed) c(seq.int(2L, n), 1L) else seq.int(2L, n)
  starts <- seq_along(ends)
  pieces <- pmax(1, pmin(64, ceiling(abs(d$x[ends] - d$x[starts]) /
    ref$span[[id]] * 32 * sqrt(max(1, abs(ref$curvature[[id]]))))))
  edge <- rep(starts, pieces)
  fraction <- sequence(pieces, from = 0L) / pieces[edge]
  next_ <- ends[edge]
  out <- d[edge, , drop = FALSE]
  out$x <- (1 - fraction) * d$x[edge] + fraction * d$x[next_]
  out$y <- (1 - fraction) * d$y[edge] + fraction * d$y[next_]
  if (!closed) out <- rbind(out, d[n, , drop = FALSE])
  out
}

ggchord_bend_local_data <- function(d, ref, component) {
  if (!is.data.frame(d) || !nrow(d) || !"accver" %in% names(d)) return(d)
  path_components <- c("reference", "axis_lines", "seq", "gene_poly", "seq_region", "track")
  if (component %in% path_components && all(c("x", "y") %in% names(d))) {
    kind <- d$.component %||% rep("", nrow(d))
    strip <- d$group %||% seq_len(nrow(d))
    if (component %in% c("reference", "axis_lines")) strip <- d$accver
    key <- interaction(d$accver, strip, kind, drop = TRUE)
    by_path <- split(seq_len(nrow(d)), match(key, unique(key)))
    d <- do.call(rbind, lapply(by_path, function(rows) {
      shape <- kind[rows[1L]]
      if (component == "track" && !shape %in% c("track_line", "track_poly", "track_area_fill"))
        return(d[rows, , drop = FALSE])
      closed <- component == "seq_region" ||
        shape %in% c("polygon", "band", "track_poly", "track_area_fill")
      ggchord_bend_subdivide(d[rows, , drop = FALSE], ref, closed)
    }))
  }
  for (pair in list(c("x", "y"), c("x0", "y0"), c("x1", "y1"), c("xend", "yend"),
      c("text_x", "text_y"), c("label_x", "label_y"), c("anchor_x", "anchor_y"))) {
    if (!all(pair %in% names(d))) next
    at <- which(d$accver %in% names(ref$curvature) & is.finite(d[[pair[1L]]]))
    d[[pair[2L]]][at] <- d[[pair[2L]]][at] +
      ggchord_backbone_bend(d[[pair[1L]]][at], as.character(d$accver[at]), ref)
  }
  d
}

ggchord_bend_front <- function(id, x, y, ref) {
  if (!isTRUE(unname(ref$curvature[id]) != 0)) return(cbind(x = x, y = y))
  d <- ggchord_bend_subdivide(data.frame(x = x, y = y, accver = id), ref)
  cbind(x = d$x, y = d$y + ggchord_backbone_bend(d$x, id, ref))
}

ggchord_bend_local_layout <- function(layout) {
  ref <- layout$sequence_reference
  for (owner in names(layout$layer_geometry)) {
    for (component in setdiff(names(layout$layer_geometry[[owner]]), c("ribbon", "link", "seq_scale"))) {
      d <- layout$layer_geometry[[owner]][[component]]
      layout$layer_geometry[[owner]][[component]] <- ggchord_bend_local_data(d, ref, component)
    }
  }
  layout$seq_arcs <- lapply(layout$seq_arcs, ggchord_bend_local_data, ref = ref, component = "reference")
  for (field in c("axis_lines", "axis_ticks", "gene_labels", "gene_label_segments",
      "seq_labels_df", "track_axes", "seq_breaks", "restriction_sites")) {
    layout[[field]] <- ggchord_bend_local_data(layout[[field]], ref, field)
  }
  ticks <- layout$axis_ticks
  if (NROW(ticks)) {
    orientation <- as.character(layout$axisLabelOrientation[as.character(ticks$accver)])
    at <- which(orientation %in% c("auto", "parallel", "perpendicular", "tangent"))
    slope <- ggchord_backbone_bend(ticks$label_x[at], as.character(ticks$accver[at]), ref, slope = TRUE)
    ticks$label_angle[at] <- (ticks$label_angle[at] + atan(slope) * 180 / pi + 90) %% 180 - 90
    layout$axis_ticks <- ticks
  }
  for (pair in list(c("gene_polys", "gene_poly"), c("region_polys", "seq_region"))) {
    layout[[pair[1L]]] <- ggchord_rbind_fill(lapply(layout$layer_geometry, `[[`, pair[2L])) %||% data.frame()
  }
  layout
}
