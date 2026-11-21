# One source backbone frame, before native scales/Coord and device conversion.
# Mirrored and radial views reuse their own explicit instance transformations.
ggchord_panel_backbone_frame <- function(layout, id, position) {
  ref <- layout$sequence_reference
  if (!isTRUE(layout$collinear)) {
    frame <- ggchord_track_frame(position, id, layout)
    return(list(xy = frame$base, normal = frame$normal))
  }
  xy <- cbind(x = ggchord_collinear_x(position, id, ref),
    y = rep(unname(ref$y[id]), length(position)))
  normal <- cbind(x = rep(0, length(position)), y = rep(1, length(position)))
  project <- if (isTRUE(layout$sankey_projected)) function(p)
    ggchord_sankey_points(p[, 1], p[, 2], id, layout$sankey_frames) else
    if (isTRUE(layout$radial_projected)) function(p)
      ggchord_radial_points(p[, 1], p[, 2], id, layout$radial_frames) else
    if (identical(layout$seq_axis, "vertical")) function(p) cbind(x = p[, 2], y = p[, 1]) else identity
  tip <- project(xy + normal); xy <- project(xy)
  list(xy = xy, normal = tip - xy)
}

# Move child exports from their local composition to the enclosing canvas.
# Only known geometric columns are shifted, never arbitrary source metadata.
ggchord_fusion_translate_table <- function(data, offset, x = character(), y = character()) {
  if (!is.data.frame(data) || !nrow(data)) return(data)
  for (name in intersect(x, names(data))) data[[name]] <- data[[name]] + offset[1L]
  for (name in intersect(y, names(data))) data[[name]] <- data[[name]] + offset[2L]
  data
}

ggchord_fusion_prefix_table <- function(data, parent, columns) {
  if (!is.data.frame(data) || !nrow(data)) return(data)
  for (name in intersect(columns, names(data))) data[[name]] <- paste0(parent, "/", data[[name]])
  data
}

ggchord_fusion_nested_export <- function(value, parent, offset) {
  prefix <- function(d, cols) if (is.null(parent)) d else ggchord_fusion_prefix_table(d, parent, cols)
  value$panels <- prefix(ggchord_fusion_translate_table(value$panels, offset,
    c("container_xmin", "container_xmax", "panel_xmin", "panel_xmax"),
    c("container_ymin", "container_ymax", "panel_ymin", "panel_ymax")), "panel_id")
  if (length(value$metadata$content_bounds_mm) == 4L)
    value$metadata$content_bounds_mm <- value$metadata$content_bounds_mm + rep(offset, each = 2L)
  if (!is.null(value$metadata$plot_decoration)) value$metadata$plot_decoration$bounds_mm <-
    ggchord_fusion_translate_table(value$metadata$plot_decoration$bounds_mm,
      offset, c("xmin", "xmax"), c("ymin", "ymax"))
  value$metadata$guide_bounds_mm <- ggchord_fusion_translate_table(value$metadata$guide_bounds_mm,
    offset, c("xmin", "xmax"), c("ymin", "ymax"))
  value$metadata$guide_layout$bounds <- ggchord_fusion_translate_table(value$metadata$guide_layout$bounds,
    offset, c("xmin_mm", "xmax_mm"), c("ymin_mm", "ymax_mm"))
  value$instances <- prefix(value$instances, "panel_id")
  value$routing_obstacles <- prefix(ggchord_fusion_translate_table(value$routing_obstacles, offset,
    c("xmin", "xmax"), c("ymin", "ymax")), "panel_id")
  prepend <- function(data, fields) {
    if (is.null(parent) || !is.data.frame(data) || !nrow(data)) return(data)
    for (name in intersect(fields, names(data))) data[[name]] <- lapply(data[[name]], function(path)
      c(ggchord_fusion_path_ids(parent), path))
    data
  }
  scope <- function(data) {
    if (is.null(parent) || !is.data.frame(data) || !nrow(data) || !"composition_id" %in% names(data)) return(data)
    data$composition_id <- ifelse(data$composition_id == "", parent,
      paste0(parent, "/", data$composition_id))
    data
  }
  value$connections <- scope(prepend(prefix(value$connections, c("qpanel", "spanel")),
    c("qpanel_path", "spanel_path")))
  for (name in c("fronts", "paths", "zoom_paths")) {
    value[[name]] <- ggchord_fusion_translate_table(value[[name]], offset, c("x", "xend"), c("y", "yend"))
    value[[name]] <- scope(prefix(value[[name]], "panel_id"))
  }
  value$zoom_regions <- scope(prepend(value$zoom_regions, c("from_panel_path", "to_panel_path")))
  value$overlaps <- prefix(ggchord_fusion_translate_table(value$overlaps, offset,
    c("xmin", "xmax"), c("ymin", "ymax")), c("under", "over"))
  value
}

# Clip to a leaf panel, then restore the full composition coordinate frame.
# Native grid rectangles work without rasterising plots or device-specific masks.
ggchord_fusion_clip_connections <- function(grobs, bounds, canvas) {
  full <- viewport(x = -bounds[1], y = -bounds[3], width = canvas[1], height = canvas[2],
    just = c("left", "bottom"), default.units = "mm", clip = "inherit")
  clip <- viewport(x = bounds[1], y = bounds[3], width = diff(bounds[1:2]), height = diff(bounds[3:4]),
    just = c("left", "bottom"), default.units = "mm", clip = "on")
  grobTree(do.call(grobTree, c(unname(grobs), list(vp = full))), vp = clip)
}

# Disjoint rectangles avoid repainting translucent connections more than once.
ggchord_fusion_rect_subtract <- function(rect, exclusions) {
  pieces <- list(unname(rect))
  for (b in exclusions) {
    next_pieces <- list()
    for (a in pieces) {
      lo <- c(max(a[1], b[1]), max(a[3], b[3])); hi <- c(min(a[2], b[2]), min(a[4], b[4]))
      if (any(lo >= hi)) { next_pieces[[length(next_pieces) + 1L]] <- a; next }
      candidates <- list(c(a[1], lo[1], a[3], a[4]), c(hi[1], a[2], a[3], a[4]),
        c(lo[1], hi[1], a[3], lo[2]), c(lo[1], hi[1], hi[2], a[4]))
      next_pieces <- c(next_pieces, Filter(function(r) r[1] < r[2] && r[3] < r[4], candidates))
    }
    pieces <- next_pieces
  }
  pieces
}
