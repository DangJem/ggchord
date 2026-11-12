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
    if (identical(layout$seq_axis, "vertical")) function(p) cbind(x = -p[, 2], y = p[, 1]) else identity
  tip <- project(xy + normal); xy <- project(xy)
  list(xy = xy, normal = tip - xy)
}

# Explicit physical or relative rectangles are visual allocations only.
ggchord_fusion_box <- function(box) {
  if (is.null(box)) return(invisible(NULL))
  physical <- grid::is.unit(box)
  if (length(box) != 4L || (!physical && !is.numeric(box)))
    ggchord_stop("Fusion: box must be xmin, xmax, ymin, ymax (four fractions or absolute units)")
  if (!is.null(names(box))) {
    if (anyDuplicated(names(box)) || !setequal(names(box), c("xmin", "xmax", "ymin", "ymax")))
      ggchord_stop("Fusion: name box coordinates xmin, xmax, ymin, ymax")
    box <- box[c("xmin", "xmax", "ymin", "ymax")]
  }
  value <- if (physical) vapply(seq_along(box), function(i) ggchord_unit_inches(box[i]) * 25.4, numeric(1)) else as.numeric(box)
  if (any(!is.finite(value)) || min(value) < 0 || value[1] >= value[2] || value[3] >= value[4] ||
      (!physical && max(value) > 1)) ggchord_stop("Fusion: box must have ordered nonnegative bounds within its allocation frame")
  invisible(box)
}

ggchord_fusion_sides <- function() c("right", "top", "left", "bottom", "topright", "topleft", "bottomright", "bottomleft")

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
  prefix <- function(d, cols) ggchord_fusion_prefix_table(d, parent, cols)
  value$panels <- prefix(ggchord_fusion_translate_table(value$panels, offset,
    c("container_xmin", "container_xmax", "panel_xmin", "panel_xmax"),
    c("container_ymin", "container_ymax", "panel_ymin", "panel_ymax")), "panel_id")
  value$instances <- prefix(value$instances, "panel_id")
  value$connections <- prefix(value$connections, c("qpanel", "spanel", "layer_id"))
  for (name in c("fronts", "paths", "focus_fronts", "leaders")) {
    value[[name]] <- ggchord_fusion_translate_table(value[[name]], offset, c("x", "xend"), c("y", "yend"))
    value[[name]] <- prefix(value[[name]], c("panel_id", "layer_id", "composition_id"))
  }
  value$focus_regions <- prefix(value$focus_regions, c("source_panel", "detail_panel", "composition_id"))
  value$overlaps <- prefix(ggchord_fusion_translate_table(value$overlaps, offset,
    c("xmin", "xmax"), c("ymin", "ymax")), c("under", "over"))
  value
}

ggchord_fusion_panel_paths <- function(plot, parent) {
  if (!inherits(plot, "ggchord_fusion")) return(parent)
  c(parent, unlist(lapply(names(plot$panels), function(id)
    ggchord_fusion_panel_paths(plot$panels[[id]]$plot, paste0(parent, "/", id))), use.names = FALSE))
}

# Clip to a leaf panel, then restore the full composition coordinate frame.
# Native grid rectangles work without rasterising plots or device-specific masks.
ggchord_fusion_clip_connections <- function(grobs, bounds, canvas) {
  full <- grid::viewport(x = -bounds[1], y = -bounds[3], width = canvas[1], height = canvas[2],
    just = c("left", "bottom"), default.units = "mm", clip = "inherit")
  clip <- grid::viewport(x = bounds[1], y = bounds[3], width = diff(bounds[1:2]), height = diff(bounds[3:4]),
    just = c("left", "bottom"), default.units = "mm", clip = "on")
  grid::grobTree(do.call(grid::grobTree, c(unname(grobs), list(vp = full))), vp = clip)
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

ggchord_fusion_box_rect <- function(box, canvas) {
  box <- ggchord_fusion_box(box)
  if (grid::is.unit(box)) return(vapply(seq_along(box), function(i)
    ggchord_unit_inches(box[i]) * 25.4, numeric(1)) / rep(canvas, each = 2L))
  unname(box)
}
