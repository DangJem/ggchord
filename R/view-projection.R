# Final view projection changes geometry, never source input tables.
ggchord_project_vertical_data <- function(data) {
  if (!is.data.frame(data) || !nrow(data)) return(data)
  for (pair in list(c("x", "y"), c("x0", "y0"), c("x1", "y1"),
      c("text_x", "text_y"), c("label_x", "label_y"), c("anchor_x", "anchor_y"))) {
    if (!all(pair %in% names(data))) next
    x <- data[[pair[1L]]]; y <- data[[pair[2L]]]
    data[[pair[1L]]] <- -y; data[[pair[2L]]] <- x
  }
  for (field in intersect(c("angle", "text_angle", "label_angle"), names(data))) {
    value <- (data[[field]] + 90 + 180) %% 360 - 180
    # Preserve readable text rather than turning labels upside down.
    value[value > 90] <- value[value > 90] - 180
    value[value < -90] <- value[value < -90] + 180
    data[[field]] <- value
  }
  data
}

ggchord_project_vertical_layout <- function(layout) {
  if (isTRUE(layout$vertical_projected)) return(layout)
  geometry <- c("gene_polys", "gene_labels", "gene_label_segments", "region_polys",
    "ribbon_polys", "link_lines", "seq_labels_df", "axis_lines", "axis_ticks",
    "sequence_scales", "track_axes", "seq_breaks")
  for (key in geometry) layout[[key]] <- ggchord_project_vertical_data(layout[[key]])
  layout$seq_arcs <- lapply(layout$seq_arcs, ggchord_project_vertical_data)
  layout$layer_geometry <- lapply(layout$layer_geometry, function(owner)
    lapply(owner, ggchord_project_vertical_data))
  # The inverse mapper retains a sequence-local frame; expose final bounds
  # separately, together with the original along-axis range.
  if (is.data.frame(layout$sequence_layout)) {
    d <- layout$sequence_layout
    d$axis_min <- d$xmin; d$axis_max <- d$xmax
    d$ymin <- d$xmin; d$ymax <- d$xmax
    d$xmin <- d$xmax <- -d$y
    d$y <- (d$ymin + d$ymax) / 2
    layout$sequence_layout <- d
  }
  ext <- layout$extremes
  if (!is.null(ext)) layout$extremes <- list(x_min = -ext$y_max, x_max = -ext$y_min,
    y_min = ext$x_min, y_max = ext$x_max)
  layout$seq_axis <- "vertical"
  layout$vertical_projected <- TRUE
  layout
}
