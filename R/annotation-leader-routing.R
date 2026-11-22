# Protect static inner/axis text after the shared exterior queue has settled.
# Positions and real feature anchors are unchanged; only collided leader pieces
# can acquire a short deterministic polyline detour.
ggchord_annotation_text_obstacles <- function(registry, layout, include_axis = TRUE,
    scale = c(1, 1), endpoints_only = FALSE, include_outer = TRUE) {
  all <- list()
  unit <- if (identical(scale, c(1, 1))) layout$text_units_per_inch %||% .35 else 1
  add <- function(text, owner) {
    if (!nrow(text)) return()
    keep <- !is.na(text$label) & nzchar(text$label)
    if (!include_outer && "annotation_side" %in% names(text))
      keep <- keep & !text$annotation_side %in% "outer"
    text <- text[keep, , drop = FALSE]
    if (!nrow(text)) return()
    text$text <- text$label
    text$text_x <- text$x / scale[1L]; text$text_y <- text$y / scale[2L]
    text$text_angle <- text$angle %||% 0
    arc <- is.finite(text$.arc_units %||% rep(NA_real_, nrow(text)))
    arc_text <- ggchord_arc_text_layout(text[arc, , drop = FALSE],
      units_per_inch = unit, boxes_only = TRUE)
    boxes <- ggchord_rbind_fill(list(
      ggchord_text_boxes(text[!arc, , drop = FALSE], units_per_inch = unit, box_padding = 0),
      arc_text$boxes))
    rows <- c(which(!arc), which(arc)[arc_text$owners])
    boxes$owner <- owner
    boxes$source_row <- (text$source_row %||% rep(NA_integer_, nrow(text)))[rows]
    object <- if ('.arc_parent_group' %in% names(text)) as.character(text$.arc_parent_group) else
      if ('source_row' %in% names(text)) as.character(text$source_row) else as.character(seq_len(nrow(text)))
    boxes$text_object <- paste(owner, object[rows], sep = '/')
    all[[length(all)+1L]] <<- boxes
  }
  for (owner in names(registry)) {
    d <- registry[[owner]]$gene_label_repel
    if (is.data.frame(d) && nrow(d)) {
      keep <- d$.component %in% c('text', 'arc_text')
      if ('.draw_as_arc' %in% names(d)) keep <- keep &
        !(d$.component %in% 'text' & d$.draw_as_arc %in% TRUE)
      add(d[keep, , drop = FALSE], owner)
    }
    d <- registry[[owner]]$restriction_site
    if (is.data.frame(d) && nrow(d)) add(d[d$restriction_component %in% 'label', , drop = FALSE], owner)
    d <- registry[[owner]]$seq_center_label
    if (is.data.frame(d) && nrow(d)) add(d, owner)
    d <- registry[[owner]]$seq_label
    if (is.data.frame(d) && nrow(d)) add(d, owner)
    d <- registry[[owner]]$domain_label
    if (is.data.frame(d) && nrow(d))
      add(d[d$.component %in% "text", , drop = FALSE], owner)
  }
  if (include_axis) {
    axis_layout <- layout
    axis_layout$layer_geometry <- registry
    axis <- ggchord_axis_geometry(axis_layout, clear = FALSE)
    keep <- axis$.component %in% 'text'
    if (endpoints_only) keep <- keep & ggchord_axis_endpoints(axis)
    if (any(keep)) add(axis[keep, , drop = FALSE], '.axis')
  }
  ggchord_rbind_fill(all) %||% data.frame()
}

ggchord_annotation_metric <- function(layout) {
  metric <- ggchord_render_context()$units_per_mm
  if (length(metric) != 2L || any(!is.finite(metric) | metric <= 0))
    metric <- rep((layout$text_units_per_inch %||% .35) / 25.4, 2L)
  setNames(metric, c("x", "y"))
}

ggchord_annotation_segments <- function(data, scale = c(1, 1)) {
  if (!NROW(data)) return(data.frame())
  valid <- is.finite(data$x) & is.finite(data$y) &
    is.finite(data$xend) & is.finite(data$yend)
  data.frame(x0 = data$x[valid] / scale[1L], y0 = data$y[valid] / scale[2L],
    x1 = data$xend[valid] / scale[1L], y1 = data$yend[valid] / scale[2L])
}

ggchord_annotation_backbones <- function(layout, scale = c(1, 1)) {
  lapply(layout$seq_arcs, function(path) {
    if (NROW(path) < 2L) return(data.frame())
    data.frame(x0 = head(path$x, -1L) / scale[1L], y0 = head(path$y, -1L) / scale[2L],
      x1 = tail(path$x, -1L) / scale[1L], y1 = tail(path$y, -1L) / scale[2L])
  })
}

ggchord_annotation_polygons <- function(registry, scale) {
  polygons <- list()
  for (owner in names(registry)) {
    d <- registry[[owner]]$gene_poly
    if (!NROW(d)) next
    d <- d[d$.component %in% "polygon", , drop = FALSE]
    if (!nrow(d)) next
    d$x <- d$x / scale[1L]; d$y <- d$y / scale[2L]
    d$.label_polygon_owner <- owner
    polygons <- c(polygons, unname(split(d, d$group)))
  }
  polygons
}

# Routing can break a leader into several segments. A crossing belongs to its
# source label, so preserve one line type through every bend.
ggchord_mark_leader_crossings <- function(registry, layout) {
  paths <- ggchord_annotation_backbones(layout)
  text_boxes <- NULL
  for (id in names(registry)) {
    data <- registry[[id]]$domain_label
    if (!NROW(data)) next
    leaders <- which(data$.component %in% "leader")
    if (!length(leaders)) next
    starts <- split(leaders, data$source_row[leaders])
    for (at in starts) {
      path <- paths[[data$accver[at[1L]]]]
      line <- ggchord_annotation_segments(data[at, , drop = FALSE])
      if (NROW(path))
        data$.leader_crosses_backbone[at] <- ggchord_inner_leaders_cross(line, path)
      if (isTRUE(layout$circular) && data$annotation_side[at[1L]] %in% "inner") {
        if (is.null(text_boxes)) text_boxes <- ggchord_annotation_text_obstacles(
          registry, layout, include_axis = FALSE, include_outer = FALSE)
        other <- !(text_boxes$owner %in% id & text_boxes$source_row %in% data$source_row[at[1L]])
        data$.leader_crosses_text[at] <- ggchord_inner_leader_hits(line,
          text_boxes[other, , drop = FALSE])
      }
    }
    registry[[id]]$domain_label <- data
  }
  registry
}

# The exterior queue already provides direct or one-shoulder routes.
# Preserve those routes instead of turning collisions into multi-bend graph
# detours. Text positions and fixed source roots remain coordinate-owned.
ggchord_detour_annotation_leaders <- function(registry, layout) {
  if (!isTRUE(layout$circular)) return(registry)
  ggchord_measure_external_leader_crossings(registry)
}
