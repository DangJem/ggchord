# Protect static inner/axis text after the shared exterior queue has settled.
# Positions and real feature anchors are unchanged; only collided leader pieces
# can acquire a short deterministic polyline detour.
ggchord_annotation_text_obstacles <- function(registry, layout, include_axis = TRUE) {
  all <- list()
  unit <- layout$text_units_per_inch %||% .35
  add <- function(text, owner) {
    if (!nrow(text)) return()
    text <- text[!is.na(text$label) & nzchar(text$label), , drop = FALSE]
    if (!nrow(text)) return()
    text$text <- text$label; text$text_x <- text$x; text$text_y <- text$y
    text$text_angle <- text$angle %||% 0
    boxes <- ggchord_text_boxes(text, units_per_inch = unit, box_padding = 0)
    boxes$owner <- owner; boxes$source_row <- text$source_row %||% rep(NA_integer_, nrow(text))
    object <- if ('.arc_parent_group' %in% names(text)) as.character(text$.arc_parent_group) else
      if ('source_row' %in% names(text)) as.character(text$source_row) else as.character(seq_len(nrow(text)))
    boxes$text_object <- paste(owner, object, sep = '/')
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
  }
  if (include_axis) {
    axis_layout <- layout
    axis_layout$layer_geometry <- registry
    axis <- ggchord_axis_geometry(axis_layout)
    if (nrow(axis)) add(axis[axis$.component %in% 'text', , drop = FALSE], '.axis')
  }
  ggchord_rbind_fill(all) %||% data.frame()
}

ggchord_detour_annotation_leaders <- function(registry, layout) {
  if (!isTRUE(layout$circular)) return(registry)
  needed <- any(vapply(registry, function(owner) {
    sites <- owner$restriction_site
    if (is.data.frame(sites) && nrow(sites) && any(sites$restriction_component %in% 'leader')) return(TRUE)
    d <- owner$gene_label_repel
    if (!is.data.frame(d) || !nrow(d)) return(FALSE)
    any(d$.component %in% 'segment')
  }, logical(1)))
  if (!needed) return(registry)
  boxes <- ggchord_annotation_text_obstacles(registry, layout)
  if (!nrow(boxes)) return(registry)
  unit <- layout$text_units_per_inch %||% .35
  clearance <- .25 / 25.4 * unit
  # Work in translated native layout units. The graph solver is independent of
  # Fusion panels and accepts a rectangle enclosing these existing annotations.
  anchors <- unlist(lapply(registry, function(owner) {
    d <- owner$gene_label_repel
    if (!is.data.frame(d)) return(NULL)
    unlist(d[intersect(c('x', 'y', 'xend', 'yend'), names(d))], use.names = FALSE)
  }), use.names = FALSE)
  extent <- range(c(boxes$xmin, boxes$xmax, boxes$ymin, boxes$ymax, anchors), finite = TRUE)
  margin <- max(clearance*4, diff(extent)*.05)
  origin <- rep(extent[1L]-margin, 2L); canvas <- rep(diff(extent)+2*margin, 2L)
  search_boxes <- data.frame(xmin = boxes$xmin-origin[1L], xmax = boxes$xmax-origin[1L],
    ymin = boxes$ymin-origin[2L], ymax = boxes$ymax-origin[2L])
  for (owner in names(registry)) {
    geometry <- registry[[owner]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry)) next
    # Feature-track labels retain their direct / single-shoulder route.
    # Crossing a feature is allowed; graph detours obscure source association.
    if ("label_layout" %in% names(geometry) &&
        any(geometry$label_layout %in% "feature")) next
    external <- unique(geometry$source_row[geometry$.component %in% 'text' & ggchord_outer_callouts(geometry)])
    for (source in external[!is.na(external)]) {
      row <- which(geometry$.component %in% 'text' & geometry$source_row == source)[1L]
      own <- boxes$owner == owner & !is.na(boxes$source_row) & boxes$source_row == source
      obstacles <- search_boxes[!own, , drop = FALSE]
      oriented <- boxes[!own, , drop = FALSE]
      obstacle_group <- boxes$text_object[!own]
      at <- which(geometry$.component %in% 'segment' & geometry$source_row == source)
      active_rows <- at[geometry$x[at] != geometry$xend[at] | geometry$y[at] != geometry$yend[at]]
      hit <- unique(unlist(lapply(active_rows, function(i) {
        possible <- which(ggchord_route_segment_hits(c(geometry$x[i], geometry$y[i])-origin,
          c(geometry$xend[i], geometry$yend[i])-origin, obstacles))
        segment <- data.frame(x0 = geometry$x[i], y0 = geometry$y[i], x1 = geometry$xend[i], y1 = geometry$yend[i])
        possible[vapply(possible, function(j) ggchord_inner_leader_hits(segment,
          oriented[j, , drop = FALSE]), logical(1))]
      }), use.names = FALSE))
      changed <- FALSE; unresolved <- FALSE
      if (length(hit)) {
        # Keep the true root and final text-edge attachment, but do not force a
        # detour back through an obsolete intermediate shoulder.
        distance <- function(x, y) (x-geometry$text_x[row])^2 + (y-geometry$text_y[row])^2
        root <- active_rows[which.max(distance(geometry$x[active_rows], geometry$y[active_rows]))]
        terminal <- active_rows[which.min(distance(geometry$xend[active_rows], geometry$yend[active_rows]))]
        start <- c(geometry$x[root], geometry$y[root])-origin
        end <- c(geometry$xend[terminal], geometry$yend[terminal])-origin
        active <- which(obstacle_group %in% obstacle_group[hit]); path <- NULL
        # Grow only the blockers actually touched by candidate detours. Every
        # accepted result is rechecked against all text, not just this subset.
        for (attempt in seq_len(length(unique(obstacle_group))+1L)) {
          path <- ggchord_route_shortest(start, end, obstacles[active, , drop = FALSE], canvas,
            clearance = clearance)
          if (is.null(path)) break
          touched <- unique(unlist(lapply(seq_len(nrow(path)-1L), function(k)
            which(ggchord_route_segment_hits(path[k, ], path[k+1L, ], obstacles))), use.names = FALSE))
          if (!length(touched)) break
          newer <- which(obstacle_group %in% obstacle_group[union(active, touched)])
          path <- NULL
          if (length(newer) == length(active)) break
          active <- newer
        }
        if (is.null(path)) unresolved <- TRUE else {
          path <- sweep(path, 2L, origin, '+')
          d <- geometry[rep(root, nrow(path)-1L), , drop = FALSE]
          d$x <- head(path[,1], -1L); d$y <- head(path[,2], -1L)
          d$xend <- tail(path[,1], -1L); d$yend <- tail(path[,2], -1L)
          for (pair in list(c('x0','x'), c('y0','y'), c('x1','xend'), c('y1','yend')))
            if (pair[1L] %in% names(d)) d[[pair[1L]]] <- d[[pair[2L]]]
          geometry <- ggchord_rbind_fill(list(geometry[-at, , drop = FALSE], d))
          changed <- TRUE
        }
      }
      rows <- !is.na(geometry$source_row) & geometry$source_row == source & geometry$.component %in% c('text', 'segment')
      if (!'leader_text_status' %in% names(geometry)) geometry$leader_text_status <- NA_character_
      geometry$leader_text_status[rows] <- if (unresolved) 'static_text_unresolved' else if (changed) 'static_text_detour' else 'static_text_clear'
    }
    # The old replacement pass leaves unused pieces with coincident endpoints.
    # They are placeholders, not leaders or marker observations; removing them
    # also prevents round caps/arrow settings from painting accidental dots.
    unused <- geometry$.component %in% 'segment' & is.finite(geometry$x) & is.finite(geometry$y) &
      is.finite(geometry$xend) & is.finite(geometry$yend) &
      geometry$x == geometry$xend & geometry$y == geometry$yend
    registry[[owner]]$gene_label_repel <- geometry[!unused, , drop = FALSE]
  }
  registry <- ggchord_detour_crossing_leaders(registry, layout)
  ggchord_measure_external_leader_crossings(registry)
}
