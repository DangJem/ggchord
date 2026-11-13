# Actual directed native paths, shared by placement and measured routing.
ggchord_annotation_leader_segments <- function(registry) {
  parts <- list()
  for (id in names(registry)) {
    d <- registry[[id]]$gene_label_repel
    if (is.data.frame(d) && nrow(d)) {
      at <- which(d$.component %in% 'segment' & is.finite(d$source_row))
      if (length(at)) parts[[length(parts)+1L]] <- data.frame(owner = id,
        component = 'gene_label_repel', source = d$source_row[at],
        x1 = d$x[at], y1 = d$y[at], x2 = d$xend[at], y2 = d$yend[at])
    }
    d <- registry[[id]]$restriction_site
    if (is.data.frame(d) && nrow(d)) {
      at <- which(d$restriction_component %in% 'leader' & is.finite(d$source_row))
      for (source in unique(d$source_row[at])) {
        rows <- at[d$source_row[at] == source]
        if (length(rows) < 2L) next
        parts[[length(parts)+1L]] <- data.frame(owner = id,
          component = 'restriction_site', source = source,
          x1 = head(d$x[rows], -1), y1 = head(d$y[rows], -1),
          x2 = tail(d$x[rows], -1), y2 = tail(d$y[rows], -1))
      }
    }
  }
  all <- ggchord_rbind_fill(parts) %||% data.frame()
  if (nrow(all)) all <- all[stats::complete.cases(all) &
    (all$x1 != all$x2 | all$y1 != all$y2), , drop = FALSE]
  all
}

# Fixed directed roots are obstacles for text placement. Mutable shoulders are
# deliberately absent so rerouting cannot change an axis obstacle afterwards.
ggchord_annotation_leader_roots <- function(registry) {
  segments <- ggchord_annotation_leader_segments(registry)
  roots <- list()
  if (nrow(segments)) {
    keys <- paste(segments$owner, segments$component, segments$source, sep = '\r')
    for (key in unique(keys)) {
      path <- segments[keys == key, , drop = FALSE]
      starts <- vapply(seq_len(nrow(path)), function(i) !any(
        (path$x2-path$x1[i])^2+(path$y2-path$y1[i])^2 < 1e-18), logical(1))
      point <- if (sum(starts) == 1L) c(path$x1[starts], path$y1[starts]) else NULL
      if (is.null(point)) {
        geometry <- registry[[path$owner[1L]]][[path$component[1L]]]
        if (all(c('anchor_x', 'anchor_y') %in% names(geometry))) {
          anchors <- unique(geometry[geometry$source_row %in% path$source[1L],
            c('anchor_x', 'anchor_y'), drop = FALSE])
          anchors <- anchors[stats::complete.cases(anchors), , drop = FALSE]
          if (nrow(anchors) == 1L) point <- unlist(anchors, use.names = FALSE)
        }
      }
      if (length(point) == 2L && all(is.finite(point))) roots[[length(roots)+1L]] <-
        data.frame(x0 = point[1L], y0 = point[2L], x1 = point[1L], y1 = point[2L])
    }
  }
  # External feature names have a fixed radial exit through the backbone.
  # Reserve that short channel, not a mutable shoulder or the entire route,
  # so automatic axis strings do not occupy the exit and cause a feedback loop.
  backbone <- unlist(lapply(registry, function(owner) {
    d <- owner[["seq"]]
    if (!is.data.frame(d) || !nrow(d)) return(numeric())
    sqrt(d$x^2+d$y^2)
  }), use.names = FALSE)
  if (length(backbone)) for (owner in registry) {
    d <- owner$gene_label_repel
    if (!is.data.frame(d) || !"feature_label_mode" %in% names(d)) next
    d <- d[d$.component %in% "text" &
      d$feature_label_mode %in% "external", , drop = FALSE]
    radius <- sqrt(d$anchor_x^2+d$anchor_y^2)
    exit <- max(backbone)+.055
    keep <- is.finite(radius) & radius>0 & radius<max(backbone)
    d <- d[keep, , drop = FALSE]; radius <- radius[keep]
    if (nrow(d)) roots[[length(roots)+1L]] <- data.frame(
      x0 = d$anchor_x, y0 = d$anchor_y,
      x1 = d$anchor_x*exit/radius, y1 = d$anchor_y*exit/radius)
  }
  ggchord_rbind_fill(roots) %||% data.frame(x0 = numeric(), y0 = numeric(),
    x1 = numeric(), y1 = numeric())
}

# Route collided annotation callouts around text and all other leaders.
# Positions, directed chain endpoints and source membership remain unchanged.
ggchord_detour_crossing_leaders <- function(registry, layout) {
  crosses <- function(path, others) nrow(others) > 0L && any(vapply(seq_len(nrow(path)), function(i)
    any(ggchord_route_segment_crossings(c(path$x1[i], path$y1[i]),
      c(path$x2[i], path$y2[i]), others)), logical(1)))
  leaders <- ggchord_annotation_leader_segments(registry)
  if (!nrow(leaders)) return(registry)
  boxes <- ggchord_annotation_text_obstacles(registry, layout)
  units <- layout$text_units_per_inch %||% .35
  clearance <- .25*units/25.4
  extent <- range(c(boxes$xmin, boxes$xmax, boxes$ymin, boxes$ymax,
    leaders$x1, leaders$y1, leaders$x2, leaders$y2), finite = TRUE)
  margin <- max(clearance*4, diff(extent)*.05)
  origin <- rep(extent[1]-margin, 2); canvas <- rep(diff(extent)+2*margin, 2)
  for (owner in names(registry)) for (component in c('gene_label_repel', 'restriction_site')) {
    geometry <- registry[[owner]][[component]]
    if (!is.data.frame(geometry) || !nrow(geometry)) next
    if (component == "gene_label_repel" && "label_layout" %in% names(geometry) &&
        any(geometry$label_layout %in% "feature")) next
    sources <- if (component == 'gene_label_repel')
      unique(geometry$source_row[geometry$.component %in% 'text']) else
      unique(geometry$source_row[geometry$restriction_component %in% 'leader'])
    for (source in sources[!is.na(sources)]) {
      leaders <- ggchord_annotation_leader_segments(registry)
      own <- leaders$owner == owner & leaders$component == component & leaders$source == source
      path <- leaders[own, , drop = FALSE]; others <- leaders[!own, , drop = FALSE]
      if (!nrow(path)) next
      own_text <- boxes$owner == owner & !is.na(boxes$source_row) & boxes$source_row == source
      native_segments <- data.frame(x0 = path$x1, y0 = path$y1, x1 = path$x2, y1 = path$y2)
      text_blocked <- ggchord_inner_leader_hits(native_segments, boxes[!own_text, , drop = FALSE])
      if (!text_blocked && !crosses(path, others)) next
      # Infer directed chain ends from adjacency, rather than guessing from
      # radius, text distance or the storage order of replacement segments.
      starts <- vapply(seq_len(nrow(path)), function(i) !any(
        (path$x2-path$x1[i])^2+(path$y2-path$y1[i])^2 < 1e-18), logical(1))
      ends <- vapply(seq_len(nrow(path)), function(i) !any(
        (path$x1-path$x2[i])^2+(path$y1-path$y2[i])^2 < 1e-18), logical(1))
      status <- if (text_blocked) 'text_unresolved' else 'crossing_unresolved'
      if (sum(starts) == 1L && sum(ends) == 1L) {
        start <- c(path$x1[starts], path$y1[starts]); end <- c(path$x2[ends], path$y2[ends])
        obstacles <- boxes[!own_text, c('xmin', 'xmax', 'ymin', 'ymax'), drop = FALSE]
        obstacles$xmin <- obstacles$xmin-origin[1]; obstacles$xmax <- obstacles$xmax-origin[1]
        obstacles$ymin <- obstacles$ymin-origin[2]; obstacles$ymax <- obstacles$ymax-origin[2]
        barriers <- others
        barriers$x1 <- barriers$x1-origin[1]; barriers$x2 <- barriers$x2-origin[1]
        barriers$y1 <- barriers$y1-origin[2]; barriers$y2 <- barriers$y2-origin[2]
        groups <- boxes$text_object[!own_text]
        active <- which(ggchord_route_segment_hits(start-origin, end-origin, obstacles))
        barrier_groups <- paste(barriers$owner, barriers$component, barriers$source, sep = '\r')
        active_barriers <- which(ggchord_route_segment_crossings(start-origin, end-origin, barriers))
        active_barriers <- which(barrier_groups %in% barrier_groups[active_barriers])
        route <- NULL
        # Add only text objects touched by trial routes. Every accepted route
        # still clears the full obstacle set; distant glyph corners need not
        # inflate the initial visibility graph.
        for (attempt in seq_len(length(unique(groups))+length(unique(barrier_groups))+1L)) {
          route <- ggchord_route_shortest(start-origin, end-origin,
            obstacles[active, , drop = FALSE], canvas,
            clearance = clearance, barriers = barriers[active_barriers, , drop = FALSE])
          if (is.null(route)) break
          touched <- unique(unlist(lapply(seq_len(nrow(route)-1L), function(k)
            which(ggchord_route_segment_hits(route[k, ], route[k+1L, ], obstacles))), use.names = FALSE))
          crossed <- unique(unlist(lapply(seq_len(nrow(route)-1L), function(k)
            which(ggchord_route_segment_crossings(route[k, ], route[k+1L, ], barriers))), use.names = FALSE))
          if (!length(touched) && !length(crossed)) break
          newer <- which(groups %in% groups[union(active, touched)])
          newer_barriers <- which(barrier_groups %in% barrier_groups[union(active_barriers, crossed)])
          route <- NULL
          if (length(newer) <= length(active) && length(newer_barriers) <= length(active_barriers)) break
          active <- newer
          active_barriers <- newer_barriers
        }
        if (!is.null(route)) {
          route <- sweep(route, 2, origin, '+')
          segments <- data.frame(x0 = head(route[, 1], -1), y0 = head(route[, 2], -1),
            x1 = tail(route[, 1], -1), y1 = tail(route[, 2], -1))
          old_length <- sum(sqrt((path$x2-path$x1)^2+(path$y2-path$y1)^2))
          new_length <- sum(sqrt((segments$x1-segments$x0)^2+(segments$y1-segments$y0)^2))
          trial <- data.frame(x1 = segments$x0, y1 = segments$y0,
            x2 = segments$x1, y2 = segments$y1)
          # Recheck all text, including the final label: contact at its edge
          # is allowed, running through its interior is not.
          if (new_length <= old_length*3+5*units/25.4 && !crosses(trial, others) &&
              !ggchord_inner_leader_hits(segments, boxes)) {
            if (component == 'gene_label_repel') {
              at <- which(geometry$.component %in% 'segment' & geometry$source_row == source)
              replacement <- geometry[rep(at[1], nrow(segments)), , drop = FALSE]
              replacement$x <- segments$x0; replacement$y <- segments$y0
              replacement$xend <- segments$x1; replacement$yend <- segments$y1
              for (pair in list(c('x0', 'x'), c('y0', 'y'), c('x1', 'xend'), c('y1', 'yend')))
                if (pair[1] %in% names(replacement)) replacement[[pair[1]]] <- replacement[[pair[2]]]
            } else {
              at <- which(geometry$restriction_component %in% 'leader' & geometry$source_row == source)
              replacement <- geometry[rep(at[1], nrow(route)), , drop = FALSE]
              replacement$x <- route[, 1]; replacement$y <- route[, 2]
            }
            geometry <- ggchord_rbind_fill(list(geometry[-at, , drop = FALSE], replacement))
            status <- if (text_blocked) 'text_detour' else 'crossing_detour'
          }
        }
      }
      if (!'leader_route_status' %in% names(geometry)) geometry$leader_route_status <- NA_character_
      geometry$leader_route_status[!is.na(geometry$source_row) & geometry$source_row == source] <- status
      registry[[owner]][[component]] <- geometry
    }
  }
  registry
}
