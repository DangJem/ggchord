# Global exterior labeling for circular maps.
#
# Every candidate is atomic: text box, attachment port and complete leader are
# selected together.  Annotation classes retain their own data and rendering,
# while the coordinate system owns the shared physical solve.

ggchord_outer_task_columns <- function(geometry) {
  for (nm in c("outer_annotation_region", "outer_spill_reason")) {
    if (!nm %in% names(geometry)) geometry[[nm]] <- NA_character_
  }
  for (nm in c("outer_track", "outer_slot", "dominant_outer_band",
      "leader_crossing_count")) {
    if (!nm %in% names(geometry)) geometry[[nm]] <- NA_integer_
  }
  for (nm in c("label_attachment_x", "label_attachment_y")) {
    if (!nm %in% names(geometry)) geometry[[nm]] <- NA_real_
  }
  if (!"outer_text_units_per_inch" %in% names(geometry)) {
    geometry$outer_text_units_per_inch <- NA_real_
  }
  geometry
}

ggchord_outer_collect_tasks <- function(registry, layout) {
  tasks <- list()
  sequence_ids <- unique(c(
    names(layout$seq_arcs),
    if (is.data.frame(layout$seq_data) && "accver" %in% names(layout$seq_data))
      as.character(layout$seq_data$accver) else character()
  ))
  sequence_ids <- sequence_ids[!is.na(sequence_ids) & nzchar(sequence_ids)]

  add <- function(id, component, row, class, geometry, root,
                  source, anchor_position) {
    accver <- if ("accver" %in% names(geometry)) {
      as.character(geometry$accver[row])
    } else NA_character_
    if (length(accver) != 1L || is.na(accver) || !nzchar(accver)) {
      accver <- if (length(sequence_ids) == 1L) sequence_ids else id
    }
    if (length(anchor_position) != 1L || !is.finite(anchor_position)) {
      anchor_position <- row
    }
    tasks[[length(tasks) + 1L]] <<- list(
      id = id, component = component, row = row, class = class,
      source = source, accver = accver,
      root = as.numeric(root), angle = atan2(root[2L], root[1L]),
      anchor_position = anchor_position
    )
  }

  for (id in names(registry)) {
    sites <- registry[[id]]$restriction_site
    if (is.data.frame(sites) && nrow(sites)) {
      sites <- ggchord_outer_task_columns(sites)
      composite_text <- rep(FALSE, nrow(sites))
      for (column in c("enzyme_label", "coordinate_label")) {
        if (column %in% names(sites)) {
          composite_text <- composite_text |
            (!is.na(sites[[column]]) & nzchar(sites[[column]]))
        }
      }
      rows <- which(sites$.component %in% "label" & (
        (!is.na(sites$label) & nzchar(sites$label)) | composite_text
      ))
      for (row in rows) {
        source <- as.integer(sites$source_row[row])
        path_rows <- which(sites$restriction_component %in% "leader" &
          is.finite(sites$source_row) & sites$source_row == source)
        if (!length(path_rows)) {
          path_rows <- which(sites$restriction_component %in% "tick" &
            is.finite(sites$source_row) & sites$source_row == source)
        }
        root <- if (length(path_rows)) {
          c(sites$x[path_rows[1L]], sites$y[path_rows[1L]])
        } else c(sites$x[row], sites$y[row])
        add(id, "restriction_site", row, "restriction", sites, root,
          source, sites$anchor_position[row] %||% row)
      }
      registry[[id]]$restriction_site <- sites
    }

    labels <- registry[[id]]$gene_label_repel
    if (!is.data.frame(labels) || !nrow(labels)) next
    labels <- ggchord_outer_task_columns(labels)
    external <- if ("feature_label_mode" %in% names(labels)) {
      labels$feature_label_mode %in% "external"
    } else rep(FALSE, nrow(labels))
    if ("annotation_class" %in% names(labels)) {
      external <- external | labels$annotation_class %in% "primer"
    }
    rows <- which(labels$.component %in% "text" & external &
      !is.na(labels$label) & nzchar(labels$label))
    for (row in rows) {
      source <- as.integer(labels$source_row[row])
      segment_rows <- which(labels$.component %in% "segment" &
        is.finite(labels$source_row) & labels$source_row == source)
      root <- if (length(segment_rows)) {
        first <- segment_rows[which.max(
          (labels$x[segment_rows] - labels$text_x[row])^2 +
            (labels$y[segment_rows] - labels$text_y[row])^2
        )]
        c(labels$x[first], labels$y[first])
      } else if (all(c("anchor_x", "anchor_y") %in% names(labels))) {
        c(labels$anchor_x[row], labels$anchor_y[row])
      } else c(labels$text_x[row], labels$text_y[row])
      class <- if ("annotation_class" %in% names(labels) &&
          identical(as.character(labels$annotation_class[row]), "primer")) {
        "primer"
      } else "feature"
      add(id, "gene_label_repel", row, class, labels, root, source,
        labels$anchor_position[row] %||% row)
    }
    registry[[id]]$gene_label_repel <- labels
  }
  list(registry = registry, tasks = tasks)
}

ggchord_outer_task_row <- function(task, registry) {
  geometry <- registry[[task$id]][[task$component]]
  row <- geometry[task$row, , drop = FALSE]
  rendered_text <- if (identical(task$class, "restriction") &&
      all(c("enzyme_label", "coordinate_label") %in% names(row))) {
    enzyme <- ifelse(is.na(row$enzyme_label[1L]), "",
      as.character(row$enzyme_label[1L]))
    coordinate <- ifelse(is.na(row$coordinate_label[1L]), "",
      as.character(row$coordinate_label[1L]))
    parts <- if (identical(as.character(row$label_order[1L]),
        "position_enzyme")) c(coordinate, enzyme) else c(enzyme, coordinate)
    paste(parts[nzchar(parts)], collapse = " ")
  } else if (identical(task$class, "primer") &&
      "feature_label" %in% names(row) &&
      !is.na(row$feature_label[1L]) && nzchar(row$feature_label[1L])) {
    as.character(row$feature_label[1L])
  } else if ("text" %in% names(row) && !is.na(row$text[1L]) &&
      nzchar(row$text[1L])) {
    as.character(row$text[1L])
  } else as.character(row$label[1L])
  row$text <- rendered_text
  row$text_x <- row$text_x %||% row$x
  row$text_y <- row$text_y %||% row$y
  row$text_angle <- 0
  row$angle <- 0
  row$hjust <- .5
  row$vjust <- .5
  row
}

ggchord_outer_translate_box <- function(box, x, y) {
  dx <- x - box$cx[1L]
  dy <- y - box$cy[1L]
  for (nm in c("x", "cx", "xmin", "xmax")) box[[nm]] <- box[[nm]] + dx
  for (nm in c("y", "cy", "ymin", "ymax")) box[[nm]] <- box[[nm]] + dy
  box
}

ggchord_outer_primer_text <- function(row, left = FALSE) {
  if (!all(c("primer_name", "primer_start", "primer_end") %in% names(row)) ||
      is.na(row$primer_name[1L])) return(as.character(row$text[1L]))
  range <- paste0("(", row$primer_start[1L], " .. ",
    row$primer_end[1L], ")")
  if (left) paste(range, row$primer_name[1L]) else
    paste(row$primer_name[1L], range)
}

ggchord_outer_measure_task <- function(task, row, units_per_inch,
                                       box_padding = .045) {
  measured <- ggchord_text_boxes(
    row, units_per_inch = units_per_inch, box_padding = box_padding
  )
  if (!identical(task$class, "restriction") ||
      !all(c("enzyme_label", "coordinate_label") %in% names(row))) {
    return(measured)
  }
  size <- row$size %||% 2.9
  enzyme <- as.character(row$enzyme_label[1L] %||% "")
  coordinate <- as.character(row$coordinate_label[1L] %||% "")
  enzyme_face <- row$enzyme_fontface %||% "plain"
  enzyme_metrics <- ggchord_restriction_text_metrics(
    enzyme, size, units_per_inch, parse = FALSE, fontface = enzyme_face
  )
  coordinate_metrics <- ggchord_restriction_text_metrics(
    coordinate, size, units_per_inch, parse = FALSE, fontface = "plain"
  )
  gap <- if (!is.na(enzyme) && nzchar(enzyme) &&
      !is.na(coordinate) && nzchar(coordinate)) {
    .65 / 25.4 * units_per_inch
  } else 0
  measured$w <- enzyme_metrics$width + coordinate_metrics$width + gap
  measured$h <- max(enzyme_metrics$height, coordinate_metrics$height)
  padding <- 2 * box_padding * units_per_inch
  measured$ow <- measured$w + padding
  measured$oh <- measured$h + padding
  measured$bw <- measured$ow
  measured$bh <- measured$oh
  measured$xmin <- measured$cx - measured$bw / 2
  measured$xmax <- measured$cx + measured$bw / 2
  measured$ymin <- measured$cy - measured$bh / 2
  measured$ymax <- measured$cy + measured$bh / 2
  measured
}

ggchord_outer_min_segment_radius <- function(segment) {
  p0 <- c(segment$x1, segment$y1)
  p1 <- c(segment$x2, segment$y2)
  delta <- p1 - p0
  denominator <- sum(delta^2)
  parameter <- if (!is.finite(denominator) || denominator <= 1e-16) 0 else
    max(0, min(1, -sum(p0 * delta) / denominator))
  sqrt(sum((p0 + parameter * delta)^2))
}

ggchord_outer_name_port <- function(task, row, full_box, root,
                                    units_per_inch, side_override = NULL) {
  centre <- c(full_box$cx[1L], full_box$cy[1L])
  if (!identical(task$class, "primer") ||
      !"primer_name" %in% names(row) || is.na(row$primer_name[1L]) ||
      !nzchar(row$primer_name[1L])) {
    direction <- root - centre
    half <- c(full_box$w[1L], full_box$h[1L]) / 2
    ratios <- half / pmax(abs(direction), 1e-12)
    return(centre + min(ratios) * direction)
  }
  name_row <- row
  name_row$text <- as.character(row$primer_name[1L])
  name_box <- ggchord_text_boxes(
    name_row, units_per_inch = units_per_inch, box_padding = .045
  )
  # The renderer keeps the name next to the plasmid: right-side labels use
  # "name (range)", left-side labels use "(range) name".  In both cases the
  # coordinate fragment is excluded from the legal attachment interval.
  name_offset <- (full_box$w[1L] - name_box$w[1L]) / 2
  name_cx <- centre[1L] + if (centre[1L] < 0) name_offset else -name_offset
  side <- side_override %||% if (root[1L] <= name_cx) -1 else 1
  c(name_cx + side * name_box$w[1L] / 2, centre[2L])
}

ggchord_outer_route <- function(root, port, backbone_outer) {
  direct <- data.frame(
    x1 = root[1L], y1 = root[2L], x2 = port[1L], y2 = port[2L],
    move_endpoint = FALSE, attachment_valid = TRUE
  )
  root_angle <- atan2(root[2L], root[1L])
  for (extra in c(.018, .04, .075, .12)) {
    elbow <- (backbone_outer + extra) * c(cos(root_angle), sin(root_angle))
    route <- data.frame(
      x1 = c(root[1L], elbow[1L]),
      y1 = c(root[2L], elbow[2L]),
      x2 = c(elbow[1L], port[1L]),
      y2 = c(elbow[2L], port[2L]),
      move_endpoint = FALSE, attachment_valid = TRUE
    )
    safe <- vapply(seq_len(nrow(route)), function(i) {
      ggchord_outer_min_segment_radius(route[i, , drop = FALSE]) >=
        backbone_outer - .003
    }, logical(1L))
    if (all(safe)) return(route)
  }
  direct
}

ggchord_outer_candidates <- function(task, registry, units_per_inch,
                                     backbone_outer, preferred_unwrapped,
                                     natural_unwrapped,
                                     max_tracks = 12L) {
  row <- ggchord_outer_task_row(task, registry)
  row_left <- row
  if (identical(task$class, "primer")) {
    row$text <- ggchord_outer_primer_text(row, left = FALSE)
    measured_right <- ggchord_outer_measure_task(
      task, row, units_per_inch, box_padding = .045
    )
    measured_left <- ggchord_outer_measure_task(
      task, row_left, units_per_inch, box_padding = .045
    )
  } else {
    measured_right <- ggchord_outer_measure_task(
      task, row, units_per_inch, box_padding = .045
    )
    measured_left <- measured_right
  }
  shifts <- if (identical(task$class, "primer")) {
    c(0, .06, -.06, .13, -.13, .22, -.22)
  } else if (identical(task$class, "restriction")) {
    c(0, .07, -.07, .14, -.14)
  } else c(0, .05, -.05, .10, -.10, .15, -.15, .22, -.22,
    .30, -.30)
  port_sides <- if (identical(task$class, "primer")) c(-1, 1) else NA_real_
  track_count <- max_tracks
  candidates <- vector(
    "list", track_count * length(shifts) * length(port_sides) + 1L
  )
  index <- 0L
  for (track in seq_len(track_count)) for (shift in shifts) for (
      port_side in port_sides) {
    control_angle <- preferred_unwrapped + shift
    direction <- c(cos(control_angle), sin(control_angle))
    measured <- if (identical(task$class, "primer") &&
        direction[1L] < 0) measured_left else measured_right
    candidate_row <- if (identical(task$class, "primer") &&
        direction[1L] < 0) row_left else row
    radial_half <- abs(direction[1L]) * measured$w[1L] / 2 +
      abs(direction[2L]) * measured$h[1L] / 2
    control_radius <- backbone_outer + .035 + (track - 1L) * .052
    centre_radius <- control_radius + radial_half
    centre <- centre_radius * direction
    box <- ggchord_outer_translate_box(measured, centre[1L], centre[2L])
    port <- ggchord_outer_name_port(
      task, candidate_row, box, task$root, units_per_inch,
      side_override = if (is.finite(port_side)) port_side else NULL
    )
    route <- ggchord_outer_route(task$root, port, backbone_outer)
    leader_length <- sum(sqrt(
      (route$x2 - route$x1)^2 + (route$y2 - route$y1)^2
    ))
    index <- index + 1L
    candidates[[index]] <- list(
      hidden = FALSE, row = candidate_row, box = box, route = route,
      centre = centre, port = port, angle = control_angle,
      track = track, slot = match(shift, shifts) - 1L,
      units_per_inch = units_per_inch,
      cost = (track - 1L)^2 * .42 +
        ((control_angle - natural_unwrapped) / .035)^2 * .055 +
        leader_length * .32 + (nrow(route) - 1L) * .18
    )
  }
  if (identical(task$class, "primer") || isTRUE(task$required)) {
    return(candidates[seq_len(index)])
  }
  index <- index + 1L
  candidates[[index]] <- list(
    hidden = TRUE, row = row, box = data.frame(), route = data.frame(),
    centre = c(NA_real_, NA_real_), port = c(NA_real_, NA_real_),
    angle = natural_unwrapped, track = NA_integer_, slot = NA_integer_,
    units_per_inch = units_per_inch,
    # Hiding is lexicographically more expensive than any visible layout.
    # Feature display priority can later refine this common fallback.
    cost = 1e5
  )
  candidates
}

ggchord_outer_route_hits_box <- function(route, box, tolerance = 1e-8) {
  if (!nrow(route) || !nrow(box)) return(FALSE)
  xmin <- box$xmin[1L] + tolerance
  xmax <- box$xmax[1L] - tolerance
  ymin <- box$ymin[1L] + tolerance
  ymax <- box$ymax[1L] - tolerance
  for (i in seq_len(nrow(route))) {
    x0 <- route$x1[i]; y0 <- route$y1[i]
    dx <- route$x2[i] - x0; dy <- route$y2[i] - y0
    p <- c(-dx, dx, -dy, dy)
    q <- c(x0 - xmin, xmax - x0, y0 - ymin, ymax - y0)
    lower <- 0; upper <- 1; hit <- TRUE
    for (edge in seq_len(4L)) {
      if (abs(p[edge]) <= tolerance) {
        if (q[edge] < 0) { hit <- FALSE; break }
      } else {
        ratio <- q[edge] / p[edge]
        if (p[edge] < 0) lower <- max(lower, ratio) else
          upper <- min(upper, ratio)
        if (lower > upper) { hit <- FALSE; break }
      }
    }
    if (hit && upper > tolerance && lower < 1 - tolerance) return(TRUE)
  }
  FALSE
}

ggchord_outer_candidate_compatible <- function(candidate, state) {
  if (isTRUE(candidate$hidden)) return(TRUE)
  if (is.finite(state$last_angle) && candidate$angle <= state$last_angle) {
    return(FALSE)
  }
  for (chosen in state$chosen) {
    if (isTRUE(chosen$hidden)) next
    if (candidate$box$xmin < chosen$box$xmax - 1e-7 &&
        candidate$box$xmax > chosen$box$xmin + 1e-7 &&
        candidate$box$ymin < chosen$box$ymax - 1e-7 &&
        candidate$box$ymax > chosen$box$ymin + 1e-7) {
      return(FALSE)
    }
    for (i in seq_len(nrow(candidate$route))) {
      for (j in seq_len(nrow(chosen$route))) {
        if (ggchord_external_leader_crosses(
            candidate$route[i, ], chosen$route[j, ])) return(FALSE)
      }
    }
    if (ggchord_outer_route_hits_box(candidate$route, chosen$box) ||
        ggchord_outer_route_hits_box(chosen$route, candidate$box)) {
      return(FALSE)
    }
  }
  TRUE
}

ggchord_outer_solve_block <- function(tasks, preferred, natural, registry,
                                      units_per_inch, backbone_outer,
                                      beam_width = 512L,
                                      occupied = list()) {
  candidate_sets <- lapply(seq_along(tasks), function(i) {
    ggchord_outer_candidates(
      tasks[[i]], registry, units_per_inch, backbone_outer,
      preferred[i], natural[i]
    )
  })
  occupied_angles <- vapply(occupied, function(candidate) {
    if (candidate$hidden) NA_real_ else candidate$angle
  }, numeric(1L))
  occupied_angles <- occupied_angles[is.finite(occupied_angles)]
  fixed_state <- list(
    chosen = occupied,
    last_angle = if (length(occupied_angles)) max(occupied_angles) else -Inf
  )
  candidate_sets <- lapply(candidate_sets, function(candidates) {
    keep <- vapply(candidates, ggchord_outer_candidate_compatible,
      logical(1L), state = fixed_state)
    candidates[keep]
  })
  if (any(!lengths(candidate_sets))) {
    ggchord_stop("Circular external-label solver exhausted a candidate set")
  }

  # Geometry is immutable inside one solve. Cache only candidate pairs that
  # the beam visits; distant pairs are audited on the complete solution.
  pair_matrix <- function(previous, current) matrix(
    NA, nrow = length(candidate_sets[[previous]]),
    ncol = length(candidate_sets[[current]])
  )
  compatibility <- vector("list", length(tasks))
  for (i in seq_along(tasks)) {
    compatibility[[i]] <- vector("list", i - 1L)
    if (i <= 1L) next
    for (previous in seq_len(i - 1L)) {
      if (natural[i] - natural[previous] > .65) {
        compatibility[[i]][[previous]] <- TRUE
        next
      }
      compatibility[[i]][[previous]] <- pair_matrix(previous, i)
    }
  }

  pair_compatible <- function(previous, current, previous_index,
                              current_index) {
    cache <- compatibility[[current]][[previous]]
    if (length(cache) == 1L && isTRUE(cache)) return(TRUE)
    value <- cache[previous_index, current_index]
    if (!is.na(value)) return(value)
    prior <- candidate_sets[[previous]][[previous_index]]
    state <- list(
      chosen = list(prior),
      last_angle = if (prior$hidden) -Inf else prior$angle
    )
    value <- ggchord_outer_candidate_compatible(
      candidate_sets[[current]][[current_index]], state
    )
    cache[previous_index, current_index] <- value
    compatibility[[current]][[previous]] <<- cache
    value
  }

  run_search <- function() {
    states <- list(list(indices = integer(), cost = 0, hidden = 0L,
      feature_hidden = 0L))
    failed_at <- NA_integer_
    # Most distant pairs are provisionally compatible and need no lookup in
    # the inner beam loop. The complete-solution audit below can promote any
    # distant conflict to a matrix before the next search pass.
    active_previous <- lapply(seq_along(tasks), function(current) {
      if (current <= 1L) return(integer())
      which(!vapply(compatibility[[current]], function(cache) {
        length(cache) == 1L && isTRUE(cache)
      }, logical(1L)))
    })
    for (i in seq_along(tasks)) {
      expanded <- list()
      for (state in states) for (candidate_index in
          seq_along(candidate_sets[[i]])) {
        compatible <- TRUE
        for (previous in active_previous[[i]]) {
          if (!pair_compatible(previous, i, state$indices[previous],
              candidate_index)) {
            compatible <- FALSE
            break
          }
        }
        if (!compatible) next
        candidate <- candidate_sets[[i]][[candidate_index]]
        next_state <- list(
          indices = c(state$indices, candidate_index),
          cost = state$cost + candidate$cost,
          hidden = state$hidden + as.integer(candidate$hidden),
          feature_hidden = state$feature_hidden +
            as.integer(candidate$hidden &&
              identical(tasks[[i]]$class, "feature"))
        )
        expanded[[length(expanded) + 1L]] <- next_state
      }
      if (!length(expanded)) {
        failed_at <- i
        break
      }
      hidden_count <- vapply(expanded, `[[`, integer(1L), "hidden")
      feature_hidden <- vapply(expanded, `[[`, integer(1L),
        "feature_hidden")
      state_cost <- vapply(expanded, `[[`, numeric(1L), "cost")
      # Preserve several partial solutions at every visibility level.  A
      # purely lexicographic prefix would discard a state that hides one early
      # low-priority restriction even when that is the only state capable of
      # accommodating a mandatory primer encountered later in circular order.
      buckets <- split(seq_along(expanded), paste(feature_hidden,
        hidden_count, sep = ":"))
      quota <- max(1L, floor(beam_width / length(buckets)))
      keep <- unlist(lapply(buckets, function(rows) {
        rows[order(state_cost[rows])][seq_len(min(length(rows), quota))]
      }), use.names = FALSE)
      expanded <- expanded[keep]
      keys <- vapply(expanded, function(state) {
        paste(tail(state$indices, 4L), collapse = ".")
      }, character(1L))
      states <- expanded[!duplicated(keys)]
    }
    if (!length(states) || length(states[[1L]]$indices) != length(tasks)) {
      detail <- if (is.finite(failed_at)) paste0(" at label ", failed_at,
        " (`", tasks[[failed_at]]$class, "`, source ",
        tasks[[failed_at]]$source, ")") else ""
      ggchord_stop(
        "Circular external-label solver failed to construct a complete state",
        detail
      )
    }
    states[[order(
      vapply(states, `[[`, integer(1L), "feature_hidden"),
      vapply(states, `[[`, integer(1L), "hidden"),
      vapply(states, `[[`, numeric(1L), "cost")
    )[1L]]]
  }

  best <- NULL
  chosen <- NULL
  for (iteration in seq_len(8L)) {
    best <- run_search()
    chosen <- lapply(seq_along(tasks), function(i) {
      candidate_sets[[i]][[best$indices[i]]]
    })
    conflicts <- list()
    if (length(chosen) > 1L) for (current in 2:length(chosen)) {
      for (previous in seq_len(current - 1L)) {
        state <- list(
          chosen = list(chosen[[previous]]),
          last_angle = if (chosen[[previous]]$hidden) -Inf else
            chosen[[previous]]$angle
        )
        if (!ggchord_outer_candidate_compatible(chosen[[current]], state)) {
          conflicts[[length(conflicts) + 1L]] <- c(previous, current)
        }
      }
    }
    if (!length(conflicts)) break
    for (pair in conflicts) {
      compatibility[[pair[2L]]][[pair[1L]]] <-
        pair_matrix(pair[1L], pair[2L])
    }
  }
  if (length(conflicts)) {
    ggchord_stop(
      "Circular external-label solver did not converge to a planar state"
    )
  }
  list(tasks = tasks, candidates = chosen,
    hidden = best$hidden, cost = best$cost)
}

ggchord_outer_solve_sequence <- function(tasks, registry, units_per_inch,
                                         backbone_outer) {
  angles <- vapply(tasks, `[[`, numeric(1L), "angle")
  ties <- vapply(tasks, `[[`, numeric(1L), "anchor_position")
  order <- ggchord_open_circular_order(angles, ties)
  mandatory <- which(vapply(tasks[order], function(task) {
    task$class %in% c("primer", "feature")
  }, logical(1L)))
  if (length(mandatory) && mandatory[1L] > 1L) {
    cut <- mandatory[1L]
    order <- c(order[seq.int(cut, length(order))], order[seq_len(cut - 1L)])
  }
  tasks <- tasks[order]
  unwrapped <- angles[order]
  if (length(unwrapped) > 1L) for (i in 2:length(unwrapped)) {
    # Equal anchors are common for enzymes sharing one cut site.  They belong
    # to one local fan and must not be mistaken for a wrap across the circular
    # origin; the ordered packing pass below gives them distinct controls.
    while (unwrapped[i] < unwrapped[i - 1L] - 1e-10) {
      unwrapped[i] <- unwrapped[i] + 2 * pi
    }
  }
  # Reserve tangential space for every label before choosing radial levels.
  # This is a global contour preference, not a final placement: the atomic
  # candidate solve below may still move a port or change its radial level.
  measured <- lapply(tasks, function(task) {
    row <- ggchord_outer_task_row(task, registry)
    ggchord_outer_measure_task(
      task, row, units_per_inch, box_padding = .045
    )
  })
  measured <- ggchord_rbind_fill(measured)
  packing_radius <- backbone_outer + .10
  preferred <- unwrapped
  for (iteration in seq_len(3L)) {
    half_extent <- abs(sin(preferred)) * measured$w / 2 +
      abs(cos(preferred)) * measured$h / 2
    preferred <- ggchord_pack_label_axis(
      unwrapped * packing_radius, seq_along(tasks),
      half_extent, half_extent, gap = .012
    ) / packing_radius
  }
  preferred <- pmax(unwrapped - .35, pmin(unwrapped + .35, preferred))
  feature_rows <- which(vapply(tasks, function(task) {
    identical(task$class, "feature")
  }, logical(1L)))
  if (length(feature_rows)) {
    # Repeated short features form a local callout fan.  Packing every class
    # from the circular cut forward accumulates displacement from unrelated
    # restriction labels, eventually pinning the entire feature fan against
    # one angular limit.  Give each nearby feature run a centred, ordered
    # contour preference; the final compatibility solve still includes all
    # annotation classes and their complete geometry.
    feature_groups <- cumsum(c(1L,
      as.integer(diff(unwrapped[feature_rows]) > .12)))
    for (group in unique(feature_groups)) {
      members <- feature_rows[feature_groups == group]
      if (length(members) < 2L) next
      if (length(members) >= 4L) for (member in members) {
        tasks[[member]]$required <- TRUE
      }
      half_extent <- abs(sin(unwrapped[members])) *
        measured$w[members] / 2 + abs(cos(unwrapped[members])) *
        measured$h[members] / 2
      local <- ggchord_pack_label_axis(
        unwrapped[members] * packing_radius, seq_along(members),
        half_extent, half_extent, gap = .012
      ) / packing_radius
      local <- local - mean(local - unwrapped[members])
      preferred[members] <- pmax(unwrapped[members] - .50,
        pmin(unwrapped[members] + .50, local))
    }
  }
  # Try to preserve dense, closely spaced feature runs as a whole. If the
  # complete joint problem cannot satisfy that extra constraint, fall back to
  # the ordinary soft visibility preference instead of failing the plot.
  solve <- function(task_set) ggchord_outer_solve_block(
    task_set, preferred, unwrapped, registry,
    units_per_inch, backbone_outer
  )
  tryCatch(solve(tasks), error = function(error) {
    has_required <- any(vapply(tasks, function(task) isTRUE(task$required),
      logical(1L)))
    solver_error <- grepl("^Circular external-label solver ",
      conditionMessage(error))
    if (!has_required || !solver_error) stop(error)
    relaxed <- lapply(tasks, function(task) {
      task$required <- FALSE
      task
    })
    solve(relaxed)
  })
}

ggchord_outer_hide_feature <- function(geometry, row) {
  source <- geometry$source_row[row]
  geometry$label[row] <- NA_character_
  if ("text" %in% names(geometry)) geometry$text[row] <- NA_character_
  segments <- geometry$.component %in% "segment" &
    is.finite(geometry$source_row) & geometry$source_row == source
  geometry$xend[segments] <- geometry$x[segments]
  geometry$yend[segments] <- geometry$y[segments]
  geometry$outer_spill_reason[c(row, which(segments))] <-
    "hidden_no_feasible_outer_layout"
  geometry
}

ggchord_outer_apply_restriction <- function(geometry, row, candidate) {
  source <- geometry$source_row[row]
  if (candidate$hidden) {
    geometry$label[row] <- NA_character_
    if ("text" %in% names(geometry)) geometry$text[row] <- NA_character_
    if ("enzyme_label" %in% names(geometry)) {
      geometry$enzyme_label[row] <- NA_character_
    }
    if ("coordinate_label" %in% names(geometry)) {
      geometry$coordinate_label[row] <- NA_character_
    }
    path_rows <- which(geometry$restriction_component %in% "leader" &
      is.finite(geometry$source_row) & geometry$source_row == source)
    if (length(path_rows) && ".path_visible" %in% names(geometry)) {
      geometry$.path_visible[path_rows] <- FALSE
    }
    geometry$outer_spill_reason[c(row, path_rows)] <-
      "hidden_no_feasible_outer_layout"
    return(geometry)
  }
  geometry$x[row] <- candidate$centre[1L]
  geometry$y[row] <- candidate$centre[2L]
  geometry$label[row] <- candidate$row$text
  if ("text" %in% names(geometry)) geometry$text[row] <- candidate$row$text
  if ("enzyme_label" %in% names(geometry)) {
    geometry$enzyme_label[row] <- candidate$row$enzyme_label[1L]
  }
  if ("coordinate_label" %in% names(geometry)) {
    geometry$coordinate_label[row] <- candidate$row$coordinate_label[1L]
  }
  if ("hjust" %in% names(geometry)) geometry$hjust[row] <- .5
  if ("vjust" %in% names(geometry)) geometry$vjust[row] <- .5
  if ("angle" %in% names(geometry)) geometry$angle[row] <- 0
  geometry$outer_annotation_region[row] <- "restriction"
  geometry$outer_track[row] <- candidate$track
  geometry$outer_slot[row] <- candidate$slot
  geometry$dominant_outer_band[row] <- 1L
  geometry$outer_spill_reason[row] <- if (candidate$track > 1L) {
    "global_contour_stack"
  } else NA_character_
  geometry$leader_crossing_count[row] <- 0L
  geometry$outer_text_units_per_inch[row] <- candidate$units_per_inch
  path_rows <- which(geometry$restriction_component %in% "leader" &
    is.finite(geometry$source_row) & geometry$source_row == source)
  if (length(path_rows) && ".path_visible" %in% names(geometry)) {
    geometry$.path_visible[path_rows] <- TRUE
  }
  points <- rbind(
    data.frame(x = candidate$route$x1[1L], y = candidate$route$y1[1L]),
    data.frame(x = candidate$route$x2, y = candidate$route$y2)
  )
  if (nrow(points) > length(path_rows)) {
    template <- geometry[if (length(path_rows)) tail(path_rows, 1L) else row,
      , drop = FALSE]
    extra <- template[rep(1L, nrow(points) - length(path_rows)), , drop = FALSE]
    geometry <- rbind(geometry, extra)
    path_rows <- c(path_rows,
      seq.int(nrow(geometry) - nrow(extra) + 1L, nrow(geometry)))
  }
  use <- seq_len(min(length(path_rows), nrow(points)))
  geometry$x[path_rows[use]] <- points$x[use]
  geometry$y[path_rows[use]] <- points$y[use]
  if (length(path_rows) > length(use)) {
    unused <- path_rows[-use]
    geometry$x[unused] <- tail(points$x, 1L)
    geometry$y[unused] <- tail(points$y, 1L)
  }
  geometry$outer_annotation_region[path_rows] <- "restriction_leader"
  geometry$outer_track[path_rows] <- candidate$track
  geometry$outer_slot[path_rows] <- candidate$slot
  geometry$dominant_outer_band[path_rows] <- 1L
  geometry$outer_spill_reason[path_rows] <- geometry$outer_spill_reason[row]
  geometry$leader_crossing_count[path_rows] <- 0L
  if ("label_attachment_x" %in% names(geometry)) {
    geometry$label_attachment_x[row] <- candidate$port[1L]
    geometry$label_attachment_y[row] <- candidate$port[2L]
  }
  geometry
}

ggchord_outer_apply_solution <- function(registry, solution) {
  for (i in seq_along(solution$tasks)) {
    task <- solution$tasks[[i]]
    candidate <- solution$candidates[[i]]
    geometry <- registry[[task$id]][[task$component]]
    if (identical(task$component, "restriction_site")) {
      geometry <- ggchord_outer_apply_restriction(
        geometry, task$row, candidate
      )
    } else if (candidate$hidden) {
      geometry <- ggchord_outer_hide_feature(geometry, task$row)
    } else {
      translated <- candidate$row
      translated$text_x <- candidate$centre[1L]
      translated$text_y <- candidate$centre[2L]
      placed <- list(
        row = translated, box = candidate$box,
        band = candidate$track, slot = candidate$slot,
        collision = FALSE, box_collision = FALSE,
        leader_collision = FALSE, leader_segments = candidate$route
      )
      geometry <- ggchord_shift_feature_callout(
        geometry, task$row, placed
      )
      # Candidate boxes are measured about their centre.  Preserve that exact
      # rendering contract when writing the solution back; carrying an older
      # side-specific justification would shift the visible glyphs outside
      # the box that participated in collision detection.
      if ("hjust" %in% names(geometry)) geometry$hjust[task$row] <- .5
      if ("vjust" %in% names(geometry)) geometry$vjust[task$row] <- .5
      if ("text" %in% names(geometry)) {
        geometry$text[task$row] <- candidate$row$text
      }
      if ("label" %in% names(geometry)) {
        geometry$label[task$row] <- candidate$row$text
      }
      if ("text_angle" %in% names(geometry)) {
        geometry$text_angle[task$row] <- 0
      }
      if ("angle" %in% names(geometry)) geometry$angle[task$row] <- 0
      geometry$label_attachment_x[task$row] <- candidate$port[1L]
      geometry$label_attachment_y[task$row] <- candidate$port[2L]
      geometry$outer_text_units_per_inch[task$row] <-
        candidate$units_per_inch
      geometry$outer_spill_reason[task$row] <- if (
        candidate$track > 1L) "global_contour_stack" else NA_character_
    }
    registry[[task$id]][[task$component]] <- geometry
  }
  # The solver's atomic routes must survive layer-specific geometry writeback
  # exactly; otherwise a planar candidate solution would not describe what is
  # actually rendered.
  for (i in seq_along(solution$tasks)) {
    task <- solution$tasks[[i]]
    candidate <- solution$candidates[[i]]
    if (candidate$hidden) next
    geometry <- registry[[task$id]][[task$component]]
    if (identical(task$component, "restriction_site")) {
      rows <- which(geometry$restriction_component %in% "leader" &
        is.finite(geometry$source_row) & geometry$source_row == task$source)
      actual <- if (length(rows) >= 2L) data.frame(
        x1 = geometry$x[rows[-length(rows)]],
        y1 = geometry$y[rows[-length(rows)]],
        x2 = geometry$x[rows[-1L]], y2 = geometry$y[rows[-1L]]
      ) else data.frame()
    } else {
      rows <- which(geometry$.component %in% "segment" &
        is.finite(geometry$source_row) & geometry$source_row == task$source)
      actual <- geometry[rows, c("x", "y", "xend", "yend"), drop = FALSE]
      names(actual) <- c("x1", "y1", "x2", "y2")
    }
    keep <- if (nrow(actual)) sqrt(
      (actual$x2 - actual$x1)^2 + (actual$y2 - actual$y1)^2
    ) > 1e-10 else logical()
    actual <- actual[keep, , drop = FALSE]
    expected <- candidate$route[, c("x1", "y1", "x2", "y2"), drop = FALSE]
    if (nrow(actual) != nrow(expected) ||
        (nrow(actual) && max(abs(as.matrix(actual) - as.matrix(expected))) >
          1e-8)) {
      ggchord_stop(
        "External-label route changed during geometry writeback for `",
        task$class, "` source ", task$source
      )
    }
  }
  registry
}

ggchord_share_external_annotations <- function(registry, layout) {
  if (!isTRUE(layout$circular) || !length(registry)) return(registry)
  collected <- ggchord_outer_collect_tasks(registry, layout)
  registry <- collected$registry
  if (!length(collected$tasks)) return(registry)
  # Exterior labels are rendered in the final panel, whose conservative
  # device-space conversion is .35 data units per inch.  An earlier layout
  # estimate may be smaller; using it here would make the solved boxes narrower
  # than the glyphs that are subsequently drawn and registered.
  units_per_inch <- max(layout$text_units_per_inch %||% .30, .35)
  backbone_outer <- if (is.data.frame(layout$backbone_bounds) &&
      nrow(layout$backbone_bounds)) {
    max(layout$backbone_bounds$outer, na.rm = TRUE) + .040
  } else 1.04
  sequence <- vapply(collected$tasks, `[[`, character(1L), "accver")
  for (id in unique(sequence)) {
    solution <- ggchord_outer_solve_sequence(
      collected$tasks[sequence == id], registry,
      units_per_inch, backbone_outer
    )
    registry <- ggchord_outer_apply_solution(registry, solution)
  }
  ggchord_measure_external_leader_crossings(registry)
}
