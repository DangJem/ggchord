# Coord-owned perimeter annotation coordination for circular maps.

ggchord_external_candidates <- function(max_band = 16L, max_slot = 20L,
                                         class = "restriction") {
  candidates <- expand.grid(
    band = seq_len(max_band),
    slot = c(0L, as.vector(rbind(seq_len(max_slot), -seq_len(max_slot))))
  )
  # A shared dominant band is the first objective. Small local tangential
  # movements are cheaper than opening a new band, but remain bounded so an
  # annotation cannot migrate around the circle.
  radial_cost <- .9
  candidates$score <- (candidates$band - 1L) * radial_cost +
    abs(candidates$slot) * 1.60
  candidates[order(candidates$score, candidates$band,
    abs(candidates$slot), candidates$slot < 0), , drop = FALSE]
}

ggchord_open_circular_order <- function(angle, tie = seq_along(angle)) {
  if (length(angle) <= 1L) return(seq_along(angle))
  normalized <- (angle + 2 * pi) %% (2 * pi)
  sorted <- order(normalized, tie)
  gaps <- diff(c(normalized[sorted], normalized[sorted[1L]] + 2 * pi))
  cut <- which.max(gaps)
  c(if (cut < length(sorted)) sorted[seq.int(cut + 1L, length(sorted))]
    else integer(), sorted[seq_len(cut)])
}

ggchord_external_place <- function(row, occupied, units_per_inch,
                                    base_edge_radius = NULL,
                                    class = "restriction",
                                    clustered = FALSE,
                                    leader_segments = data.frame(),
                                    occupied_leaders = data.frame(),
                                    preferred_centres = data.frame(),
                                    rail_leaders = data.frame()) {
  x <- row$text_x %||% row$x
  y <- row$text_y %||% row$y
  radius <- sqrt(x^2 + y^2)
  if (!is.finite(radius) || radius <= 1e-8) {
    return(list(row = row, box = data.frame(), band = NA_integer_,
      slot = NA_integer_, collision = FALSE))
  }
  anchor_vector <- if (nrow(leader_segments) && all(is.finite(
      unlist(leader_segments[1L, c("x1", "y1")]))) ) {
    c(leader_segments$x1[1L], leader_segments$y1[1L])
  } else c(x, y)
  anchor_radius <- sqrt(sum(anchor_vector^2))
  direction <- if (is.finite(anchor_radius) && anchor_radius > 1e-8) {
    anchor_vector / anchor_radius
  } else c(x, y) / radius
  tangent <- c(-direction[2L], direction[1L])
  row$text <- if (identical(class, "primer") &&
      "feature_label" %in% names(row) &&
      !is.na(row$feature_label[1L]) && nzchar(row$feature_label[1L])) {
    row$feature_label
  } else row$label %||% row$text
  row$text_x <- x
  row$text_y <- y
  row$text_angle <- row$angle %||% row$text_angle %||% 0
  measured <- ggchord_text_boxes(
    row, units_per_inch = units_per_inch, box_padding = .045
  )
  occupied_leader_matrix <- if (nrow(occupied_leaders)) {
    as.matrix(occupied_leaders[, c("x1", "y1", "x2", "y2"), drop = FALSE])
  } else matrix(numeric(), ncol = 4L)
  radial_half <- abs(direction[1L]) * measured$w[1L] / 2 +
    abs(direction[2L]) * measured$h[1L] / 2
  if (is.null(base_edge_radius) || !is.finite(base_edge_radius)) {
    base_edge_radius <- max(0, radius - radial_half)
  }
  step <- max(.052, measured$h[1L] * 1.02)
  slot_step <- max(.038, measured$h[1L] * 1.05)
  assess_candidate <- function(candidate, band, slot,
                               segments = leader_segments,
                               shift_segments = TRUE) {
    candidate$x <- candidate$text_x
    candidate$y <- candidate$text_y
    # Candidate movement changes only the anchor. Reuse the exact measured
    # glyph dimensions instead of opening a graphics device for every band
    # and slot in dense mixed-annotation maps.
    box <- measured
    box_dx <- candidate$text_x - measured$x
    box_dy <- candidate$text_y - measured$y
    for (column in c("x", "cx", "xmin", "xmax")) {
      box[[column]] <- box[[column]] + box_dx
    }
    for (column in c("y", "cy", "ymin", "ymax")) {
      box[[column]] <- box[[column]] + box_dy
    }
    box_collision <- nrow(occupied) &&
      any(ggchord_oriented_box_overlaps(box, occupied))
    backbone_collision <- FALSE
    if (!is.null(base_edge_radius) && is.finite(base_edge_radius)) {
      nearest_x <- if (box$xmin[1L] > 0) box$xmin[1L] else if (
        box$xmax[1L] < 0) box$xmax[1L] else 0
      nearest_y <- if (box$ymin[1L] > 0) box$ymin[1L] else if (
        box$ymax[1L] < 0) box$ymax[1L] else 0
      backbone_collision <- sqrt(nearest_x^2 + nearest_y^2) <
        base_edge_radius
      box_collision <- box_collision || backbone_collision
    }
    candidate_leaders <- segments
    leader_collision <- FALSE
    if (nrow(candidate_leaders)) {
      if (isTRUE(shift_segments)) {
        dx <- candidate$text_x - x
        dy <- candidate$text_y - y
        moving <- if ("move_endpoint" %in% names(candidate_leaders)) {
          candidate_leaders$move_endpoint %in% TRUE
        } else rep(TRUE, nrow(candidate_leaders))
        candidate_leaders$x2[moving] <- candidate_leaders$x2[moving] + dx
        candidate_leaders$y2[moving] <- candidate_leaders$y2[moving] + dy
      }
      if (nrow(occupied_leaders)) {
        leader_collision <- ggchord_external_route_crosses(
          candidate_leaders, occupied_leader_matrix
        )
      }
      if (!box_collision && leader_collision && nrow(candidate_leaders) &&
          nrow(occupied_leaders)) {
        # A label can have a collision-free bbox while its straight leader
        # cuts across a dense restriction fan. Search a single angular elbow
        # outside the backbone before rejecting that circular-track slot.
        # This keeps the route visibly straight/折线 and never approximates a
        # circular arc with many sampled points.
        root <- c(candidate_leaders$x1[1L], candidate_leaders$y1[1L])
        endpoint <- c(tail(candidate_leaders$x2, 1L),
          tail(candidate_leaders$y2, 1L))
        root_radius <- sqrt(sum(root^2))
        endpoint_radius <- sqrt(sum(endpoint^2))
        root_angle <- atan2(root[2L], root[1L])
        endpoint_angle <- atan2(endpoint[2L], endpoint[1L])
        delta_angle <- atan2(sin(endpoint_angle - root_angle),
          cos(endpoint_angle - root_angle))
        segment_radius <- function(p0, p1) {
          delta <- p1 - p0
          denominator <- sum(delta^2)
          parameter <- if (denominator <= 1e-16) 0 else
            max(0, min(1, -sum(p0 * delta) / denominator))
          sqrt(sum((p0 + parameter * delta)^2))
        }
        route_crosses <- function(route) {
          ggchord_external_route_crosses(route, occupied_leader_matrix)
        }
        safe_radius <- max(0, min(root_radius,
          base_edge_radius %||% root_radius) - .004)
        elbow_angles <- root_angle + delta_angle *
          c(1, .875, .75, .625, .5, .375, .25, .125, 0)
        elbow_radii <- max(root_radius, endpoint_radius,
          base_edge_radius %||% 0) + c(0, .03, .065, .11, .17, .25, .35, .48)
        routed <- NULL
        for (elbow_radius in elbow_radii) {
          for (elbow_angle in elbow_angles) {
            elbow <- elbow_radius * c(cos(elbow_angle), sin(elbow_angle))
            if (segment_radius(root, elbow) < safe_radius ||
                segment_radius(elbow, endpoint) < safe_radius) next
            route <- data.frame(
              x1 = c(root[1L], elbow[1L]),
              y1 = c(root[2L], elbow[2L]),
              x2 = c(elbow[1L], endpoint[1L]),
              y2 = c(elbow[2L], endpoint[2L]),
              move_endpoint = c(FALSE, FALSE)
            )
            if (!route_crosses(route)) {
              routed <- route
              break
            }
          }
          if (!is.null(routed)) break
        }
        if (!is.null(routed)) {
          candidate_leaders <- routed
          leader_collision <- FALSE
        }
      }
    }
    collision <- box_collision || leader_collision
    list(row = candidate, box = box, band = band, slot = slot,
      collision = collision, box_collision = box_collision,
      backbone_collision = backbone_collision,
      leader_collision = leader_collision,
      leader_segments = candidate_leaders, radial_half = radial_half)
  }
  chosen <- NULL
  # Cartesian perimeter rails are intentionally not candidates.  The retained
  # argument is schema-compatible with earlier development builds only.
  if (FALSE && nrow(leader_segments) && nrow(rail_leaders)) {
    root <- c(leader_segments$x1[1L], leader_segments$y1[1L])
    root_side <- if (abs(root[1L]) >= abs(root[2L])) {
      if (root[1L] >= 0) "right" else "left"
    } else if (root[2L] >= 0) "top" else "bottom"
    horizontal <- root_side %in% c("top", "bottom")
    fixed_endpoint <- if (horizontal) rail_leaders$y2 else rail_leaders$x2
    endpoint_key <- round(fixed_endpoint, 3L)
    counts <- table(endpoint_key)
    eligible_keys <- as.numeric(names(counts)[counts >= 3L])
    eligible_keys <- eligible_keys[is.finite(eligible_keys) &
      if (root_side %in% c("right", "top")) {
        eligible_keys > (base_edge_radius %||% 0)
      } else eligible_keys < -(base_edge_radius %||% 0)]
    if (length(eligible_keys)) {
      rail_value <- eligible_keys[which.max(abs(eligible_keys))]
      rail_rows <- which(abs(fixed_endpoint - rail_value) <= .0015)
      root_axis <- if (horizontal) rail_leaders$x1[rail_rows] else
        rail_leaders$y1[rail_rows]
      root_fixed <- if (horizontal) rail_leaders$y1[rail_rows] else
        rail_leaders$x1[rail_rows]
      endpoint_axis <- if (horizontal) rail_leaders$x2[rail_rows] else
        rail_leaders$y2[rail_rows]
      ord <- order(root_axis, endpoint_axis)
      root_axis <- root_axis[ord]
      root_fixed <- root_fixed[ord]
      endpoint_axis <- endpoint_axis[ord]
      source_axis <- if (horizontal) root[1L] else root[2L]
      target_axis <- stats::approx(
        root_axis, endpoint_axis,
        xout = source_axis,
        rule = 2, ties = "ordered"
      )$y
      entry_fixed <- stats::approx(
        root_axis, root_fixed, xout = source_axis,
        rule = 2, ties = "ordered"
      )$y
      moving <- if ("move_endpoint" %in% names(leader_segments)) {
        which(leader_segments$move_endpoint %in% TRUE)[1L]
      } else 1L
      old_endpoint <- c(leader_segments$x2[moving],
        leader_segments$y2[moving])
      target <- if (horizontal) c(target_axis, rail_value) else
        c(rail_value, target_axis)
      entry <- if (horizontal) c(source_axis, entry_fixed) else
        c(entry_fixed, source_axis)
      base_centre <- c(x, y) + target - old_endpoint
      outward <- switch(root_side, right = c(1, 0), left = c(-1, 0),
        top = c(0, 1), bottom = c(0, -1))
      for (band in 0:8) {
        candidate <- row
        centre <- base_centre + outward * band * step
        candidate$text_x <- centre[1L]
        candidate$text_y <- centre[2L]
        label_endpoint <- old_endpoint + centre - c(x, y)
        routed <- data.frame(
          x1 = c(root[1L], entry[1L], target[1L]),
          y1 = c(root[2L], entry[2L], target[2L]),
          x2 = c(entry[1L], target[1L], label_endpoint[1L]),
          y2 = c(entry[2L], target[2L], label_endpoint[2L]),
          move_endpoint = c(FALSE, FALSE, FALSE)
        )
        trial <- assess_candidate(
          candidate, band + 1L, 0L, segments = routed,
          shift_segments = FALSE
        )
        if (!trial$collision) {
          chosen <- trial
          break
        }
      }
    }
  }
  if ((is.null(chosen) || isTRUE(chosen$collision)) &&
      nrow(preferred_centres)) for (peer in seq_len(nrow(preferred_centres))) {
    peer_radius <- sqrt(preferred_centres$x[peer]^2 +
      preferred_centres$y[peer]^2)
    if (!is.finite(peer_radius) || peer_radius <= 1e-8) next
    peer_direction <- c(preferred_centres$x[peer],
      preferred_centres$y[peer]) / peer_radius
    peer_tangent <- c(-peer_direction[2L], peer_direction[1L])
    bands <- if (isTRUE(clustered)) 0:7 else seq_len(8L)
    slots <- if (isTRUE(clustered))
      c(0L, -1L, 1L, -2L, 2L, -3L, 3L, -4L, 4L) else
      c(0L, -1L, 1L, -2L, 2L)
    for (band in bands) for (slot in slots) {
      candidate <- row
      candidate_radius <- peer_radius + band * step
      candidate$text_x <- peer_direction[1L] * candidate_radius +
        peer_tangent[1L] * slot * slot_step
      candidate$text_y <- peer_direction[2L] * candidate_radius +
        peer_tangent[2L] * slot * slot_step
      trial <- assess_candidate(candidate, band + 1L, slot)
      if (!trial$collision) {
        chosen <- trial
        break
      }
    }
    if (!is.null(chosen) && !chosen$collision) break
  }
  candidates <- ggchord_external_candidates(class = class)
  best_failed <- NULL
  best_failed_score <- Inf
  if (is.null(chosen) || isTRUE(chosen$collision)) for (i in seq_len(nrow(candidates))) {
    band <- candidates$band[i]
    slot <- candidates$slot[i]
    candidate_radius <- base_edge_radius + radial_half + (band - 1L) * step
    candidate <- row
    candidate$text_x <- direction[1L] * candidate_radius +
      tangent[1L] * slot * slot_step
    candidate$text_y <- direction[2L] * candidate_radius +
      tangent[2L] * slot * slot_step
    trial <- assess_candidate(candidate, band, slot)
    if (!trial$collision) {
      chosen <- trial
      break
    }
    # A dense mixed fan may have no completely free slot yet. Preserve the
    # least harmful visible candidate instead of accidentally using the last
    # band/slot in the search order.
    crossing_count <- if (nrow(trial$leader_segments) &&
        nrow(occupied_leaders)) {
      ggchord_external_route_crossing_count(
        trial$leader_segments, occupied_leader_matrix
      )
    } else 0L
    overlap_count <- if (nrow(occupied)) {
      sum(ggchord_oriented_box_overlaps(trial$box, occupied))
    } else 0L
    failed_score <- as.numeric(trial$backbone_collision) * 1e6 +
      as.numeric(trial$box_collision) * 10000 +
      overlap_count * 1000 + crossing_count * 100 +
      (band - 1L) * .9 + abs(slot) * 1.6
    if (failed_score < best_failed_score) {
      best_failed <- trial
      best_failed_score <- failed_score
    }
  }
  if (!is.null(best_failed) &&
      (is.null(chosen) || isTRUE(chosen$collision))) chosen <- best_failed
  # Never fall back to a Cartesian side wall.  Horizontal text is allowed to
  # extend beyond the circular anchor track; only its measured box participates
  # in collision and local outward stacking.
  if (FALSE && isTRUE(chosen$collision) && nrow(occupied)) {
    box_cx <- (occupied$xmin + occupied$xmax) / 2
    box_cy <- (occupied$ymin + occupied$ymax) / 2
    side_specs <- list()
    right <- which(box_cx > .45)
    left <- which(box_cx < -.45)
    top <- which(box_cy > .45)
    bottom <- which(box_cy < -.45)
    if (length(right)) side_specs$right <- c(min(occupied$xmin[right]), 0)
    if (length(left)) side_specs$left <- c(max(occupied$xmax[left]), 0)
    if (length(top)) side_specs$top <- c(0, min(occupied$ymin[top]))
    if (length(bottom)) side_specs$bottom <- c(0, max(occupied$ymax[bottom]))
    side_vectors <- list(right = c(1, 0), left = c(-1, 0),
      top = c(0, 1), bottom = c(0, -1))
    side_names <- names(side_specs)
    if (length(side_names)) {
      side_names <- side_names[order(vapply(side_names, function(side) {
        -sum(direction * side_vectors[[side]])
      }, numeric(1L)))]
      rail_candidates <- ggchord_external_candidates(
        max_band = 4L, max_slot = 12L, class = class
      )
      for (side in side_names) {
        for (i in seq_len(nrow(rail_candidates))) {
          band <- rail_candidates$band[i]
          slot <- rail_candidates$slot[i]
          candidate <- row
          if (side == "right") {
            candidate$text_x <- side_specs[[side]][1L] +
              measured$w[1L] / 2 + .050 + (band - 1L) * step
            candidate$text_y <- y + slot * slot_step
          } else if (side == "left") {
            candidate$text_x <- side_specs[[side]][1L] -
              measured$w[1L] / 2 - .050 - (band - 1L) * step
            candidate$text_y <- y + slot * slot_step
          } else if (side == "top") {
            candidate$text_x <- x + slot * slot_step
            candidate$text_y <- side_specs[[side]][2L] +
              measured$h[1L] / 2 + .050 + (band - 1L) * step
          } else {
            candidate$text_x <- x + slot * slot_step
            candidate$text_y <- side_specs[[side]][2L] -
              measured$h[1L] / 2 - .050 - (band - 1L) * step
          }
          trial <- assess_candidate(candidate, band, slot)
          if (!trial$collision) {
            chosen <- trial
            break
          }
        }
        if (!chosen$collision) break
      }
    }
  }
  chosen
}

ggchord_shift_feature_callout <- function(geometry, row, candidate) {
  dx <- candidate$row$text_x - geometry$text_x[row]
  dy <- candidate$row$text_y - geometry$text_y[row]
  geometry$text_x[row] <- candidate$row$text_x
  geometry$text_y[row] <- candidate$row$text_y
  geometry$x[row] <- candidate$row$text_x
  geometry$y[row] <- candidate$row$text_y
  geometry$outer_annotation_region[row] <- if (
    "annotation_class" %in% names(geometry) &&
      identical(as.character(geometry$annotation_class[row]), "primer")) {
    "primer"
  } else "feature"
  geometry$outer_track[row] <- candidate$band
  geometry$outer_slot[row] <- candidate$slot
  geometry$dominant_outer_band[row] <- 1L
  geometry$outer_spill_reason[row] <- if (isTRUE(candidate$collision)) {
    if (isTRUE(candidate$box_collision)) {
      "bbox_collision_unresolved"
    } else if (isTRUE(candidate$leader_collision)) {
      "leader_crossing_unresolved"
    } else "unresolved"
  } else if (candidate$band > 1L) {
    "bbox_collision"
  } else NA_character_
  geometry$leader_crossing_count[row] <- 0L
  segment_rows <- which(geometry$.component %in% "segment" &
    geometry$source_row == geometry$source_row[row])
  if (length(segment_rows)) {
    geometry$xend[segment_rows] <- geometry$x[segment_rows]
    geometry$yend[segment_rows] <- geometry$y[segment_rows]
    routed <- candidate$leader_segments
    if (nrow(routed) > length(segment_rows)) {
      template <- geometry[tail(segment_rows, 1L), , drop = FALSE]
      extra <- template[rep(1L, nrow(routed) - length(segment_rows)), ,
        drop = FALSE]
      geometry <- rbind(geometry, extra)
      segment_rows <- c(segment_rows,
        seq.int(nrow(geometry) - nrow(extra) + 1L, nrow(geometry)))
    }
    use <- seq_len(min(length(segment_rows), nrow(routed)))
    if (length(use)) {
      target_rows <- segment_rows[use]
      geometry$x[target_rows] <- routed$x1[use]
      geometry$y[target_rows] <- routed$y1[use]
      geometry$xend[target_rows] <- routed$x2[use]
      geometry$yend[target_rows] <- routed$y2[use]
    }
    if ("x1" %in% names(geometry)) {
      geometry$x1[segment_rows] <- geometry$xend[segment_rows]
      geometry$y1[segment_rows] <- geometry$yend[segment_rows]
    }
    geometry$outer_annotation_region[segment_rows] <-
      paste0(geometry$outer_annotation_region[row], "_leader")
    geometry$outer_track[segment_rows] <- candidate$band
    geometry$outer_slot[segment_rows] <- candidate$slot
    geometry$dominant_outer_band[segment_rows] <- 1L
    geometry$outer_spill_reason[segment_rows] <-
      geometry$outer_spill_reason[row]
    geometry$leader_crossing_count[segment_rows] <- 0L
  }
  geometry
}

ggchord_external_route_crosses <- function(route, occupied,
                                            tolerance = 1e-9) {
  if (!nrow(route) || !nrow(occupied)) return(FALSE)
  ax <- occupied[, 1L]
  ay <- occupied[, 2L]
  bx <- occupied[, 3L] - ax
  by <- occupied[, 4L] - ay
  for (i in seq_len(nrow(route))) {
    px <- route$x1[i]
    py <- route$y1[i]
    rx <- route$x2[i] - px
    ry <- route$y2[i] - py
    denominator <- rx * by - ry * bx
    valid <- is.finite(denominator) & abs(denominator) > tolerance
    if (!any(valid)) next
    qx <- ax[valid] - px
    qy <- ay[valid] - py
    denom <- denominator[valid]
    t <- (qx * by[valid] - qy * bx[valid]) / denom
    u <- (qx * ry - qy * rx) / denom
    if (any(t > tolerance & t < 1 - tolerance &
            u > tolerance & u < 1 - tolerance)) return(TRUE)
  }
  FALSE
}

ggchord_external_route_crossing_count <- function(route, occupied,
                                                   tolerance = 1e-9) {
  if (!nrow(route) || !nrow(occupied)) return(0L)
  ax <- occupied[, 1L]
  ay <- occupied[, 2L]
  bx <- occupied[, 3L] - ax
  by <- occupied[, 4L] - ay
  crossed <- rep(FALSE, nrow(occupied))
  for (i in seq_len(nrow(route))) {
    px <- route$x1[i]
    py <- route$y1[i]
    rx <- route$x2[i] - px
    ry <- route$y2[i] - py
    denominator <- rx * by - ry * bx
    valid <- is.finite(denominator) & abs(denominator) > tolerance
    if (!any(valid)) next
    qx <- ax[valid] - px
    qy <- ay[valid] - py
    denom <- denominator[valid]
    t <- (qx * by[valid] - qy * bx[valid]) / denom
    u <- (qx * ry - qy * rx) / denom
    crossed[which(valid)] <- crossed[which(valid)] |
      (t > tolerance & t < 1 - tolerance &
        u > tolerance & u < 1 - tolerance)
  }
  sum(crossed)
}

ggchord_external_leader_crosses <- function(a, b, tolerance = 1e-9) {
  cross2 <- function(x, y) x[1L] * y[2L] - x[2L] * y[1L]
  p <- c(a$x1, a$y1)
  r <- c(a$x2 - a$x1, a$y2 - a$y1)
  q <- c(b$x1, b$y1)
  s <- c(b$x2 - b$x1, b$y2 - b$y1)
  denominator <- cross2(r, s)
  if (!is.finite(denominator) || abs(denominator) <= tolerance) return(FALSE)
  t <- cross2(q - p, s) / denominator
  u <- cross2(q - p, r) / denominator
  # Shared roots and label attachments are intentional contacts, not crossings.
  t > tolerance && t < 1 - tolerance &&
    u > tolerance && u < 1 - tolerance
}

ggchord_measure_external_leader_crossings <- function(registry) {
  leaders <- list()
  add_leader <- function(id, component, source, segments) {
    finite <- stats::complete.cases(segments[, c("x1", "y1", "x2", "y2")])
    segments <- segments[finite, , drop = FALSE]
    if (nrow(segments)) leaders[[length(leaders) + 1L]] <<- list(
      id = id, component = component, source = source,
      segments = segments, crossings = 0L
    )
  }
  for (id in names(registry)) {
    geometry <- registry[[id]]$gene_label_repel
    if (is.data.frame(geometry) && nrow(geometry)) {
      external_leader <- if ("outer_annotation_region" %in% names(geometry)) {
        !is.na(geometry$outer_annotation_region) &
          geometry$outer_annotation_region %in%
            c("feature_leader", "primer_leader")
      } else rep(FALSE, nrow(geometry))
      rows <- which(geometry$.component %in% "segment" &
        is.finite(geometry$source_row) & external_leader)
      for (source in unique(geometry$source_row[rows])) {
        index <- rows[geometry$source_row[rows] == source]
        add_leader(id, "gene_label_repel", source, data.frame(
          x1 = geometry$x[index], y1 = geometry$y[index],
          x2 = geometry$xend[index], y2 = geometry$yend[index]
        ))
      }
    }
    sites <- registry[[id]]$restriction_site
    if (is.data.frame(sites) && nrow(sites)) {
      rows <- which(sites$restriction_component %in% "leader" &
        is.finite(sites$source_row))
      for (source in unique(sites$source_row[rows])) {
        index <- rows[sites$source_row[rows] == source]
        if (length(index) < 2L) next
        # `append_path()` stores every polyline in draw order.  Radius is not
        # monotone for four-sided rails (an exterior shoulder may be farther
        # out than the final text-edge attachment), so sorting by radius
        # invents segments that are never rendered and reports false
        # crossings.  Measure exactly the consecutive path rows instead.
        add_leader(id, "restriction_site", source, data.frame(
          x1 = sites$x[index[-length(index)]],
          y1 = sites$y[index[-length(index)]],
          x2 = sites$x[index[-1L]], y2 = sites$y[index[-1L]]
        ))
      }
    }
  }
  if (length(leaders) > 1L) {
    pairs <- utils::combn(seq_along(leaders), 2L)
    for (column in seq_len(ncol(pairs))) {
      i <- pairs[1L, column]
      j <- pairs[2L, column]
      crosses <- any(vapply(seq_len(nrow(leaders[[i]]$segments)), function(a) {
        any(vapply(seq_len(nrow(leaders[[j]]$segments)), function(b) {
          ggchord_external_leader_crosses(
            leaders[[i]]$segments[a, ], leaders[[j]]$segments[b, ]
          )
        }, logical(1L)))
      }, logical(1L)))
      if (crosses) {
        leaders[[i]]$crossings <- leaders[[i]]$crossings + 1L
        leaders[[j]]$crossings <- leaders[[j]]$crossings + 1L
      }
    }
  }
  for (leader in leaders) {
    geometry <- registry[[leader$id]][[leader$component]]
    rows <- is.finite(geometry$source_row) &
      geometry$source_row == leader$source
    geometry$leader_crossing_count[rows] <- leader$crossings
    registry[[leader$id]][[leader$component]] <- geometry
  }
  registry
}

# Resolve every exterior annotation against one measured occupancy table.
# Restriction sites first establish circular anchor tracks. Feature and primer
# callouts then use the same measured occupancy table and the same circular
# band semantics. Horizontal text may leave the nominal annulus; its complete
# box is a collision body, not a requirement to join a Cartesian side rail.
ggchord_share_external_annotations <- function(registry, layout) {
  if (!isTRUE(layout$circular) || !length(registry)) return(registry)
  units_per_inch <- layout$text_units_per_inch %||% .30
  backbone_outer <- if (is.data.frame(layout$backbone_bounds) &&
      nrow(layout$backbone_bounds)) {
    max(layout$backbone_bounds$outer, na.rm = TRUE) + .040
  } else NULL
  occupied <- data.frame()
  occupied_leaders <- data.frame()
  restriction_leaders <- data.frame()
  placed_external <- data.frame(
    root_x = numeric(), root_y = numeric(), x = numeric(), y = numeric(),
    cluster_key = character()
  )

  # Register the already solved restriction contours first. Their exact
  # rendered boxes (including justification) become obstacles for both other
  # annotation classes.
  for (id in names(registry)) {
    sites <- registry[[id]]$restriction_site
    if (!is.data.frame(sites) || !nrow(sites)) next
    rows <- which(sites$.component %in% "label" &
      !is.na(sites$label) & nzchar(sites$label))
    if (!length(rows)) next
    for (nm in c("outer_annotation_region", "outer_spill_reason")) {
      if (!nm %in% names(sites)) sites[[nm]] <- NA_character_
    }
    for (nm in c("outer_track", "outer_slot", "dominant_outer_band",
        "leader_crossing_count")) {
      if (!nm %in% names(sites)) sites[[nm]] <- NA_integer_
    }
    sites$outer_annotation_region[rows] <- "restriction"
    sites$dominant_outer_band[rows] <- 1L
    sites$leader_crossing_count[rows] <- 0L
    used_bands <- sort(unique(sites$outer_track[rows]))
    used_bands <- used_bands[is.finite(used_bands)]
    if (length(used_bands)) {
      remap <- stats::setNames(seq_along(used_bands), used_bands)
      affected <- is.finite(sites$outer_track)
      sites$outer_track[affected] <- unname(remap[
        as.character(sites$outer_track[affected])
      ])
    }
    measured <- sites[rows, , drop = FALSE]
    measured$text <- measured$label
    measured$text_x <- measured$x
    measured$text_y <- measured$y
    measured$text_angle <- measured$angle %||% 0
    boxes <- ggchord_text_boxes(
      measured, units_per_inch = units_per_inch, box_padding = .045
    )
    occupied <- if (!nrow(occupied)) boxes else
      ggchord_rbind_fill(list(occupied, boxes))
    leader_rows <- which(sites$restriction_component %in% "leader" &
      is.finite(sites$source_row))
    for (source in unique(sites$source_row[leader_rows])) {
      index <- leader_rows[sites$source_row[leader_rows] == source]
      if (length(index) < 2L) next
      segments <- data.frame(
        x1 = sites$x[index[-length(index)]],
        y1 = sites$y[index[-length(index)]],
        x2 = sites$x[index[-1L]], y2 = sites$y[index[-1L]]
      )
      occupied_leaders <- if (!nrow(occupied_leaders)) segments else
        rbind(occupied_leaders, segments)
      restriction_leaders <- if (!nrow(restriction_leaders)) segments else
        rbind(restriction_leaders, segments)
      label_row <- rows[match(source, sites$source_row[rows])]
      if (length(label_row) == 1L && is.finite(label_row)) {
        placed_external <- rbind(placed_external, data.frame(
          root_x = segments$x1[1L], root_y = segments$y1[1L],
          x = sites$x[label_row], y = sites$y[label_row],
          cluster_key = NA_character_
        ))
      }
    }
    registry[[id]]$restriction_site <- sites
  }

  # Collect all external feature and primer callouts before placing any of
  # them.  Layer order is not genomic order: a primer layer inserted after a
  # feature layer can otherwise claim a slot behind a later root and force
  # their leaders to cross.  The shared pass uses each real leader root.
  tasks <- list()
  for (id in names(registry)) {
    geometry <- registry[[id]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry)) next
    external <- if ("feature_label_mode" %in% names(geometry)) {
      geometry$feature_label_mode %in% "external"
    } else rep(FALSE, nrow(geometry))
    if ("annotation_class" %in% names(geometry)) {
      external <- external | geometry$annotation_class %in% "primer"
    }
    rows <- which(geometry$.component %in% "text" & external &
      !is.na(geometry$label) & nzchar(geometry$label))
    if (!length(rows)) next
    for (nm in c("outer_annotation_region", "outer_spill_reason")) {
      if (!nm %in% names(geometry)) geometry[[nm]] <- NA_character_
    }
    for (nm in c("outer_track", "outer_slot", "dominant_outer_band",
        "leader_crossing_count")) {
      if (!nm %in% names(geometry)) geometry[[nm]] <- NA_integer_
    }
    for (row in rows) {
      segment_rows <- which(geometry$.component %in% "segment" &
        geometry$source_row == geometry$source_row[row])
      if (length(segment_rows)) {
        root <- segment_rows[which.max(
          (geometry$x[segment_rows] - geometry$text_x[row])^2 +
            (geometry$y[segment_rows] - geometry$text_y[row])^2
        )]
        root_xy <- c(geometry$x[root], geometry$y[root])
      } else {
        root_xy <- if (all(c("anchor_x", "anchor_y") %in%
            names(geometry))) c(geometry$anchor_x[row],
          geometry$anchor_y[row]) else numeric()
      }
      if (length(root_xy) != 2L || any(!is.finite(root_xy))) root_xy <-
        c(geometry$text_x[row], geometry$text_y[row])
      position <- if ("anchor_position" %in% names(geometry))
        geometry$anchor_position[row] else row
      if (length(position) != 1L || !is.finite(position)) position <- row
      tasks[[length(tasks) + 1L]] <- data.frame(
        id = id, row = row,
        accver = as.character(geometry$accver[row]),
        angle = atan2(root_xy[2L], root_xy[1L]),
        anchor_position = position,
        stringsAsFactors = FALSE
      )
    }
    registry[[id]]$gene_label_repel <- geometry
  }
  if (length(tasks)) {
    tasks <- do.call(rbind, tasks)
    ordered <- integer()
    for (sid in unique(tasks$accver)) {
      rows <- which(tasks$accver == sid)
      ordered <- c(ordered, rows[ggchord_open_circular_order(
        tasks$angle[rows], tasks$anchor_position[rows]
      )])
    }
    for (task_index in ordered) {
      id <- tasks$id[task_index]
      row <- tasks$row[task_index]
      geometry <- registry[[id]]$gene_label_repel
      class <- if ("annotation_class" %in% names(geometry) &&
        identical(as.character(geometry$annotation_class[row]), "primer")) {
        "primer"
      } else "feature"
      segment_rows <- which(geometry$.component %in% "segment" &
        geometry$source_row == geometry$source_row[row])
      leader_segments <- data.frame()
      if (length(segment_rows)) {
        terminal <- segment_rows[which.min(
          (geometry$xend[segment_rows] - geometry$text_x[row])^2 +
            (geometry$yend[segment_rows] - geometry$text_y[row])^2
        )]
        root <- segment_rows[which.max(
          (geometry$x[segment_rows] - geometry$text_x[row])^2 +
            (geometry$y[segment_rows] - geometry$text_y[row])^2
        )]
        leader_segments <- data.frame(
          x1 = geometry$x[root], y1 = geometry$y[root],
          x2 = geometry$xend[terminal], y2 = geometry$yend[terminal],
          move_endpoint = TRUE
        )
      }
      preferred_centres <- data.frame()
      cluster_key <- if ("outer_cluster_key" %in% names(geometry))
        as.character(geometry$outer_cluster_key[row]) else NA_character_
      clustered <- !is.na(cluster_key) && nzchar(cluster_key)
      if (clustered && nrow(placed_external)) {
        nearby <- which(!is.na(placed_external$cluster_key) &
          placed_external$cluster_key == cluster_key)
        if (length(nearby)) preferred_centres <- placed_external[nearby,
          c("x", "y"), drop = FALSE]
      } else if (nrow(leader_segments) && nrow(placed_external) &&
          !clustered) {
        root_distance <- sqrt(
          (placed_external$root_x - leader_segments$x1[1L])^2 +
            (placed_external$root_y - leader_segments$y1[1L])^2
        )
        nearby <- which(root_distance <= .080)
        if (length(nearby)) preferred_centres <- placed_external[nearby,
          c("x", "y"), drop = FALSE]
      }
      candidate <- ggchord_external_place(
        geometry[row, , drop = FALSE], occupied, units_per_inch,
        base_edge_radius = backbone_outer,
        class = class, clustered = clustered,
        leader_segments = leader_segments,
        occupied_leaders = occupied_leaders,
        preferred_centres = preferred_centres,
        rail_leaders = restriction_leaders
      )
      geometry <- ggchord_shift_feature_callout(geometry, row, candidate)
      occupied <- if (!nrow(occupied)) candidate$box else
        ggchord_rbind_fill(list(occupied, candidate$box))
      if (nrow(candidate$leader_segments)) {
        placed_leaders <- candidate$leader_segments[
          c("x1", "y1", "x2", "y2")
        ]
        occupied_leaders <- if (!nrow(occupied_leaders)) {
          placed_leaders
        } else rbind(occupied_leaders, placed_leaders)
        placed_external <- rbind(placed_external, data.frame(
          root_x = candidate$leader_segments$x1[1L],
          root_y = candidate$leader_segments$y1[1L],
          x = candidate$row$text_x, y = candidate$row$text_y,
          cluster_key = cluster_key
        ))
      }
      registry[[id]]$gene_label_repel <- geometry
    }
  }
  for (id in names(registry)) {
    geometry <- registry[[id]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry) ||
        !"outer_track" %in% names(geometry)) next
    used_bands <- sort(unique(geometry$outer_track))
    used_bands <- used_bands[is.finite(used_bands)]
    if (length(used_bands)) {
      remap <- stats::setNames(seq_along(used_bands), used_bands)
      affected <- is.finite(geometry$outer_track)
      geometry$outer_track[affected] <- unname(remap[
        as.character(geometry$outer_track[affected])
      ])
    }
    registry[[id]]$gene_label_repel <- geometry
  }
  ggchord_measure_external_leader_crossings(registry)
}
