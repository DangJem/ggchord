# Horizontal external labels on offsets of the actual sequence curves.
# Parameters are arc lengths, not angles around the global origin, so the
# same packing works for unequal radii, rotated curves and open sequences.
ggchord_radial_label_lanes <- function(gl, seq_arcs, side = "outside",
                                        units_per_inch = 0.35,
                                        box_padding = 0.04,
                                        point_padding = 0.05,
                                        repel_boxes = NULL, .attempt = 0L, .fixed_paths = NULL,
                                        .minimum_levels = numeric(), .balance = TRUE) {
  n <- nrow(gl)
  active <- which(!is.na(gl$text) & nzchar(gl$text))
  original <- gl
  source <- ggchord_label_curve_frame(gl, seq_arcs)
  signs <- switch(side, outside = rep(1, n), inside = rep(-1, n),
                  ifelse(source$signed_distance < 0, -1, 1))
  gl$text_angle <- 0
  gl$hjust <- 0.5
  gl$vjust <- 0.5
  boxes <- ggchord_text_boxes(gl, units_per_inch = units_per_inch,
                             box_padding = box_padding)
  gl$.radial_bend_x <- gl$anchor_x + signs * source$outward_x *
    max(0.025, 0.08 * units_per_inch)
  gl$.radial_bend_y <- gl$anchor_y + signs * source$outward_y *
    max(0.025, 0.08 * units_per_inch)
  tracks <- rep(NA_integer_, n)
  directions <- rep(NA_character_, n)
  parameter <- rep(NA_real_, n)
  occupied <- ggchord_text_boxes(data.frame())
  paths <- data.frame(x0 = numeric(), y0 = numeric(), x1 = numeric(),
                      y1 = numeric(), group = integer())
  if (!is.null(.fixed_paths)) paths <- .fixed_paths
  groups <- split(active, paste(gl$seq_id[active], signs[active], sep = "\r"))
  if (length(groups) > 1L) {
    order_groups <- seq_along(groups)
    if (.attempt %% 2L) order_groups <- rev(order_groups)
    rotation <- (.attempt %/% 2L) %% length(groups)
    if (rotation) order_groups <- c(tail(order_groups, -rotation), utils::head(order_groups, rotation))
    groups <- groups[order_groups]
  }
  arc_ids <- vapply(seq_arcs, function(a) as.character(a$seq_id[1]), character(1))

  for (rows in groups) {
    arc <- seq_arcs[[match(gl$seq_id[rows[1]], arc_ids)]]
    s <- c(0, cumsum(sqrt(diff(arc$x)^2 + diff(arc$y)^2)))
    keep <- !duplicated(s)
    arc <- arc[keep, , drop = FALSE]
    s <- s[keep]
    if (length(s) < 2L) ggchord_stop("Sequence curve must have positive length")
    tx <- c(diff(arc$x), tail(diff(arc$x), 1))
    ty <- c(diff(arc$y), tail(diff(arc$y), 1))
    norm <- sqrt(tx^2 + ty^2)
    path_direction <- attr(arc, "ggchord_path_direction") %||% 1
    nx <- path_direction * ty / norm
    ny <- -path_direction * tx / norm
    nx <- nx * signs[rows[1]]
    ny <- ny * signs[rows[1]]
    preferred <- vapply(rows, function(i) {
      s[which.min((arc$x - gl$anchor_x[i])^2 + (arc$y - gl$anchor_y[i])^2)]
    }, numeric(1))
    names(preferred) <- rows
    rows <- rows[order(preferred, rows)]
    interpolate <- function(v, at) stats::approx(s, v, at, rule = 2)$y
    anchor_offset <- (gl$anchor_x[rows] - source$curve_x[rows]) *
      signs[rows] * source$outward_x[rows] +
      (gl$anchor_y[rows] - source$curve_y[rows]) * signs[rows] * source$outward_y[rows]
    clearance <- max(anchor_offset, 0) + point_padding + 0.18 * units_per_inch
    step <- max(0.015, 0.05 * units_per_inch)
    # Keep a small physical shoulder at each endpoint, rather than extending
    # the contour indefinitely into a neighbouring sequence's space.
    endpoint_padding <- 0.25 * units_per_inch
    shifts <- seq(0, max(diff(range(s)) / 2, sum(boxes$bw[rows]) / 2), by = step)
    offsets <- sort(unique(c(shifts, -shifts)))
    offsets <- offsets[order(abs(offsets), offsets)]
    best <- NULL
    group_key <- paste(gl$seq_id[rows[1]], signs[rows[1]], sep = "\r")
    minimum_level <- .minimum_levels[group_key]
    if (!length(minimum_level) || is.na(minimum_level)) minimum_level <- 0
    levels <- c(0:12, 16, 24, 40)
    for (level in levels[levels >= minimum_level]) {
      distance <- clearance + level * max(0.08, 0.25 * units_per_inch)
      for (reverse in c(FALSE, TRUE)) {
        for (spread in c(1, 1.25, 1.5, 2)) {
          for (start_shift in c(0, 0.3)) {
            ordered <- if (reverse) rev(rows) else rows
            candidate <- gl
            local_boxes <- occupied
            local_paths <- paths
            previous <- if (reverse) Inf else -Inf
            cost <- 0
            solved <- TRUE
            local_parameter <- parameter
            local_direction <- directions
            for (index in seq_along(ordered)) {
              i <- ordered[index]
              selected <- FALSE
              bias <- (preferred[as.character(i)] - mean(preferred)) * (spread - 1)
              if (index == 1L) bias <- bias + if (reverse) start_shift else -start_shift
              local_offsets <- offsets[order(abs(offsets - bias), offsets)]
              for (shift in local_offsets) {
                at <- preferred[as.character(i)] + shift
                if (at < -endpoint_padding || at > tail(s, 1) + endpoint_padding) next
                if ((!reverse && at <= previous + 1e-8) ||
                    (reverse && at >= previous - 1e-8)) next
                normal <- c(interpolate(nx, at), interpolate(ny, at))
                normal <- normal / sqrt(sum(normal^2))
                point <- c(interpolate(arc$x, at), interpolate(arc$y, at))
                if (at < 0) point <- point + at * c(tx[1], ty[1]) / norm[1]
                if (at > tail(s, 1)) point <- point + (at - tail(s, 1)) *
                  c(tail(tx, 1), tail(ty, 1)) / tail(norm, 1)
                direction <- if (abs(normal[1]) > 0.35) {
                  if (normal[1] < 0) "left" else "right"
                } else if (normal[2] < 0) "bottom" else "top"
                hjust <- c(left = 1, right = 0, top = 0.5, bottom = 0.5)[direction]
                vjust <- c(left = 0.5, right = 0.5, top = 0, bottom = 1)[direction]
                extent <- abs(normal[1]) * boxes$w[i] / 2 +
                  abs(normal[2]) * boxes$h[i] / 2
                centre <- (0.5 - hjust) * boxes$w[i] * normal[1] +
                  (0.5 - vjust) * boxes$h[i] * normal[2]
                position <- point + normal * (distance + extent - centre)
                candidate$text_x[i] <- position[1]
                candidate$text_y[i] <- position[2]
                candidate$hjust[i] <- hjust
                candidate$vjust[i] <- vjust
                box <- boxes[i, , drop = FALSE]
                box$x <- position[1]
                box$y <- position[2]
                box$cx <- position[1] + (0.5 - hjust) * box$w
                box$cy <- position[2] + (0.5 - vjust) * box$h
                box$xmin <- box$cx - box$bw / 2
                box$xmax <- box$cx + box$bw / 2
                box$ymin <- box$cy - box$bh / 2
                box$ymax <- box$cy + box$bh / 2
                if (any(ggchord_oriented_box_overlaps(box, local_boxes)) ||
                    any(ggchord_oriented_box_overlaps(box, repel_boxes))) next
                segments <- data.frame(
                  x0 = c(gl$anchor_x[i], gl$.radial_bend_x[i]),
                  y0 = c(gl$anchor_y[i], gl$.radial_bend_y[i]),
                  x1 = c(gl$.radial_bend_x[i], position[1]),
                  y1 = c(gl$.radial_bend_y[i], position[2]), group = i)
                # Unmoved labels are a single normal leader. Only displaced
                # labels retain the short normal departure and long fan segment.
                a <- c(segments$x1[1] - segments$x0[1], segments$y1[1] - segments$y0[1])
                b <- c(position[1] - gl$anchor_x[i], position[2] - gl$anchor_y[i])
                if (sum(a * b) > 0 && abs(a[1] * b[2] - a[2] * b[1]) <=
                    sin(3 * pi / 180) * sqrt(sum(a^2) * sum(b^2))) {
                  segments <- data.frame(x0 = gl$anchor_x[i], y0 = gl$anchor_y[i],
                    x1 = position[1], y1 = position[2], group = i)
                }
                crossed <- TRUE
                # Shorten the departure only if a neighbouring sequence blocks it.
                # At the limiting zero-length departure use a validated direct
                # connector, rather than failing an otherwise drawable map.
                for (fraction in c(1, 0.5, 0.2, 0.05, 0)) {
                  routed <- segments
                  if (nrow(routed) == 2L) {
                    bx <- gl$anchor_x[i] + fraction * (segments$x1[1] - gl$anchor_x[i])
                    by <- gl$anchor_y[i] + fraction * (segments$y1[1] - gl$anchor_y[i])
                    routed$x1[1] <- routed$x0[2] <- bx
                    routed$y1[1] <- routed$y0[2] <- by
                    if (fraction == 0) routed <- data.frame(
                      x0 = gl$anchor_x[i], y0 = gl$anchor_y[i],
                      x1 = position[1], y1 = position[2], group = i)
                  }
                  crossed <- any(vapply(seq_len(nrow(routed)), function(k) any(
                    ggchord_segments_cross(routed$x0[k], routed$y0[k],
                      routed$x1[k], routed$y1[k], local_paths$x0, local_paths$y0,
                      local_paths$x1, local_paths$y1)), logical(1)))
                  if (!crossed) {
                    segments <- routed
                    break
                  }
                }
                if (crossed) next
                candidate$.radial_bend_x[i] <- if (nrow(segments) == 1L) gl$anchor_x[i] else segments$x1[1]
                candidate$.radial_bend_y[i] <- if (nrow(segments) == 1L) gl$anchor_y[i] else segments$y1[1]
                local_boxes <- rbind(local_boxes, box)
                local_paths <- rbind(local_paths, segments)
                previous <- at
                local_parameter[i] <- at
                local_direction[i] <- direction
                cost <- cost + sum(sqrt((segments$x1 - segments$x0)^2 +
                                        (segments$y1 - segments$y0)^2))
                selected <- TRUE
                break
              }
              if (!selected) {
                solved <- FALSE
                break
              }
            }

            if (solved && (is.null(best) || cost < best$cost)) best <- list(
              labels = candidate, boxes = local_boxes, paths = local_paths,
              parameter = local_parameter, directions = local_direction,
              level = level, cost = cost)
          }
        }
      }
      if (!is.null(best)) break
    }
    if (is.null(best)) {
      if (.attempt < max(1L, 2L * length(groups) - 1L)) return(
        ggchord_radial_label_lanes(original, seq_arcs, side, units_per_inch,
          box_padding, point_padding, repel_boxes, .attempt + 1L, .fixed_paths,
          .minimum_levels, .balance))
      ggchord_stop("Cannot fit radial labels for ", gl$seq_id[rows[1]],
                   " on this device; enlarge the output or reduce label size")
    }
    gl <- best$labels
    occupied <- best$boxes
    paths <- best$paths
    parameter <- best$parameter
    directions <- best$directions
    tracks[rows] <- best$level + 1L
  }
  # Balance opposite horizontal sectors after finding their feasible compact
  # contours. Half the denser sector's extra clearance is a soft lower bound
  # for the other sector; all collisions and actual leaders are solved again.
  if (.balance && length(groups) > 1L) {
    sector <- vapply(groups, function(rows) {
      x <- mean(signs[rows] * source$outward_x[rows])
      y <- mean(signs[rows] * source$outward_y[rows])
      if (abs(y) < abs(x)) return("side")
      if (y > 0) "top" else "bottom"
    }, character(1))
    if (all(c("top", "bottom") %in% sector)) {
      extra <- vapply(groups, function(rows) max(tracks[rows] - 1L), numeric(1))
      target <- floor(max(extra[sector != "side"]) / 2)
      need <- sector != "side" & extra < target
      if (any(need)) {
        minimum <- extra * 0
        minimum[need] <- target
        return(ggchord_radial_label_lanes(original, seq_arcs, side, units_per_inch,
          box_padding, point_padding, repel_boxes, 0L, .fixed_paths,
          minimum, .balance = FALSE))
      }
    }
  }
  gl$.radial_parameter <- parameter
  list(labels = gl, lanes = paste(gl$seq_id, signs, sep = "\r"),
       directions = directions, tracks = tracks,
       draw_segment = seq_len(n) %in% active)
}

# Auto retains shared vertical columns only for locally left/right features.
# The other features use exactly the radial solver and its final leader paths.
ggchord_auto_label_lanes <- function(gl, seq_arcs, side = "outside",
                                      units_per_inch = 0.35,
                                      box_padding = 0.05,
                                      point_padding = 0.08,
                                      repel_boxes = NULL, .column_attempt = 0L) {
  frame <- ggchord_label_curve_frame(gl, seq_arcs)
  signs <- switch(side, outside = rep(1, nrow(gl)), inside = rep(-1, nrow(gl)),
                  ifelse(frame$signed_distance < 0, -1, 1))
  horizontal <- abs(frame$outward_x) >= abs(frame$outward_y)
  directions <- ifelse(signs * frame$outward_x < 0, "left", "right")
  columns <- gl
  columns$text[!horizontal] <- ""
  result <- ggchord_side_label_columns(columns, seq_arcs, side,
    units_per_inch, box_padding, point_padding + .column_attempt * 0.15, repel_boxes,
    directions = directions)
  boxes <- ggchord_text_boxes(result$labels,
    units_per_inch = units_per_inch, box_padding = box_padding)
  bends <- ggchord_elbow_bends(result$labels, boxes$w, boxes$h, directions)
  ids <- which(horizontal & !is.na(gl$text) & nzchar(gl$text))
  paths <- data.frame(
    x0 = c(gl$anchor_x[ids], bends$x[ids]),
    y0 = c(gl$anchor_y[ids], bends$y[ids]),
    x1 = c(bends$x[ids], result$labels$text_x[ids]),
    y1 = c(bends$y[ids], result$labels$text_y[ids]), group = rep(-ids, 2L))
  paths <- ggchord_collapse_crossed_elbows(
    transform(paths, group = -group), result$lanes)
  paths$group <- -paths$group
  # Carry the post-collapse bend into the shared renderer.
  for (i in ids) {
    first <- which(paths$group == -i)[1]
    bends$x[i] <- paths$x1[first]
    bends$y[i] <- paths$y1[first]
  }
  labels <- gl
  labels[horizontal, ] <- result$labels[horizontal, ]
  labels$.radial_bend_x <- bends$x
  labels$.radial_bend_y <- bends$y
  labels$.radial_parameter <- NA_real_
  radial_rows <- which(!horizontal)
  if (length(radial_rows)) {
    radial <- tryCatch(ggchord_radial_label_lanes(gl[radial_rows, , drop = FALSE],
      seq_arcs, side, units_per_inch, box_padding, point_padding,
      rbind(repel_boxes, boxes[ids, ]), .fixed_paths = paths), error = identity)
    if (inherits(radial, "error")) {
      if (.column_attempt < 6L && grepl("Cannot fit radial", conditionMessage(radial))) {
        return(ggchord_auto_label_lanes(gl, seq_arcs, side, units_per_inch,
          box_padding, point_padding, repel_boxes, .column_attempt + 1L))
      }
      stop(radial)
    }
    labels[radial_rows, names(radial$labels)] <- radial$labels
    result$tracks[radial_rows] <- radial$tracks
    result$directions[radial_rows] <- radial$directions
  }
  result$labels <- labels
  result$draw_segment <- !is.na(labels$text) & nzchar(labels$text)
  result
}
