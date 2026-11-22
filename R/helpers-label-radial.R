# Horizontal external labels on offsets of the actual sequence curves.
# Parameters are arc lengths, not angles around the global origin, so the
# same packing works for unequal radii, rotated curves and open sequences.
ggchord_radial_label_lanes <- function(gl, seq_arcs, side = "outside",
                                        units_per_inch = 0.35,
                                        box_padding = 0.04,
                                        point_padding = 0.05,
                                        repel_boxes = NULL, .attempt = 0L, .fixed_paths = NULL,
                                        .minimum_levels = numeric(), .balance = TRUE,
                                        .tangent_limit = NULL) {
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
  groups <- split(active, paste(gl$accver[active], signs[active], sep = "\r"))
  if (length(groups) > 1L) {
    order_groups <- seq_along(groups)
    if (.attempt %% 2L) order_groups <- rev(order_groups)
    rotation <- (.attempt %/% 2L) %% length(groups)
    if (rotation) order_groups <- c(tail(order_groups, -rotation), head(order_groups, rotation))
    groups <- groups[order_groups]
  }
  arc_ids <- vapply(seq_arcs, function(a) as.character(a$accver[1]), character(1))

  for (rows in groups) {
    arc <- seq_arcs[[match(gl$accver[rows[1]], arc_ids)]]
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
    interpolate <- function(v, at) approx(s, v, at, rule = 2)$y
    anchor_offset <- (gl$anchor_x[rows] - source$curve_x[rows]) *
      signs[rows] * source$outward_x[rows] +
      (gl$anchor_y[rows] - source$curve_y[rows]) * signs[rows] * source$outward_y[rows]
    clearance <- max(anchor_offset, 0) + point_padding + 0.18 * units_per_inch
    step <- max(0.015, 0.05 * units_per_inch)
    # Keep a small physical shoulder at each endpoint, rather than extending
    # the contour indefinitely into a neighbouring sequence's space.
    endpoint_padding <- 0.25 * units_per_inch
    shift_limit <- max(diff(range(s)) / 2, sum(boxes$bw[rows]) / 2)
    if (!is.null(.tangent_limit)) {
      shift_limit <- min(diff(range(s)) / 2,
        max(sum(boxes$bw[rows]) / 2, .tangent_limit))
    }
    shifts <- seq(0, shift_limit, by = step)
    offsets <- sort(unique(c(shifts, -shifts)))
    offsets <- offsets[order(abs(offsets), offsets)]
    # The normal and curve point depend only on the source label and offset.
    # Interpolate each row's allowed offsets once for all levels and orders.
    curve_end <- tail(s, 1L)
    slots <- lapply(preferred, function(origin) {
      at <- unname(origin + offsets)
      keep <- which(at >= -endpoint_padding & at <= curve_end + endpoint_padding)
      at <- at[keep]
      normal_x <- interpolate(nx, at)
      normal_y <- interpolate(ny, at)
      norm_xy <- sqrt(normal_x^2 + normal_y^2)
      normal_x <- normal_x / norm_xy
      normal_y <- normal_y / norm_xy
      point_x <- interpolate(arc$x, at)
      point_y <- interpolate(arc$y, at)
      low <- at < 0
      high <- at > curve_end
      point_x[low] <- point_x[low] + at[low] * tx[1L] / norm[1L]
      point_y[low] <- point_y[low] + at[low] * ty[1L] / norm[1L]
      point_x[high] <- point_x[high] + (at[high] - curve_end) * tail(tx, 1L) / tail(norm, 1L)
      point_y[high] <- point_y[high] + (at[high] - curve_end) * tail(ty, 1L) / tail(norm, 1L)
      direction <- ifelse(abs(normal_x) > .35, ifelse(normal_x < 0, "left", "right"),
        ifelse(normal_y < 0, "bottom", "top"))
      list(at = at, offset = offsets[keep], nx = normal_x, ny = normal_y,
        x = point_x, y = point_y, direction = direction,
        hjust = unname(c(left = 1, right = 0, top = .5, bottom = .5)[direction]),
        vjust = unname(c(left = .5, right = .5, top = 0, bottom = 1)[direction]))
    })
    best <- NULL
    group_key <- paste(gl$accver[rows[1]], signs[rows[1]], sep = "\r")
    minimum_level <- .minimum_levels[group_key]
    if (!length(minimum_level) || is.na(minimum_level)) minimum_level <- 0
    levels <- c(0:12, 16, 24, 40)
    for (level in levels[levels >= minimum_level]) {
      distance <- clearance + level * max(0.08, 0.25 * units_per_inch)
      # At one distance all layout trials share the exact same candidate
      # rectangles. Only the accepted obstacles and the search order change.
      rectangles <- lapply(rows, function(i) {
        slot <- slots[[as.character(i)]]
        extent <- abs(slot$nx) * boxes$w[i] / 2 + abs(slot$ny) * boxes$h[i] / 2
        centre <- (.5 - slot$hjust) * boxes$w[i] * slot$nx +
          (.5 - slot$vjust) * boxes$h[i] * slot$ny
        x <- slot$x + slot$nx * (distance + extent - centre)
        y <- slot$y + slot$ny * (distance + extent - centre)
        box <- boxes[rep.int(i, length(x)), , drop = FALSE]
        box$x <- x
        box$y <- y
        box$cx <- x + (.5 - slot$hjust) * box$w
        box$cy <- y + (.5 - slot$vjust) * box$h
        box$xmin <- box$cx - box$bw / 2
        box$xmax <- box$cx + box$bw / 2
        box$ymin <- box$cy - box$bh / 2
        box$ymax <- box$cy + box$bh / 2
        box
      })
      names(rectangles) <- rows
      for (reverse in c(FALSE, TRUE)) {
        spreads <- if (length(rows) == 1L) 1 else c(1, 1.25, 1.5, 2)
        for (spread in spreads) {
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
              slot <- slots[[as.character(i)]]
              ranked <- order(abs(slot$offset - bias), slot$offset)
              ordered_slots <- if (reverse) ranked[slot$at[ranked] < previous - 1e-8] else
                ranked[slot$at[ranked] > previous + 1e-8]
              for (j in ordered_slots) {
                at <- slot$at[j]
                box <- rectangles[[as.character(i)]]
                if (any(ggchord_oriented_box_overlaps(box, local_boxes, i = j)) ||
                    any(ggchord_oriented_box_overlaps(box, repel_boxes, i = j))) next
                x <- box$x[j]
                y <- box$y[j]
                ax <- gl$anchor_x[i]
                ay <- gl$anchor_y[i]
                dx <- gl$.radial_bend_x[i] - ax
                dy <- gl$.radial_bend_y[i] - ay
                tx_end <- x - ax
                ty_end <- y - ay
                direct <- dx * tx_end + dy * ty_end > 0 &&
                  abs(dx * ty_end - dy * tx_end) <= sin(3 * pi / 180) *
                  sqrt((dx^2 + dy^2) * (tx_end^2 + ty_end^2))
                clear <- FALSE
                # The first collision-free departure wins. The direct route
                # is identical for every fraction and needs just one test.
                for (fraction in if (direct) 0 else c(1, .5, .2, .05, 0)) {
                  bx <- ax + fraction * dx
                  by <- ay + fraction * dy
                  x0 <- if (fraction == 0) ax else c(ax, bx)
                  y0 <- if (fraction == 0) ay else c(ay, by)
                  x1 <- if (fraction == 0) x else c(bx, x)
                  y1 <- if (fraction == 0) y else c(by, y)
                  clear <- !any(vapply(seq_along(x0), function(k) any(
                    ggchord_segments_cross(x0[k], y0[k], x1[k], y1[k],
                      local_paths$x0, local_paths$y0, local_paths$x1, local_paths$y1)), logical(1)))
                  if (clear) break
                }
                if (!clear) next
                candidate$text_x[i] <- x
                candidate$text_y[i] <- y
                candidate$hjust[i] <- slot$hjust[j]
                candidate$vjust[i] <- slot$vjust[j]
                candidate$.radial_bend_x[i] <- bx
                candidate$.radial_bend_y[i] <- by
                local_boxes <- rbind(local_boxes, box[j, , drop = FALSE])
                local_paths <- rbind(local_paths, list2DF(list(x0 = x0, y0 = y0,
                  x1 = x1, y1 = y1, group = rep.int(i, length(x0)))))
                previous <- at
                local_parameter[i] <- at
                local_direction[i] <- slot$direction[j]
                cost <- cost + sum(sqrt((x1 - x0)^2 + (y1 - y0)^2))
                selected <- TRUE
                break
              }
              # Later segments add nonnegative length; equal cost keeps the
              # earlier solution in the original deterministic search order.
              if (!selected || (!is.null(best) && cost >= best$cost)) {
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
          .minimum_levels, .balance, .tangent_limit))
      ggchord_stop("Cannot fit radial labels for ", gl$accver[rows[1]],
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
          minimum, .balance = FALSE, .tangent_limit = .tangent_limit))
      }
    }
  }
  gl$.radial_parameter <- parameter
  list(labels = gl, lanes = paste(gl$accver, signs, sep = "\r"),
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
      ggchord_rethrow(radial)
    }
    labels[radial_rows, names(radial$labels)] <- radial$labels
    result$tracks[radial_rows] <- radial$tracks
    result$directions[radial_rows] <- radial$directions
  }
  result$labels <- labels
  result$draw_segment <- !is.na(labels$text) & nzchar(labels$text)
  result
}
