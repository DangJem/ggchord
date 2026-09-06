ggchord_uncross_labels <- function(gl, lanes = NULL, max_swaps = NULL,
                                   endpoint_fun = NULL) {
  n <- nrow(gl)
  if (n < 2) return(list(labels = gl, swaps = 0L))
  if (is.null(lanes)) lanes <- gl$accver %||% rep("all", n)
  if (length(lanes) != n) {
    ggchord_stop("lanes must have one value per gene label")
  }

  active <- !is.na(gl$text) & nzchar(gl$text)
  lane_rows <- split(which(active), as.character(lanes[active]), drop = TRUE)
  total_swaps <- 0L

  for (rows in lane_rows) {
    nr <- length(rows)
    if (nr < 2) next
    swap_limit <- max_swaps %||% max(100L, 4L * nr^2)
    changed <- TRUE
    lane_swaps <- 0L

    while (changed && lane_swaps < swap_limit) {
      changed <- FALSE
      endpoints <- if (is.null(endpoint_fun)) {
        data.frame(x = gl$text_x, y = gl$text_y)
      } else {
        endpoint_fun(gl)
      }
      for (ii in seq_len(nr - 1L)) {
        for (jj in (ii + 1L):nr) {
          i <- rows[ii]
          j <- rows[jj]
          crossed <- ggchord_segments_cross(
            gl$anchor_x[i], gl$anchor_y[i], endpoints$x[i], endpoints$y[i],
            gl$anchor_x[j], gl$anchor_y[j], endpoints$x[j], endpoints$y[j]
          )
          if (!crossed) next

          if (is.null(endpoint_fun)) {
            candidate_x_i <- gl$text_x[j]
            candidate_y_i <- gl$text_y[j]
            candidate_x_j <- gl$text_x[i]
            candidate_y_j <- gl$text_y[i]
            candidate_endpoints <- NULL
          } else {
            candidate <- gl
            candidate$text_x[c(i, j)] <- gl$text_x[c(j, i)]
            candidate$text_y[c(i, j)] <- gl$text_y[c(j, i)]
            candidate_endpoints <- endpoint_fun(candidate)
          }
          old_length <- sqrt((endpoints$x[i] - gl$anchor_x[i])^2 +
                               (endpoints$y[i] - gl$anchor_y[i])^2) +
            sqrt((endpoints$x[j] - gl$anchor_x[j])^2 +
                   (endpoints$y[j] - gl$anchor_y[j])^2)
          new_length <- if (is.null(endpoint_fun)) {
            sqrt((candidate_x_i - gl$anchor_x[i])^2 +
                   (candidate_y_i - gl$anchor_y[i])^2) +
              sqrt((candidate_x_j - gl$anchor_x[j])^2 +
                     (candidate_y_j - gl$anchor_y[j])^2)
          } else {
            sqrt((candidate_endpoints$x[i] - gl$anchor_x[i])^2 +
                   (candidate_endpoints$y[i] - gl$anchor_y[i])^2) +
              sqrt((candidate_endpoints$x[j] - gl$anchor_x[j])^2 +
                     (candidate_endpoints$y[j] - gl$anchor_y[j])^2)
          }
          if (new_length >= old_length - 1e-10) next

          if (is.null(endpoint_fun)) {
            gl$text_x[c(i, j)] <- c(candidate_x_i, candidate_x_j)
            gl$text_y[c(i, j)] <- c(candidate_y_i, candidate_y_j)
          } else {
            gl <- candidate
          }
          lane_swaps <- lane_swaps + 1L
          total_swaps <- total_swaps + 1L
          changed <- TRUE
          break
        }
        if (changed) break
      }
    }
  }

  list(labels = gl, swaps = total_swaps)
}

ggchord_elbow_bends <- function(gl, text_widths, text_heights = NULL,
                                directions = NULL) {
  n <- nrow(gl)
  if (n == 0) return(data.frame(x = numeric(0), y = numeric(0)))
  if (is.null(text_heights)) text_heights <- text_widths
  if (is.null(directions)) {
    directions <- ifelse(gl$hjust < 0.5, "right", "left")
  }
  horiz <- abs(gl$text_x - gl$anchor_x)
  vert <- abs(gl$text_y - gl$anchor_y)
  horizontal_stub <- directions %in% c("left", "right")
  span <- ifelse(horizontal_stub, horiz, vert)
  extent <- ifelse(horizontal_stub, text_widths, text_heights)
  stub_len <- pmin(pmax(0.02, 0.3 * span), pmax(0.3 * extent, 0.04))

  bx <- gl$text_x
  by <- gl$text_y
  left <- directions == "left"
  right <- directions == "right"
  bottom <- directions == "bottom"
  top <- directions == "top"
  bx[left] <- pmin(gl$text_x[left] + stub_len[left], gl$anchor_x[left])
  bx[right] <- pmax(gl$text_x[right] - stub_len[right], gl$anchor_x[right])
  by[bottom] <- pmin(gl$text_y[bottom] + stub_len[bottom], gl$anchor_y[bottom])
  by[top] <- pmax(gl$text_y[top] - stub_len[top], gl$anchor_y[top])
  data.frame(x = bx, y = by)
}

ggchord_label_curve_frame <- function(gl, seq_arcs) {
  n <- nrow(gl)
  frame <- data.frame(
    curve_x = numeric(n), curve_y = numeric(n),
    outward_x = numeric(n), outward_y = numeric(n),
    signed_distance = numeric(n),
    curve_index = integer(n)
  )
  if (n == 0) return(frame)

  arc_ids <- vapply(seq_arcs, function(a) {
    if (nrow(a) == 0) "" else as.character(unique(a$accver)[1])
  }, character(1))

  for (sid in unique(gl$accver)) {
    rows <- which(gl$accver == sid)
    arc_pos <- match(sid, arc_ids)
    if (is.na(arc_pos) || nrow(seq_arcs[[arc_pos]]) < 2) next
    arc <- seq_arcs[[arc_pos]]

    for (i in rows) {
      # Locate the label by its fixed gene anchor, not by the repelled text.
      # This prevents a far-moving label from snapping to another part of a
      # highly curved sequence path.
      d2 <- (arc$x - gl$anchor_x[i])^2 + (arc$y - gl$anchor_y[i])^2
      k <- which.min(d2)
      k0 <- max(1L, k - 1L)
      k1 <- min(nrow(arc), k + 1L)
      tx <- arc$x[k1] - arc$x[k0]
      ty <- arc$y[k1] - arc$y[k0]
      tangent_length <- sqrt(tx^2 + ty^2)
      if (!is.finite(tangent_length) || tangent_length < 1e-12) next

      # Reference paths follow increasing genomic angle. Their right normal
      # is the outside track, including concave curves and paths crossing the
      # origin; an origin-dot-product test would flip sides mid-sequence.
      path_direction <- attr(arc, "ggchord_path_direction") %||% 1
      nx <- path_direction * ty / tangent_length
      ny <- -path_direction * tx / tangent_length
      frame$curve_x[i] <- arc$x[k]
      frame$curve_y[i] <- arc$y[k]
      frame$outward_x[i] <- nx
      frame$outward_y[i] <- ny
      frame$curve_index[i] <- k
      frame$signed_distance[i] <-
        (gl$text_x[i] - arc$x[k]) * nx +
        (gl$text_y[i] - arc$y[k]) * ny
    }
  }
  frame
}

ggchord_enforce_label_side <- function(gl, seq_arcs, side = "auto") {
  if (nrow(gl) == 0 || identical(side, "auto")) return(gl)
  frame <- ggchord_label_curve_frame(gl, seq_arcs)
  want_inside <- identical(side, "inside")
  flip <- (frame$signed_distance < 0) != want_inside
  flip[!is.finite(frame$signed_distance)] <- FALSE
  if (!any(flip)) return(gl)

  # Reflect across the actual local tangent of the sequence curve. Unlike a
  # radius-from-origin approximation, this remains correct when seq_radius,
  # seq_curvature, seq_gap, seq_order or global rotation change the geometry.
  correction <- 2 * frame$signed_distance[flip]
  gl$text_x[flip] <- gl$text_x[flip] -
    correction * frame$outward_x[flip]
  gl$text_y[flip] <- gl$text_y[flip] -
    correction * frame$outward_y[flip]
  gl
}

# Pack one-dimensional label anchors while preserving their gene order.
ggchord_pack_label_axis <- function(preferred, order_value,
                                    before, after, gap = 0) {
  n <- length(preferred)
  if (n < 2) return(preferred)
  ord <- order(order_value, seq_len(n))
  position <- preferred[ord]
  before <- before[ord]
  after <- after[ord]

  # A forward pass creates the minimum legal spacing. Translating the whole
  # lane afterwards is the least-squares fit back to the preferred positions
  # and does not change any of those spacings.
  for (i in 2:n) {
    position[i] <- max(
      position[i], position[i - 1L] + after[i - 1L] + before[i] + gap
    )
  }
  position <- position + mean(preferred[ord] - position)
  out <- numeric(n)
  out[ord] <- position
  out
}

# Pack away from a shared lane centre. The label nearest that centre remains
# close to its gene, while congestion is absorbed by labels farther towards
# either end of the sequence. This avoids giving adjacent labels parallel,
# similarly offset leaders that can intersect on a curved arc.
ggchord_pack_label_axis_outward <- function(preferred, order_value,
                                            before, after, centre,
                                            gap = 0) {
  n <- length(preferred)
  if (n < 2) return(preferred)
  ord <- order(order_value, seq_len(n))
  position <- preferred[ord]
  before <- before[ord]
  after <- after[ord]
  pivot <- which.min(abs(order_value[ord] - centre))

  if (pivot > 1) {
    for (i in seq.int(pivot - 1L, 1L)) {
      position[i] <- min(
        position[i], position[i + 1L] - after[i] - before[i + 1L] - gap
      )
    }
  }
  if (pivot < n) {
    for (i in seq.int(pivot + 1L, n)) {
      position[i] <- max(
        position[i], position[i - 1L] + after[i - 1L] + before[i] + gap
      )
    }
  }
  out <- numeric(n)
  out[ord] <- position
  out
}

# Put horizontal labels on compact cardinal lanes around their own sequence.
# Shared left/right columns used by the auto layout.
ggchord_side_label_columns <- function(gl, seq_arcs,
                                        side = "outside",
                                        units_per_inch = 0.35,
                                        box_padding = 0.25,
                                        point_padding = 0.1,
                                        repel_boxes = NULL,
                                        max_iter = 100, directions = NULL) {
  n <- nrow(gl)
  if (n == 0) return(list(labels = gl, lanes = character(0)))

  active <- !is.na(gl$text) & nzchar(gl$text)
  source_frame <- ggchord_label_curve_frame(gl, seq_arcs)
  anchor_labels <- gl
  anchor_labels$text_x <- anchor_labels$anchor_x
  anchor_labels$text_y <- anchor_labels$anchor_y
  anchor_frame <- ggchord_label_curve_frame(anchor_labels, seq_arcs)
  side_sign <- if (identical(side, "outside")) {
    rep(1, n)
  } else if (identical(side, "inside")) {
    rep(-1, n)
  } else {
    ifelse(source_frame$signed_distance < 0, -1, 1)
  }
  side_sign[!is.finite(side_sign)] <- 1
  desired_x <- side_sign * anchor_frame$outward_x
  desired_y <- side_sign * anchor_frame$outward_y
  # Give every sequence-side group one primary cardinal rail. This keeps the
  # labels belonging to one sequence visually close and prevents neighbouring
  # sequences from competing for the same corners. Parallel rows absorb genuine
  # congestion without treating unused sides as a target to distribute toward.
  if (is.null(directions)) directions <- ifelse(desired_x < 0, "left", "right")
  directions[!active] <- NA_character_
  base_lanes <- paste(gl$accver, directions, sep = "\r")

  gl$text_angle[active] <- 0
  gl$hjust[active] <- c(left = 1, right = 0, top = 0.5, bottom = 0.5)[
    directions[active]
  ]
  gl$vjust[active] <- c(left = 0.5, right = 0.5, top = 0, bottom = 1)[
    directions[active]
  ]
  axis_gap <- max(0.02, 0.025 * units_per_inch)
  base_lane_rows <- split(which(active), base_lanes[active], drop = TRUE)
  rail_gap <- point_padding + box_padding * units_per_inch +
    max(0.04, 0.04 * units_per_inch)
  lanes <- base_lanes
  tracks <- rep(NA_integer_, n)

  # All labels assigned to the same side share one compact vertical column.
  # Packing from feature anchors, rather than previous-pass label positions,
  # avoids layout feedback, preserves their global vertical order and prevents
  # two sequence-specific columns from sending leaders across one another.
  vertical_rows_all <- which(active & directions %in% c("left", "right"))
  vertical_groups <- split(
    vertical_rows_all, directions[vertical_rows_all], drop = TRUE
  )
  template_boxes <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = box_padding
  )
  for (rows in vertical_groups) {
    direction <- directions[rows[1]]
    before <- template_boxes$y[rows] - template_boxes$ymin[rows]
    after <- template_boxes$ymax[rows] - template_boxes$y[rows]
    gl$text_y[rows] <- ggchord_pack_label_axis(
      gl$anchor_y[rows], gl$anchor_y[rows], before, after, gap = axis_gap
    )
    gl$text_x[rows] <- if (direction == "left") {
      min(anchor_frame$curve_x[rows]) - rail_gap
    } else {
      max(anchor_frame$curve_x[rows]) + rail_gap
    }
    tracks[rows] <- 1L
    lanes[rows] <- direction
  }

  # A shared side column must also use a shared endpoint order. Curved and
  # rotated sequences can place two feature anchors in an order that differs
  # from their first packed label positions even when those positions do not
  # overlap. Swap only crossing endpoints, then repack the resulting slots for
  # their actual text heights. Repeating this small 2-opt pass produces a
  # deterministic, compact column whose direct leaders do not cross.
  if (length(vertical_rows_all) > 1L) {
    for (pass in seq_len(8L)) {
      candidate <- gl[vertical_rows_all, , drop = FALSE]
      uncrossed <- ggchord_uncross_labels(
        candidate, lanes = directions[vertical_rows_all]
      )
      if (uncrossed$swaps == 0L) break
      candidate <- uncrossed$labels
      candidate_boxes <- ggchord_text_boxes(
        candidate, units_per_inch = units_per_inch,
        box_padding = box_padding
      )
      candidate_groups <- split(
        seq_along(vertical_rows_all),
        directions[vertical_rows_all], drop = TRUE
      )
      for (local_rows in candidate_groups) {
        before <- candidate_boxes$y[local_rows] -
          candidate_boxes$ymin[local_rows]
        after <- candidate_boxes$ymax[local_rows] -
          candidate_boxes$y[local_rows]
        candidate$text_y[local_rows] <- ggchord_pack_label_axis(
          candidate$text_y[local_rows], candidate$text_y[local_rows],
          before, after, gap = axis_gap
        )
      }
      gl$text_x[vertical_rows_all] <- candidate$text_x
      gl$text_y[vertical_rows_all] <- candidate$text_y
    }
  }

  move_lane_outward <- function(rows, amount) {
    if (!is.finite(amount) || amount <= 0) return()
    direction <- directions[rows[1]]
    if (direction == "left") gl$text_x[rows] <<- gl$text_x[rows] - amount
    if (direction == "right") gl$text_x[rows] <<- gl$text_x[rows] + amount
    if (direction == "bottom") gl$text_y[rows] <<- gl$text_y[rows] - amount
    if (direction == "top") gl$text_y[rows] <<- gl$text_y[rows] + amount
  }

  ensure_vertical_side <- function() {
    for (rows in vertical_groups) {
      direction <- directions[rows[1]]
      dx <- gl$text_x[rows] - anchor_frame$curve_x[rows]
      dy <- gl$text_y[rows] - anchor_frame$curve_y[rows]
      signed <- side_sign[rows] *
        (dx * anchor_frame$outward_x[rows] +
           dy * anchor_frame$outward_y[rows])
      projection <- switch(
        direction,
        left = -side_sign[rows] * anchor_frame$outward_x[rows],
        right = side_sign[rows] * anchor_frame$outward_x[rows],
        bottom = -side_sign[rows] * anchor_frame$outward_y[rows],
        top = side_sign[rows] * anchor_frame$outward_y[rows]
      )
      usable <- is.finite(projection) & projection > 0.25
      if (!any(usable)) next
      target <- max(0.015, 0.025 * units_per_inch)
      amount <- max((target - signed[usable]) / projection[usable], 0)
      # The rail coordinates already place it beyond the relevant curve
      # extreme. This cap is only a local-side correction, never a second
      # layout offset.
      move_lane_outward(rows, min(amount, rail_gap))
    }
  }
  ensure_vertical_side()

  # Sequence and axis text are fixed obstacles. Resolve those conflicts by
  # moving the complete rail outwards, which preserves its alignment and the
  # order of every leader instead of pushing individual labels off the rail.
  if (!is.null(repel_boxes) && nrow(repel_boxes) > 0) {
    for (pass in seq_len(8)) {
      moved <- FALSE
      boxes <- ggchord_text_boxes(
        gl, units_per_inch = units_per_inch, box_padding = 0
      )
      for (rows in vertical_groups) {
        overlap_x <- outer(
          boxes$xmin[rows], repel_boxes$xmax,
          function(a, b) a < b - 1e-7
        ) & outer(
          boxes$xmax[rows], repel_boxes$xmin,
          function(a, b) a > b + 1e-7
        )
        overlap_y <- outer(
          boxes$ymin[rows], repel_boxes$ymax,
          function(a, b) a < b - 1e-7
        ) & outer(
          boxes$ymax[rows], repel_boxes$ymin,
          function(a, b) a > b + 1e-7
        )
        hits <- which(overlap_x & overlap_y, arr.ind = TRUE)
        if (nrow(hits) == 0) next
        direction <- directions[rows[1]]
        amount <- switch(
          direction,
          left = max(boxes$xmax[rows[hits[, 1]]] -
                       repel_boxes$xmin[hits[, 2]] + axis_gap),
          right = max(repel_boxes$xmax[hits[, 2]] -
                        boxes$xmin[rows[hits[, 1]]] + axis_gap),
          bottom = max(boxes$ymax[rows[hits[, 1]]] -
                         repel_boxes$ymin[hits[, 2]] + axis_gap),
          top = max(repel_boxes$ymax[hits[, 2]] -
                      boxes$ymin[rows[hits[, 1]]] + axis_gap)
        )
        move_lane_outward(rows, amount)
        moved <- TRUE
      }
      if (!moved) break
    }
  }
  ensure_vertical_side()


  list(
    labels = gl,
    lanes = lanes,
    directions = directions,
    tracks = tracks
  )
}

# Put labels on the nearest collision-free local offset track. Unlike a
# radius-from-origin layout, every candidate position is measured from the
# actual sequence curve and its local outward normal, so straight, strongly
# curved and differently sized sequences use the same algorithm.
ggchord_offset_label_tracks <- function(gl, seq_arcs,
                                        side = "outside",
                                        orientation = c("horizontal", "arc"),
                                        units_per_inch = 0.35,
                                        box_padding = 0.18,
                                        point_padding = 0.08,
                                        repel_boxes = NULL) {
  orientation <- match.arg(orientation)
  n <- nrow(gl)
  if (n == 0) {
    return(list(labels = gl, lanes = character(0),
                directions = character(0), tracks = integer(0),
                draw_segment = logical(0)))
  }

  active <- !is.na(gl$text) & nzchar(gl$text)
  source_frame <- ggchord_label_curve_frame(gl, seq_arcs)
  anchor_labels <- gl
  anchor_labels$text_x <- anchor_labels$anchor_x
  anchor_labels$text_y <- anchor_labels$anchor_y
  frame <- ggchord_label_curve_frame(anchor_labels, seq_arcs)
  side_sign <- if (identical(side, "outside")) {
    rep(1, n)
  } else if (identical(side, "inside")) {
    rep(-1, n)
  } else {
    ifelse(source_frame$signed_distance < 0, -1, 1)
  }
  side_sign[!is.finite(side_sign)] <- 1
  normal_x <- side_sign * frame$outward_x
  normal_y <- side_sign * frame$outward_y
  directions <- ifelse(
    abs(normal_x) >= abs(normal_y),
    ifelse(normal_x < 0, "left", "right"),
    ifelse(normal_y < 0, "bottom", "top")
  )

  if (identical(orientation, "horizontal")) {
    gl$text_angle[active] <- 0
    gl$hjust[active] <- c(left = 1, right = 0, top = 0.5, bottom = 0.5)[
      directions[active]
    ]
    gl$vjust[active] <- c(left = 0.5, right = 0.5, top = 0, bottom = 1)[
      directions[active]
    ]
  } else {
    tangent_x <- -frame$outward_y
    tangent_y <- frame$outward_x
    angle <- (atan2(tangent_y, tangent_x) * 180 / pi + 360) %% 360
    upside_down <- angle > 90 & angle < 270
    angle[upside_down] <- (angle[upside_down] + 180) %% 360
    gl$text_angle[active] <- angle[active]
    gl$hjust[active] <- 0.5
    gl$vjust[active] <- 0.5
  }

  # Measure the anchor-relative box offset once. The innermost edge of every
  # first-track label is then placed at the same visual clearance from its
  # own sequence curve, even when text justification differs by quadrant.
  gl$text_x <- frame$curve_x
  gl$text_y <- frame$curve_y
  templates <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = box_padding
  )
  centre_offset <-
    (templates$cx - templates$x) * normal_x +
    (templates$cy - templates$y) * normal_y
  text_angle <- gl$text_angle * pi / 180
  text_x_axis_x <- cos(text_angle)
  text_x_axis_y <- sin(text_angle)
  text_y_axis_x <- -sin(text_angle)
  text_y_axis_y <- cos(text_angle)
  # Project the oriented rectangle itself, not its axis-aligned bounding box.
  # Projecting the latter makes a long tangent-aligned arc label look long in
  # the normal direction as well and can push it several unnecessary tracks
  # away from the sequence.
  normal_extent <-
    templates$w * abs(text_x_axis_x * normal_x +
                        text_x_axis_y * normal_y) / 2 +
    templates$h * abs(text_y_axis_x * normal_x +
                        text_y_axis_y * normal_y) / 2 +
    box_padding * units_per_inch
  clearance <- point_padding + max(0.035, 0.04 * units_per_inch)
  base_distance <- pmax(
    clearance - centre_offset + normal_extent,
    clearance
  )

  base_lanes <- paste(gl$accver, side_sign, sep = "\r")
  lane_rows <- split(which(active), base_lanes[active], drop = TRUE)
  tracks <- rep(NA_integer_, n)
  placed_boxes <- NULL
  axis_gap <- max(0.015, 0.02 * units_per_inch)

  overlaps_boxes <- function(candidate, other) {
    any(ggchord_oriented_box_overlaps(candidate, other))
  }

  for (rows in lane_rows) {
    rows <- rows[order(frame$curve_index[rows], rows)]
    lane_step <- max(
      2 * normal_extent[rows] + axis_gap,
      0.08 + 0.08 * units_per_inch,
      na.rm = TRUE
    )
    for (i in rows) {
      selected <- FALSE
      # At most one new track per active label is needed when boxes are
      # finite, with a small reserve for fixed sequence/axis obstacles.
      for (track in seq_len(length(rows) + 8L)) {
        distance <- base_distance[i] + (track - 1L) * lane_step
        # Keep the text centre on its feature's local normal. Tangential
        # nudging can reverse two neighbouring labels and consequently make
        # their leaders cross, even when both text boxes remain disjoint.
        candidate_label <- gl[i, , drop = FALSE]
        candidate_label$text_x <- frame$curve_x[i] + normal_x[i] * distance
        candidate_label$text_y <- frame$curve_y[i] + normal_y[i] * distance
        candidate_box <- ggchord_text_boxes(
          candidate_label,
          units_per_inch = units_per_inch,
          box_padding = box_padding
        )
        if (!overlaps_boxes(candidate_box, placed_boxes) &&
            !overlaps_boxes(candidate_box, repel_boxes)) {
          gl$text_x[i] <- candidate_label$text_x
          gl$text_y[i] <- candidate_label$text_y
          tracks[i] <- track
          placed_boxes <- rbind(placed_boxes, candidate_box)
          selected <- TRUE
        }
        if (selected) break
      }
      if (!selected) {
        track <- length(rows) + 9L
        distance <- base_distance[i] + (track - 1L) * lane_step
        gl$text_x[i] <- frame$curve_x[i] + normal_x[i] * distance
        gl$text_y[i] <- frame$curve_y[i] + normal_y[i] * distance
        tracks[i] <- track
        placed_boxes <- rbind(
          placed_boxes,
          ggchord_text_boxes(
            gl[i, , drop = FALSE], units_per_inch = units_per_inch,
            box_padding = box_padding
          )
        )
      }
    }
  }

  lanes <- paste(base_lanes, tracks, sep = "\r")
  draw_segment <- active
  if (identical(orientation, "arc")) {
    draw_segment <- active & !is.na(tracks) & tracks > 1L
  }
  list(
    labels = gl,
    lanes = lanes,
    directions = directions,
    tracks = tracks,
    draw_segment = draw_segment
  )
}

# Hide only labels that remain conflicted after a deterministic layout. The
# default max_overlaps = Inf therefore retains every label, while a finite
# value behaves as a final decluttering threshold without influencing any
# successfully placed label coordinates.
ggchord_hide_conflicted_labels <- function(gl, max_overlaps = Inf,
                                            units_per_inch = 0.35,
                                            repel_boxes = NULL) {
  if (!is.finite(max_overlaps) || nrow(gl) == 0) return(gl)
  counts <- ggchord_label_conflict_counts(
    gl, units_per_inch = units_per_inch, repel_boxes = repel_boxes
  )
  gl$text[counts > max_overlaps] <- NA_character_
  gl
}

ggchord_label_box_conflicts <- function(gl, units_per_inch = 0.35,
                                        box_padding = 0.25,
                                        repel_boxes = NULL,
                                        tol = 1e-7) {
  active <- !is.na(gl$text) & nzchar(gl$text)
  gl <- gl[active, , drop = FALSE]
  n <- nrow(gl)
  if (n == 0) return(FALSE)

  boxes <- ggchord_text_boxes(
    gl, units_per_inch = units_per_inch, box_padding = box_padding
  )
  if (n > 1) {
    dx <- abs(outer(boxes$cx, boxes$cx, "-"))
    dy <- abs(outer(boxes$cy, boxes$cy, "-"))
    overlap <- upper.tri(dx) &
      dx < outer(boxes$bw, boxes$bw, "+") / 2 - tol &
      dy < outer(boxes$bh, boxes$bh, "+") / 2 - tol
    if (any(overlap)) return(TRUE)
  }

  if (!is.null(repel_boxes) && nrow(repel_boxes) > 0) {
    overlap_x <- outer(boxes$xmin, repel_boxes$xmax, function(a, b) a < b - tol) &
      outer(boxes$xmax, repel_boxes$xmin, function(a, b) a > b + tol)
    overlap_y <- outer(boxes$ymin, repel_boxes$ymax, function(a, b) a < b - tol) &
      outer(boxes$ymax, repel_boxes$ymin, function(a, b) a > b + tol)
    if (any(overlap_x & overlap_y)) return(TRUE)
  }

  FALSE
}

# Collapse only elbow stubs that take part in a crossing. This includes corner
# crossings between different cardinal rails; the corresponding leader becomes
# straight, while all conflict-free elbows retain their bend and stub lengths.
ggchord_collapse_crossed_elbows <- function(segments, lanes,
                                            max_passes = NULL) {
  if (nrow(segments) < 2 || length(lanes) == 0) return(segments)
  groups <- unique(segments$group)
  max_passes <- max_passes %||% length(groups)

  for (pass in seq_len(max_passes)) {
    bad <- integer(0)
    for (i in seq_len(nrow(segments) - 1L)) {
      gi <- segments$group[i]
      for (j in (i + 1L):nrow(segments)) {
        gj <- segments$group[j]
        if (gi == gj || is.na(lanes[gi]) || is.na(lanes[gj])) next
        if (ggchord_segments_cross(
          segments$x0[i], segments$y0[i], segments$x1[i], segments$y1[i],
          segments$x0[j], segments$y0[j], segments$x1[j], segments$y1[j]
        )) {
          bad <- c(bad, gi, gj)
        }
      }
    }
    bad <- unique(bad)
    if (length(bad) == 0) break

    changed <- FALSE
    for (g in bad) {
      rows <- which(segments$group == g)
      if (length(rows) < 2) next
      first <- rows[1]
      stub <- rows[length(rows)]
      if (segments$x0[stub] == segments$x1[stub] &&
          segments$y0[stub] == segments$y1[stub]) next
      # The final row ends at the text anchor. Move the preceding bend there
      # and collapse the horizontal stub to a zero-length segment.
      segments$x1[first] <- segments$x1[stub]
      segments$y1[first] <- segments$y1[stub]
      segments$x0[stub] <- segments$x1[stub]
      segments$y0[stub] <- segments$y1[stub]
      changed <- TRUE
    }
    if (!changed) break
  }

  segments
}

#' Rebuild straight leader-line segments after a final label de-overlap pass.
#' @noRd
ggchord_repel_segments <- function(gl, min_segment_length = 0.5) {
  n <- nrow(gl)
  empty <- data.frame(x0 = numeric(0), y0 = numeric(0),
                      x1 = numeric(0), y1 = numeric(0),
                      group = integer(0), stringsAsFactors = FALSE)
  if (n == 0) return(empty)

  if ("anchor_x" %in% names(gl)) {
    ax <- gl$anchor_x
    ay <- gl$anchor_y
  } else {
    ax <- gl$text_x
    ay <- gl$text_y
  }
  seg_dist <- sqrt((gl$text_x - ax)^2 + (gl$text_y - ay)^2)
  visible <- if ("text" %in% names(gl)) {
    !is.na(gl$text) & nzchar(gl$text)
  } else {
    rep(TRUE, n)
  }
  keep_seg <- seg_dist > min_segment_length & visible
  data.frame(
    x0 = ax[keep_seg], y0 = ay[keep_seg],
    x1 = gl$text_x[keep_seg], y1 = gl$text_y[keep_seg],
    group = which(keep_seg),
    stringsAsFactors = FALSE
  )
}

# Split leader lines at the real oriented rectangles of other labels. Covered
# pieces can be faded, clipped, or shown in full. The target label itself is a
# hard boundary when include_own is TRUE; it is never represented by a faded
# line running through its own text.
