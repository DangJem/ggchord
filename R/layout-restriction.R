# Restriction labels use the same radial principles as gene-label repel, with
# local site clusters so a dense MCS may rise one level without moving every
# sparse callout. Most inner label edges therefore follow one common offset
# contour; only collision-bound clusters use a nearby outer contour.
ggchord_restriction_label_lanes <- function(
    gl, seq_arcs, side = "outside", units_per_inch = .3,
    label_offset = .13, min_label_gap = NULL) {
  frame <- ggchord_label_curve_frame(gl, seq_arcs)
  labels <- gl
  labels$.radial_bend_x <- NA_real_
  labels$.radial_bend_y <- NA_real_
  labels$.radial_parameter <- NA_real_
  directions <- rep(NA_character_, nrow(gl))
  tracks <- rep(NA_integer_, nrow(gl))
  provisional_boxes <- ggchord_text_boxes(
    labels, units_per_inch = units_per_inch
  )
  collision_boxes <- provisional_boxes
  # The final restriction label is drawn as separate bold enzyme and regular
  # coordinate grobs.  Reserve a modest composite-text margin during the
  # circular solve so the later exact metrics do not turn an apparently valid
  # main-track fan into many unnecessary radial spill tracks.
  provisional_boxes$w <- provisional_boxes$w * 1.18
  provisional_boxes$h <- provisional_boxes$h * 1.08
  arc_ids <- vapply(seq_arcs, function(arc) as.character(arc$accver[1L]),
    character(1))

  # Split the circle at genuine genomic gaps. Each local fan is packed near its
  # own sites, then neighbouring fans are reconciled without making a dense MCS
  # drag labels from a distant part of the plasmid.
  groups <- list()
  for (id in unique(as.character(gl$accver))) {
    rows <- which(as.character(gl$accver) == id)
    arc <- seq_arcs[[match(id, arc_ids)]]
    distance <- c(0, cumsum(sqrt(diff(arc$x)^2 + diff(arc$y)^2)))
    keep <- !duplicated(distance)
    arc <- arc[keep, , drop = FALSE]
    distance <- distance[keep]
    total <- tail(distance, 1L)
    preferred <- vapply(rows, function(i) distance[which.min(
      (arc$x - gl$anchor_x[i])^2 + (arc$y - gl$anchor_y[i])^2
    )], numeric(1))
    ord <- order(preferred, rows)
    if (length(ord) > 1L) {
      sorted <- preferred[ord]
      gaps <- c(diff(sorted), total - tail(sorted, 1L) + sorted[1L])
      cut <- which.max(gaps)
      ord <- c(if (cut < length(ord)) ord[seq.int(cut + 1L, length(ord))]
        else integer(), ord[seq_len(cut)])
    }
    opened <- preferred[ord]
    if (length(opened) > 1L) {
      for (j in 2:length(opened)) if (opened[j] < opened[j - 1L]) {
        opened[j:length(opened)] <- opened[j:length(opened)] + total
      }
    }
    # Neighbouring labels can collide even when their sites are several percent
    # of a plasmid apart, especially around 12/6 o'clock where horizontal text
    # occupies substantial tangential width. Solve such neighbours as one
    # ordered contour fan instead of allowing independently centred fans to
    # overlap.
    base_fan_gap <- max(.025, min_label_gap %||% .014) * total
    theta <- atan2(gl$anchor_y[rows[ord]], gl$anchor_x[rows[ord]])
    tangent_half_extent <-
      abs(sin(theta)) * provisional_boxes$w[rows[ord]] / 2 +
      abs(cos(theta)) * provisional_boxes$h[rows[ord]] / 2
    adaptive_gap <- if (length(opened) > 1L) {
      pmax(
        base_fan_gap,
        tangent_half_extent[-length(tangent_half_extent)] +
          tangent_half_extent[-1L] + .025
      )
    } else numeric()
    groups[[length(groups) + 1L]] <- rows[ord]
  }
  # Groups remain genomic/circular clusters.  Do not split them at Cartesian
  # top/right/bottom/left boundaries: the exterior resource is a circular
  # anchor track, not four independent text rails.

  # The final restriction layout has its own circular packing pass below.  It
  # does not need the gene-label solver's independent radial tracks; starting
  # from the real site angles avoids both needless outer levels and device-size
  # failures on otherwise drawable dense plasmid maps.
  labels$.radial_bend_x <- gl$anchor_x
  labels$.radial_bend_y <- gl$anchor_y
  labels$.radial_parameter <- if ("sequence_length" %in% names(gl)) {
    gl$genomic_position / gl$sequence_length
  } else rep(NA_real_, nrow(gl))
  directions <- ifelse(gl$anchor_x < 0, "left", "right")
  tracks[] <- 1L
  labels$text <- gl$text
  boxes <- provisional_boxes
  labels$.restriction_column <- FALSE
  labels$.restriction_column_edge <- NA_real_
  labels$.restriction_column_radius <- NA_real_
  labels$.restriction_column_side <- NA_character_
  labels$.restriction_column_group <- NA_integer_
  labels$.restriction_column_index <- NA_real_
  labels$.restriction_column_center_y <- NA_real_
  labels$.restriction_rail_side <- NA_character_
  labels$.restriction_rail_edge_x <- NA_real_
  labels$.restriction_rail_edge_y <- NA_real_

  # Retained fields keep the layout schema stable.  Four-sided Cartesian rails
  # were an intermediate experiment and are deliberately disabled: horizontal
  # text may extend beyond its circular anchor track, while its measured box is
  # used only for collision detection and local outward stacking.
  if (FALSE && !identical(side, "inside")) for (column_group in seq_along(groups)) {
    rows <- groups[[column_group]]
    if (length(rows) < 3L) next
    radius <- sqrt(gl$anchor_x[rows]^2 + gl$anchor_y[rows]^2)
    radius[radius <= 1e-10] <- 1
    x_score <- mean(abs(gl$anchor_x[rows]) / radius)
    y_score <- mean(abs(gl$anchor_y[rows]) / radius)
    vertical_cut <- .70
    vertical <- all(abs(gl$anchor_x[rows]) >= vertical_cut * radius)
    if (vertical && length(rows) < 5L) next
    cardinal <- if (vertical) x_score else y_score
    half_coordinate <- if (vertical) gl$anchor_x[rows] else gl$anchor_y[rows]
    one_half <- length(unique(sign(half_coordinate))) == 1L
    if (!one_half || cardinal < .55) next
    rail_side <- if (vertical) {
      if (mean(gl$anchor_x[rows]) >= 0) "right" else "left"
    } else if (mean(gl$anchor_y[rows]) >= 0) "top" else "bottom"

    genomic <- if ("genomic_position" %in% names(gl)) {
      gl$genomic_position[rows]
    } else seq_along(rows)
    ord <- order(genomic, rows)
    ordered <- rows[ord]
    axis_anchor <- if (vertical) gl$anchor_y[ordered] else
      gl$anchor_x[ordered]
    half_extent <- if (vertical) boxes$h[ordered] / 2 else
      boxes$w[ordered] / 2
    axis_position <- ggchord_pack_label_axis_outward(
      axis_anchor, axis_anchor, half_extent, half_extent,
      centre = 0,
      gap = if (vertical) .023 else .045
    )

    displacement <- abs(axis_position - axis_anchor)
    displacement_scale <- max(displacement)
    if (!is.finite(displacement_scale) || displacement_scale <= 1e-10) {
      displacement_scale <- 1
    }
    # Vary the radial departure with travel distance so the fan opens
    # gradually instead of creating one visible shared bend wall.
    stub <- .014 + .030 * sqrt(displacement / displacement_scale)
    stub <- pmax(.014, pmin(.046, stub))
    labels$.radial_bend_x[ordered] <- gl$anchor_x[ordered] +
      frame$outward_x[ordered] * stub
    labels$.radial_bend_y[ordered] <- gl$anchor_y[ordered] +
      frame$outward_y[ordered] * stub
    boundary <- max(radius, na.rm = TRUE) + max(.075, .55 * label_offset)
    if (rail_side == "right") {
      edge_x <- rep(boundary, length(ordered))
      edge_y <- axis_position
      labels$text_x[ordered] <- edge_x + boxes$w[ordered] / 2
      labels$text_y[ordered] <- edge_y
    } else if (rail_side == "left") {
      edge_x <- rep(-boundary, length(ordered))
      edge_y <- axis_position
      labels$text_x[ordered] <- edge_x - boxes$w[ordered] / 2
      labels$text_y[ordered] <- edge_y
    } else if (rail_side == "top") {
      edge_x <- axis_position
      edge_y <- rep(boundary, length(ordered))
      labels$text_x[ordered] <- edge_x
      labels$text_y[ordered] <- edge_y + boxes$h[ordered] / 2
    } else {
      edge_x <- axis_position
      edge_y <- rep(-boundary, length(ordered))
      labels$text_x[ordered] <- edge_x
      labels$text_y[ordered] <- edge_y - boxes$h[ordered] / 2
    }
    labels$.restriction_column[rows] <- TRUE
    labels$.restriction_column_edge[ordered] <- if (vertical) edge_x else edge_y
    labels$.restriction_column_radius[rows] <- boundary
    labels$.restriction_column_side[rows] <- rail_side
    labels$.restriction_column_group[rows] <- column_group
    labels$.restriction_column_index[ordered] <- seq_along(ordered)
    labels$.restriction_column_center_y[rows] <- mean(labels$text_y[ordered])
    labels$.restriction_rail_side[ordered] <- rail_side
    labels$.restriction_rail_edge_x[ordered] <- edge_x
    labels$.restriction_rail_edge_y[ordered] <- edge_y
    directions[rows] <- rail_side
    tracks[rows] <- 1L
  }
  # Adjacent genomic fans on the same physical side must share one measured
  # rail.  Repack them together after local classification; otherwise two
  # individually valid top fans can overlap at the origin seam.
  if (FALSE) for (id in unique(as.character(gl$accver))) {
    for (rail_side in c("left", "right", "top", "bottom")) {
      rows <- which(as.character(gl$accver) == id &
        labels$.restriction_rail_side %in% rail_side)
      if (length(rows) < 2L) next
      vertical <- rail_side %in% c("left", "right")
      anchor_axis <- if (vertical) gl$anchor_y[rows] else gl$anchor_x[rows]
      # The rail order follows the roots' projected physical order.  Around a
      # circular origin, sorting by the raw genomic coordinate can reverse one
      # side of the fan even though every local group is internally ordered.
      ordered <- rows[order(anchor_axis, gl$genomic_position[rows], rows)]
      half_extent <- if (vertical) boxes$h[ordered] / 2 else
        boxes$w[ordered] / 2
      axis_position <- ggchord_pack_label_axis_outward(
        anchor_axis[order(anchor_axis, gl$genomic_position[rows], rows)],
        anchor_axis[order(anchor_axis, gl$genomic_position[rows], rows)],
        half_extent, half_extent,
        centre = 0,
        gap = if (vertical) .023 else .045
      )
      boundary <- max(labels$.restriction_column_radius[ordered], na.rm = TRUE)
      if (rail_side == "right") {
        edge_x <- rep(boundary, length(ordered))
        edge_y <- axis_position
        labels$text_x[ordered] <- edge_x + boxes$w[ordered] / 2
        labels$text_y[ordered] <- edge_y
      } else if (rail_side == "left") {
        edge_x <- rep(-boundary, length(ordered))
        edge_y <- axis_position
        labels$text_x[ordered] <- edge_x - boxes$w[ordered] / 2
        labels$text_y[ordered] <- edge_y
      } else if (rail_side == "top") {
        edge_x <- axis_position
        edge_y <- rep(boundary, length(ordered))
        labels$text_x[ordered] <- edge_x
        labels$text_y[ordered] <- edge_y + boxes$h[ordered] / 2
      } else {
        edge_x <- axis_position
        edge_y <- rep(-boundary, length(ordered))
        labels$text_x[ordered] <- edge_x
        labels$text_y[ordered] <- edge_y - boxes$h[ordered] / 2
      }
      labels$.restriction_rail_edge_x[ordered] <- edge_x
      labels$.restriction_rail_edge_y[ordered] <- edge_y
      labels$.restriction_column_edge[ordered] <- if (vertical) edge_x else edge_y
      labels$.restriction_column_radius[ordered] <- boundary
    }
  }
  # Put every enzyme-facing text edge on one compact circular contour.  The
  # labels are packed only along that contour, in genomic order.  Horizontal
  # text then produces the characteristic stepped silhouette without a rigid
  # Cartesian column, and real sites remain encoded solely by leader roots.
  anchor_radius <- sqrt(gl$anchor_x^2 + gl$anchor_y^2)
  anchor_radius[anchor_radius <= 1e-10] <- 1
  labels$.restriction_contour <- FALSE
  labels$.restriction_contour_x <- NA_real_
  labels$.restriction_contour_y <- NA_real_
  labels$.restriction_contour_radius <- NA_real_
  labels$.restriction_base_contour_radius <- NA_real_
  labels$.restriction_fan_group <- NA_integer_
  labels$.restriction_fan_index <- NA_integer_
  labels$.restriction_fan_size <- NA_integer_
  labels$.restriction_orientation <- NA_real_
  labels$.restriction_text_sign <- NA_real_
  labels$.restriction_cut_spill <- FALSE
  radial_sign <- if (identical(side, "inside")) -1 else 1
  contour_radius_by_id <- vapply(unique(as.character(gl$accver)), function(id) {
    rows <- which(as.character(gl$accver) == id)
    base <- if (radial_sign > 0) max(anchor_radius[rows], na.rm = TRUE) else
      min(anchor_radius[rows], na.rm = TRUE)
    nominal <- base + radial_sign * max(.040, .30 * label_offset)
    if (radial_sign < 0 || length(rows) < 2L) return(nominal)
    # Horizontal labels can rotate from top/bottom into lateral sectors during
    # packing, where their full width becomes tangential. Size the dominant
    # track conservatively from that complete measured width; the bbox may
    # still extend beyond the track and remains the collision authority.
    anchor_theta <- atan2(gl$anchor_y[rows], gl$anchor_x[rows])
    tangential_width <- abs(sin(anchor_theta)) * collision_boxes$w[rows] +
      abs(cos(anchor_theta)) * collision_boxes$h[rows]
    required_arc <- sum(tangential_width + .010)
    measured <- required_arc / (2 * pi * .90)
    # Moderate maps retain a close dominant ring and let local conflicts
    # spill; very heavily annotated maps expand one common ring first. This
    # is a load class, not an identifier-specific layout branch.
    full_measured <- sum(collision_boxes$w[rows] + .010) /
      (2 * pi * .90)
    compact_radius <- min(base + .32, full_measured)
    anchor_angles <- (anchor_theta + 2 * pi) %% (2 * pi)
    quarter_fraction <- max(vapply(anchor_angles, function(angle) {
      sum(((anchor_angles - angle) %% (2 * pi)) < pi / 2)
    }, integer(1L))) / length(rows)
    dominant_radius <- if (full_measured < 2.48) compact_radius else if (
        quarter_fraction > .55) full_measured else measured
    max(nominal, dominant_radius)
  }, numeric(1))
  for (g in seq_along(groups)) {
    rows <- groups[[g]]
    if (!length(rows)) next
    id <- as.character(gl$accver[rows[1L]])
    ordered <- rows
    theta <- atan2(gl$anchor_y[ordered], gl$anchor_x[ordered])
    arc <- seq_arcs[[match(id, arc_ids)]]
    orientation <- attr(arc, "ggchord_path_direction") %||% -1
    if (length(ordered) > 1L) {
      delta <- atan2(sin(diff(theta)), cos(diff(theta)))
      observed <- sign(stats::median(delta[abs(delta) > 1e-10]))
      if (is.finite(observed) && observed != 0) orientation <- observed
    }
    total <- if ("sequence_length" %in% names(gl)) {
      unique(gl$sequence_length[ordered])
    } else numeric()
    if (length(total) != 1L || !is.finite(total) || total <= 0) {
      total <- max(gl$genomic_position, na.rm = TRUE)
    }
    contour_radius <- contour_radius_by_id[id]
    preferred_angle <- numeric(length(ordered))
    if (length(ordered) > 1L) {
      # `ordered` is opened in the physical arc order.  Derive its unwrapped
      # angular parameter from that physical order rather than from raw genomic
      # subtraction: after rotation or reversed path direction the latter can
      # jump at the origin and silently reverse a subset of packed endpoints.
      angular_step <- orientation * atan2(
        sin(diff(theta)), cos(diff(theta))
      )
      angular_step[angular_step < 0] <- angular_step[angular_step < 0] +
        2 * pi
      preferred_angle <- c(0, cumsum(angular_step))
    }

    dense_fan <- all(!is.na(labels$.restriction_rail_side[ordered]))
    if (dense_fan) {
      # Four-sided rails are already packed in measured text units.  Their
      # connection edge, rather than the text centre, is the leader target.
      edge_x <- labels$.restriction_rail_edge_x[ordered]
      edge_y <- labels$.restriction_rail_edge_y[ordered]
      contour_radius <- max(sqrt(edge_x^2 + edge_y^2))
      packed_theta <- atan2(edge_y, edge_x)
      packed <- orientation * (packed_theta - theta[1L]) * contour_radius
    } else {
      # The projected tangent footprint is the amount of circular arc occupied
      # by a horizontal text box. Iterate because it changes slightly as a
      # crowded label moves around the circle.
      packed <- preferred_angle * contour_radius
      for (iteration in seq_len(4L)) {
        packed_theta <- theta[1L] + orientation * packed / contour_radius
        half_extent <- abs(sin(packed_theta)) * boxes$w[ordered] / 2 +
          abs(cos(packed_theta)) * boxes$h[ordered] / 2
        packed <- ggchord_pack_label_axis(
          preferred_angle * contour_radius, seq_along(ordered),
          half_extent, half_extent, gap = .010
        )
      }
      crowding_ratio <- sum(collision_boxes$w[ordered]) /
        (2 * pi * contour_radius)
      if (crowding_ratio > 1.25) {
        preferred_distance <- preferred_angle * contour_radius
        max_tangent_shift <- .35 * contour_radius
        packed <- pmax(preferred_distance - max_tangent_shift,
          pmin(preferred_distance + max_tangent_shift, packed))
      }
      packed_theta <- theta[1L] + orientation * packed / contour_radius
      edge_x <- contour_radius * cos(packed_theta)
      edge_y <- contour_radius * sin(packed_theta)
    }
    text_side <- if (dense_fan) labels$.restriction_rail_side[ordered] else
      ifelse(edge_x < 0, "left", "right")
    text_sign <- radial_sign * ifelse(edge_x < 0, -1, 1)
    if (!dense_fan) {
      labels$text_x[ordered] <- edge_x +
        text_sign * boxes$w[ordered] / 2
      labels$text_y[ordered] <- edge_y
    }
    labels$.restriction_contour[ordered] <- TRUE
    labels$.restriction_contour_x[ordered] <- edge_x
    labels$.restriction_contour_y[ordered] <- edge_y
    labels$.restriction_contour_radius[ordered] <- sqrt(edge_x^2 + edge_y^2)
    labels$.restriction_base_contour_radius[ordered] <- contour_radius
    labels$.restriction_fan_group[ordered] <- g
    labels$.restriction_fan_index[ordered] <- seq_along(ordered)
    labels$.restriction_fan_size[ordered] <- length(ordered)
    labels$.restriction_orientation[ordered] <- orientation
    labels$.restriction_text_sign[ordered] <- text_sign
    dense <- labels$.restriction_column[ordered]
    vertical_dense <- dense & labels$.restriction_rail_side[ordered] %in%
      c("left", "right")
    labels$.restriction_column_edge[ordered[vertical_dense]] <-
      edge_x[vertical_dense]
    directions[ordered] <- text_side

    angular_displacement <- abs(packed / contour_radius - preferred_angle)
    displacement_scale <- max(angular_displacement)
    if (!is.finite(displacement_scale) || displacement_scale <= 1e-10) {
      displacement_scale <- 1
    }
    stub <- .012 + .022 * sqrt(angular_displacement / displacement_scale)
    stub <- pmax(.012, pmin(.034, stub))
    labels$.radial_bend_x[ordered] <- gl$anchor_x[ordered] +
      frame$outward_x[ordered] * stub
    labels$.radial_bend_y[ordered] <- gl$anchor_y[ordered] +
      frame$outward_y[ordered] * stub
  }

  # Reconcile neighbouring local fans on the same contour. A fan moves only
  # when its actual horizontal text boxes overlap a previously placed fan, and
  # it moves forward as one rigid unit. This preserves both the internal rhythm
  # and genomic order without letting one dense region drag distant labels.
  if (FALSE) for (id in unique(as.character(gl$accver))) {
    fan_ids <- which(vapply(groups, function(rows) {
      identical(as.character(gl$accver[rows[1L]]), id)
    }, logical(1)))
    placed <- integer()
    for (fan_id in fan_ids) {
      members <- which(labels$.restriction_fan_group == fan_id)
      if (length(members) && all(!is.na(
          labels$.restriction_rail_side[members]))) {
        placed <- c(placed, members)
        next
      }
      if (!length(placed)) {
        placed <- members
        next
      }
      orientation <- labels$.restriction_orientation[members[1L]]
      for (attempt in 0:0) {
        side_sign <- labels$.restriction_text_sign[members]
        centre_x <- labels$.restriction_contour_x[members] +
          side_sign * boxes$w[members] / 2
        current <- data.frame(
          xmin = centre_x - boxes$w[members] / 2 - .015,
          xmax = centre_x + boxes$w[members] / 2 + .015,
          ymin = labels$.restriction_contour_y[members] -
            boxes$h[members] / 2 - .015,
          ymax = labels$.restriction_contour_y[members] +
            boxes$h[members] / 2 + .015
        )
        placed_sign <- labels$.restriction_text_sign[placed]
        placed_cx <- labels$.restriction_contour_x[placed] +
          placed_sign * boxes$w[placed] / 2
        previous <- data.frame(
          xmin = placed_cx - boxes$w[placed] / 2 - .015,
          xmax = placed_cx + boxes$w[placed] / 2 + .015,
          ymin = labels$.restriction_contour_y[placed] -
            boxes$h[placed] / 2 - .015,
          ymax = labels$.restriction_contour_y[placed] +
            boxes$h[placed] / 2 + .015
        )
        overlap <- any(vapply(seq_len(nrow(current)), function(i) any(
          current$xmin[i] < previous$xmax &
            current$xmax[i] > previous$xmin &
            current$ymin[i] < previous$ymax &
            current$ymax[i] > previous$ymin
        ), logical(1)))
        if (!overlap) break
        theta <- atan2(labels$.restriction_contour_y[members],
          labels$.restriction_contour_x[members]) + orientation * .003
        radius <- labels$.restriction_contour_radius[members]
        labels$.restriction_contour_x[members] <- radius * cos(theta)
        labels$.restriction_contour_y[members] <- radius * sin(theta)
      }
      side_sign <- radial_sign * ifelse(
        labels$.restriction_contour_x[members] < 0, -1, 1
      )
      labels$.restriction_text_sign[members] <- side_sign
      labels$text_x[members] <- labels$.restriction_contour_x[members] +
        side_sign * boxes$w[members] / 2
      labels$text_y[members] <- labels$.restriction_contour_y[members]
      directions[members] <- ifelse(side_sign < 0, "left", "right")
      dense <- labels$.restriction_column[members]
      labels$.restriction_column_edge[members[dense]] <-
        labels$.restriction_contour_x[members[dense]]
      placed <- c(placed, members)
    }
  }

  # Reconcile the remaining boundary collisions with local exterior tracks.
  # Moving only the later genomic label radially outward retains each fan's
  # polar position and is much more compact than forcing all labels onto one
  # ever-growing global contour.
  # Restore the real horizontal text extent before collision/track assignment.
  # The enlarged provisional boxes above are packing slack, not visible ink.
  labels$text_x <- labels$.restriction_contour_x +
    labels$.restriction_text_sign * collision_boxes$w / 2
  labels$text_y <- labels$.restriction_contour_y
  boxes <- collision_boxes
  labels$.restriction_outer_track <- 1L
  shifted_box <- function(row, text_x, text_y,
                          padding = .045 * units_per_inch) {
    # Text metrics do not change when only the circular anchor moves. Reuse
    # the measured glyph box instead of opening a grid device for every
    # candidate track of every dense restriction label.
    box <- collision_boxes[row, , drop = FALSE]
    dx <- text_x - box$x
    dy <- text_y - box$y
    for (column in c("x", "cx", "xmin", "xmax")) {
      box[[column]] <- box[[column]] + dx
    }
    for (column in c("y", "cy", "ymin", "ymax")) {
      box[[column]] <- box[[column]] + dy
    }
    for (column in c("ow", "oh", "bw", "bh")) {
      box[[column]] <- box[[column]] + 2 * padding
    }
    box$xmin <- box$xmin - padding
    box$xmax <- box$xmax + padding
    box$ymin <- box$ymin - padding
    box$ymax <- box$ymax + padding
    box
  }
  for (id in unique(as.character(gl$accver))) {
    rows <- which(as.character(gl$accver) == id)
    if (length(rows) < 2L) next
    base_radius <- unique(labels$.restriction_contour_radius[rows])
    base_radius <- base_radius[is.finite(base_radius)][1L]
    crowding_ratio <- sum(collision_boxes$w[rows]) /
      (2 * pi * base_radius)
    if (!is.finite(crowding_ratio) || crowding_ratio <= 1.25) next
    ordered <- rows[order(gl$genomic_position[rows], rows)]
    placed <- data.frame()
    for (row in ordered) {
      if (!is.na(labels$.restriction_rail_side[row])) {
        candidate_box <- shifted_box(row, labels$text_x[row],
          labels$text_y[row])
        placed <- if (!nrow(placed)) candidate_box else
          rbind(placed, candidate_box)
        next
      }
      base_radius <- labels$.restriction_contour_radius[row]
      theta <- atan2(labels$.restriction_contour_y[row],
        labels$.restriction_contour_x[row])
      chosen <- NULL
      for (track in 1:12) {
        radius <- base_radius + (track - 1L) * max(.050, boxes$h[row] * .82)
        candidate_theta <- if (track > 1L &&
            labels$.restriction_fan_index[row] == 1L) {
          members <- which(labels$.restriction_fan_group ==
            labels$.restriction_fan_group[row])
          next_row <- members[which.min(abs(
            labels$.restriction_fan_index[members] - 2L
          ))]
          last_row <- members[which.max(
            labels$.restriction_fan_index[members]
          )]
          next_theta <- atan2(labels$.restriction_contour_y[next_row],
            labels$.restriction_contour_x[next_row])
          last_theta <- atan2(labels$.restriction_contour_y[last_row],
            labels$.restriction_contour_x[last_row])
          orientation <- labels$.restriction_orientation[row]
          gap_angle <- (orientation * (next_theta - last_theta)) %% (2 * pi)
          last_theta + orientation * gap_angle * .5
        } else theta
        edge_x <- radius * cos(candidate_theta)
        edge_y <- radius * sin(candidate_theta)
        text_sign <- radial_sign * ifelse(edge_x < 0, -1, 1)
        candidate_x <- edge_x + text_sign * boxes$w[row] / 2
        candidate_y <- edge_y
        candidate_box <- shifted_box(row, candidate_x, candidate_y)
        collision <- nrow(placed) &&
          any(ggchord_oriented_box_overlaps(candidate_box, placed))
        chosen <- list(radius = radius, edge_x = edge_x, edge_y = edge_y,
          text_sign = text_sign, text_x = candidate_x,
          text_y = candidate_y, box = candidate_box, track = track)
        if (!collision) break
      }
      labels$.restriction_contour_x[row] <- chosen$edge_x
      labels$.restriction_contour_y[row] <- chosen$edge_y
      labels$.restriction_contour_radius[row] <- chosen$radius
      labels$.restriction_text_sign[row] <- chosen$text_sign
      labels$text_x[row] <- chosen$text_x
      labels$text_y[row] <- chosen$text_y
      labels$.restriction_outer_track[row] <- chosen$track
      directions[row] <- if (chosen$text_sign < 0) "left" else "right"
      if (labels$.restriction_column[row]) {
        labels$.restriction_column_edge[row] <- chosen$edge_x
      }
      placed <- rbind(placed, chosen$box)
    }
  }
  boxes <- ggchord_text_boxes(labels, units_per_inch = units_per_inch)
  labels$.text_width <- boxes$w
  labels$.text_height <- boxes$h
  labels$.text_center_x <- boxes$cx
  labels$.text_center_y <- boxes$cy
  list(labels = labels, directions = directions, tracks = tracks,
    lanes = paste(gl$accver, side, sep = "\r"),
    draw_segment = !is.na(gl$text) & nzchar(gl$text))
}

ggchord_restriction_text_metrics <- function(text, size, units_per_inch,
                                             parse = TRUE,
                                             fontface = "plain",
                                             family = "", lineheight = 1.2) {
  n <- length(text)
  size <- rep_len(size, n)
  fontface <- rep_len(fontface, n)
  family <- rep_len(family, n)
  lineheight <- rep_len(lineheight, n)
  width <- height <- numeric(n)
  valid <- !is.na(text) & nzchar(text)
  if (!any(valid)) return(data.frame(width = width, height = height))
  close_device <- if (grDevices::dev.cur() > 1L) function() invisible(NULL) else
    ggchord_measurement_device(text = text[valid])
  on.exit(close_device())
  for (i in which(valid)) {
    label <- if (isTRUE(parse)) tryCatch(parse(text = text[i])[[1L]],
      error = function(e) text[i]) else text[i]
    grob <- grid::textGrob(
      label, gp = grid::gpar(
        fontsize = size[i] * (72.27 / 25.4), fontface = fontface[i],
        fontfamily = family[i], lineheight = lineheight[i]
      )
    )
    width[i] <- grid::convertWidth(
      grid::grobWidth(grob), "inches", valueOnly = TRUE
    ) * units_per_inch
    height[i] <- grid::convertHeight(
      grid::grobHeight(grob), "inches", valueOnly = TRUE
    ) * units_per_inch
  }
  data.frame(width = width, height = height)
}

ggchord_restriction_geometry <- function(data, params, layout, seq_data) {
  ggchord_require_columns(
    data, c("accver", "position", "enzyme"), "geom_restriction_site()"
  )
  lens <- stats::setNames(seq_data$length, seq_data$accver)
  if (any(!data$accver %in% names(lens))) {
    ggchord_stop("geom_restriction_site(): unknown accver")
  }
  if (!is.numeric(data$position) || any(!is.finite(data$position))) {
    ggchord_stop("geom_restriction_site(): position must contain finite numbers")
  }

  output <- list(); gid <- 0L
  append_path <- function(points, component, source_row = NA_integer_,
                          anchor_position = NA_real_, slot_position = NA_real_,
                          source_rows = integer(), cluster_id = NA_character_,
                          junction_id = NA_character_,
                          label_direction = NA_character_,
                          label_order = NA_character_,
                          label_connection_side = NA_character_) {
    gid <<- gid + 1L
    points$.component <- "path"
    points$restriction_component <- component
    points$group <- gid
    points$label <- NA_character_
    points$source_row <- source_row
    points$anchor_position <- anchor_position
    points$slot_position <- slot_position
    points$source_rows <- I(rep(list(as.integer(source_rows)), nrow(points)))
    points$cluster_id <- cluster_id
    points$junction_id <- junction_id
    points$label_direction <- label_direction
    points$label_order <- label_order
    points$label_connection_side <- label_connection_side
    output[[length(output) + 1L]] <<- points
  }

  for (id in unique(as.character(data$accver))) {
    idx <- which(as.character(data$accver) == id)
    pattern_order <- if ("pattern_id" %in% names(data)) {
      as.character(data$pattern_id[idx])
    } else rep("", length(idx))
    source_order <- if ("pattern_source_row" %in% names(data)) {
      data$pattern_source_row[idx]
    } else idx
    idx <- idx[order(data$position[idx], pattern_order, source_order, idx)]
    # Rendering collapses enzymes only when their displayed anchor coordinates
    # are exactly identical. The 17-digit key introduces no neighbourhood
    # tolerance: adjacent bp always remain separate callouts. Enzymes sharing
    # one anchor may have different second-strand cuts and still describe the
    # same displayed restriction site.
    all_idx <- idx
    site_key <- sprintf("%.17g", data$position[all_idx])
    site_members <- unname(split(idx,
      factor(site_key, levels = unique(site_key))))
    idx <- vapply(site_members, `[`, integer(1), 1L)
    site_enzyme <- vapply(site_members, function(rows) paste(
      unique(as.character(data$enzyme[rows])), collapse = " - "
    ), character(1))
    site_warning <- vapply(site_members, function(rows) {
      "display_warning" %in% names(data) &&
        any(data$display_warning[rows] %in% "methylation_blocked", na.rm = TRUE)
    }, logical(1L))
    site_enzyme[site_warning] <- paste0(site_enzyme[site_warning], " *")
    unique_cutters <- if ("enzyme_site_count" %in% names(data)) {
      unique(as.character(data$enzyme[all_idx])[
        data$enzyme_site_count[all_idx] == 1L])
    } else names(which(vapply(
      split(data$position[all_idx], as.character(data$enzyme[all_idx])),
      function(position) length(unique(position)) == 1L,
      logical(1L)
    )))
    site_enzyme_fontface <- vapply(site_members, function(rows) {
      enzymes <- unique(as.character(data$enzyme[rows]))
      if (length(enzymes) && all(enzymes %in% unique_cutters)) "bold" else "plain"
    }, character(1L))
    arc <- layout$seq_arcs[[id]]
    n <- nrow(arc)
    frac <- data$position[idx] / lens[id]

    point_at <- function(fraction, offset = 0) {
      fraction <- fraction %% 1
      continuous_index <- 1 + fraction * (n - 1)
      k <- pmax(1L, pmin(n - 1L, floor(continuous_index)))
      knext <- pmin(n, k + 1L)
      weight <- continuous_index - k
      base_x <- arc$x[k] + weight * (arc$x[knext] - arc$x[k])
      base_y <- arc$y[k] + weight * (arc$y[knext] - arc$y[k])
      kp <- pmax(1L, k - 1L)
      kn <- pmin(n, knext + 1L)
      tx <- arc$x[kn] - arc$x[kp]
      ty <- arc$y[kn] - arc$y[kp]
      tangent_length <- sqrt(tx^2 + ty^2)
      tangent_length[tangent_length <= 1e-12] <- 1
      nx <- -ty / tangent_length
      ny <- tx / tangent_length
      flip <- nx * base_x + ny * base_y < 0
      nx[flip] <- -nx[flip]
      ny[flip] <- -ny[flip]
      data.frame(x = base_x + offset * nx, y = base_y + offset * ny)
    }

    side_sign <- if (params$label_side == "inside") -1 else 1
    position_text <- format(data$position[idx], big.mark = "",
      scientific = FALSE, trim = TRUE)
    raw_labels <- switch(params$label_style,
      enzyme_position = paste0(site_enzyme, " (", position_text, ")"),
      enzyme = site_enzyme, position = position_text)
    # Cluster IDs remain useful export metadata, but they do not imply shared
    # geometry. Open at the largest circular gap so origin-adjacent sites keep
    # stable neighbouring IDs.
    cluster <- seq_along(idx)
    if (length(frac) > 1L) {
      gap_limit <- max(.025, params$min_label_gap %||% .014)
      opened <- frac
      if (isTRUE(layout$circular)) {
        gaps <- diff(c(frac, frac[1L] + 1))
        cut <- which.max(gaps)
        ord <- c(if (cut < length(frac)) seq.int(cut + 1L, length(frac)) else
          integer(), seq_len(cut))
        opened <- frac[ord]
        opened[opened < opened[1L]] <- opened[opened < opened[1L]] + 1
        opened_cluster <- cumsum(c(TRUE, diff(opened) > gap_limit))
        cluster[ord] <- opened_cluster
      } else {
        cluster <- cumsum(c(TRUE, diff(frac) > gap_limit))
      }
    }
    cluster_names <- paste0(id, ":cluster:", sprintf("%03d", cluster))
    junction_names <- paste0(id, ":site:", sprintf("%06d", idx))
    # Restriction leaders attach to the visible outer boundary of the sequence
    # backbone, rather than disappearing underneath a double line or band.
    bases <- point_at(frac, params$backbone_outer_offset %||% 0)
    tips <- point_at(frac, side_sign * params$tick_length)

    for (member in seq_along(idx)) {
      row <- idx[member]
      append_path(
        rbind(bases[member, , drop = FALSE], tips[member, , drop = FALSE]),
        "tick", row, data$position[row], data$position[row],
        site_members[[member]],
        cluster_names[member], junction_names[member]
      )
    }
    if (!isTRUE(params$label)) next

    gl <- data.frame(
      text = raw_labels, text_x = tips$x, text_y = tips$y,
      text_angle = 0, hjust = .5, vjust = .5,
      size = params$label_size %||% 2.9, accver = id,
      group = seq_along(idx), source_row = idx,
      genomic_position = data$position[idx],
      sequence_length = unname(lens[id]),
      anchor_x = tips$x, anchor_y = tips$y,
      stringsAsFactors = FALSE
    )
    device <- if (grDevices::dev.cur() == 1L) c(8, 8) else
      grDevices::dev.size("in")
    span <- c(diff(range(arc$x)), diff(range(arc$y)))
    units_per_inch <- max(span / pmax(device - c(2.2, 1.5), device * .45))
    units_per_inch <- max(.20, min(.42, units_per_inch))
    solved <- ggchord_restriction_label_lanes(
      gl, layout$seq_arcs, side = params$label_side,
      units_per_inch = units_per_inch,
      label_offset = params$label_offset %||% .13,
      min_label_gap = params$min_label_gap
    )
    labels <- solved$labels
    directions <- solved$directions
    tracks <- solved$tracks
    bend_x <- labels$.radial_bend_x
    bend_y <- labels$.radial_bend_y
    invalid_bend <- !is.finite(bend_x) | !is.finite(bend_y)
    bend_x[invalid_bend] <- tips$x[invalid_bend]
    bend_y[invalid_bend] <- tips$y[invalid_bend]
    centre_x <- labels$.text_center_x
    centre_y <- labels$.text_center_y
    # Text order and attachment edge follow the label's visual half of the
    # plasmid, not merely the centre-to-bend vector.  The latter is unstable
    # near twelve and six o'clock because a wide label can straddle its bend.
    visual_x <- centre_x
    visual_x[abs(visual_x) < .010] <- gl$anchor_x[abs(visual_x) < .010]
    connection_side <- ifelse(visual_x < 0, "right", "left")
    contour_pre <- if (".restriction_contour" %in% names(labels)) {
      labels$.restriction_contour %in% TRUE
    } else rep(FALSE, length(idx))
    rail <- if (".restriction_rail_side" %in% names(labels)) {
      !is.na(labels$.restriction_rail_side)
    } else rep(FALSE, length(idx))
    connection_side[contour_pre] <- ifelse(
      labels$.restriction_text_sign[contour_pre] > 0, "left", "right"
    )
    rail_side <- labels$.restriction_rail_side
    connection_side[rail] <- ifelse(
      rail_side[rail] == "right", "left",
      ifelse(rail_side[rail] == "left", "right", rail_side[rail])
    )
    order_values <- rep(params$label_order, length(idx))
    if (params$label_order == "auto") {
      order_values <- ifelse(connection_side == "right",
        "position_enzyme", "enzyme_position")
      horizontal_rail <- rail & connection_side %in% c("top", "bottom")
      order_values[horizontal_rail] <- ifelse(
        centre_x[horizontal_rail] < 0,
        "position_enzyme", "enzyme_position"
      )
    }
    quote_text <- function(x) paste0("'", gsub("'", "\\\\'", x,
      fixed = TRUE), "'")
    coordinate <- format(data$position[idx], big.mark = "",
      scientific = FALSE, trim = TRUE)
    label_text <- vapply(seq_along(idx), function(member) {
      switch(params$label_style,
        enzyme_position = if (order_values[member] == "position_enzyme")
          paste0("(", coordinate[member], ") ", site_enzyme[member]) else
          paste0(site_enzyme[member], " (", coordinate[member], ")"),
        enzyme = site_enzyme[member], position = coordinate[member])
    }, character(1))
    plotmath <- vapply(seq_along(idx), function(member) {
      enzyme_expression <- if (site_enzyme_fontface[member] == "bold") {
        paste0("bold(", quote_text(site_enzyme[member]), ")")
      } else quote_text(site_enzyme[member])
      switch(params$label_style,
        enzyme_position = if (order_values[member] == "position_enzyme")
          paste0(quote_text(paste0("(", coordinate[member], ")")),
            "~", enzyme_expression) else
          paste0(enzyme_expression, "~",
            quote_text(paste0("(", coordinate[member], ")"))),
        enzyme = enzyme_expression,
        position = quote_text(coordinate[member]))
    }, character(1))
    enzyme_label <- if (params$label_style %in% c("enzyme_position", "enzyme"))
      site_enzyme else rep("", length(idx))
    coordinate_label <- if (params$label_style == "enzyme_position") {
      paste0("(", coordinate, ")")
    } else if (params$label_style == "position") coordinate else
      rep("", length(idx))
    # Use the renderer's fragments and physical gap, including bold names.
    metrics <- if (isTRUE(params$composite_labels %||% TRUE)) {
      ggchord_site_label_metrics(enzyme_label, coordinate_label,
        params$label_size %||% 2.9, params$label_family %||% '',
        site_enzyme_fontface, params$label_lineheight %||% 1.2, units_per_inch)
    } else ggchord_restriction_text_metrics(label_text,
      params$label_size %||% 2.9, units_per_inch, parse = FALSE,
      fontface = params$label_fontface %||% 'plain', family = params$label_family %||% '',
      lineheight = params$label_lineheight %||% 1.2)
    column <- if (".restriction_column" %in% names(labels)) {
      labels$.restriction_column %in% TRUE
    } else rep(FALSE, length(idx))
    contour <- if (".restriction_contour" %in% names(labels)) {
      labels$.restriction_contour %in% TRUE
    } else rep(FALSE, length(idx))
    if (any(contour)) {
      ordinary_contour <- contour & !rail
      positive_text <- ordinary_contour & labels$.restriction_text_sign > 0
      negative_text <- ordinary_contour & labels$.restriction_text_sign < 0
      connection_side[positive_text] <- "left"
      connection_side[negative_text] <- "right"
      centre_x[ordinary_contour] <-
        labels$.restriction_contour_x[ordinary_contour] +
        labels$.restriction_text_sign[ordinary_contour] *
          metrics$width[ordinary_contour] / 2
      centre_y[ordinary_contour] <-
        labels$.restriction_contour_y[ordinary_contour]
      right_rail <- rail & rail_side == "right"
      left_rail <- rail & rail_side == "left"
      top_rail <- rail & rail_side == "top"
      bottom_rail <- rail & rail_side == "bottom"
      centre_x[right_rail] <- labels$.restriction_rail_edge_x[right_rail] +
        metrics$width[right_rail] / 2
      centre_x[left_rail] <- labels$.restriction_rail_edge_x[left_rail] -
        metrics$width[left_rail] / 2
      centre_y[c(top_rail | bottom_rail)] <-
        labels$.restriction_rail_edge_y[c(top_rail | bottom_rail)] +
        ifelse(top_rail[c(top_rail | bottom_rail)], 1, -1) *
          metrics$height[c(top_rail | bottom_rail)] / 2
      centre_x[top_rail | bottom_rail] <-
        labels$.restriction_rail_edge_x[top_rail | bottom_rail]
      centre_y[right_rail | left_rail] <-
        labels$.restriction_rail_edge_y[right_rail | left_rail]
    }
    routing_half_width <- metrics$width / 2 + .003
    routing_half_height <- metrics$height / 2 + .003
    lateral_connection <- connection_side %in% c("left", "right")
    endpoint_x <- pmax(
      centre_x - routing_half_width,
      pmin(bend_x, centre_x + routing_half_width)
    )
    endpoint_x[lateral_connection] <- centre_x[lateral_connection] + ifelse(
      connection_side[lateral_connection] == "left",
      -routing_half_width[lateral_connection],
      routing_half_width[lateral_connection]
    )
    # A selected left edge must lie to the right of its bend, and a selected
    # right edge to the left.  Bottom/top packing can otherwise leave a wide
    # bbox straddling the bend, so the nominally correct edge is approached
    # backwards. Move only as far as needed to restore the actual approach.
    route_clearance <- .012
    shift_x <- ifelse(
      connection_side == "left",
      pmax(0, bend_x + route_clearance - endpoint_x),
      pmin(0, bend_x - route_clearance - endpoint_x)
    )
    shift_x[!lateral_connection] <- 0
    # The shared contour is the layout invariant. Do not move a label off that
    # contour merely to regularise the final approach angle; the selected edge
    # is already the enzyme-facing edge and remains the correct attachment.
    shift_x[contour] <- 0
    centre_x <- centre_x + shift_x
    endpoint_x <- endpoint_x + shift_x
    labels$.text_center_x <- centre_x
    labels$text_x <- labels$text_x + shift_x
    endpoint_y <- pmax(
      centre_y - routing_half_height,
      pmin(bend_y, centre_y + routing_half_height)
    )
    # Horizontal exterior labels attach at the midpoint of their circle-facing
    # edge. Clamping that endpoint to an old bend y-coordinate can exchange the
    # order of two adjacent labels after one spills outward, creating a tiny
    # but real crossing at the dominant track boundary.
    endpoint_y[lateral_connection] <- centre_y[lateral_connection]
    endpoint_y[connection_side == "top"] <-
      centre_y[connection_side == "top"] -
      routing_half_height[connection_side == "top"]
    endpoint_y[connection_side == "bottom"] <-
      centre_y[connection_side == "bottom"] +
      routing_half_height[connection_side == "bottom"]
    endpoint_x[rail] <- labels$.restriction_rail_edge_x[rail]
    endpoint_y[rail] <- labels$.restriction_rail_edge_y[rail]
    # Keep the rendered label centred on the box used by the routing solver.
    # The leader endpoint is a separate bbox intersection, so justification no
    # longer shifts the visible edge after routing has been computed.
    labels$hjust <- .5
    labels$vjust <- .5
    labels$text_x <- centre_x
    labels$text_y <- centre_y
    slot_positions <- if (".radial_parameter" %in% names(labels)) {
      labels$.radial_parameter * lens[id]
    } else vapply(seq_along(idx), function(member) {
      k <- which.min((arc$x - labels$text_x[member])^2 +
        (arc$y - labels$text_y[member])^2)
      (k - 1L) / max(1L, nrow(arc) - 1L) * lens[id]
    }, numeric(1))

    leader_frame <- ggchord_label_curve_frame(gl, layout$seq_arcs)
    point_segment_radius <- function(p0, p1) {
      delta <- p1 - p0
      denom <- sum(delta^2)
      parameter <- if (denom <= 1e-16) 0 else
        max(0, min(1, -sum(p0 * delta) / denom))
      sqrt(sum((p0 + parameter * delta)^2))
    }
    fan_size <- tabulate(labels$.restriction_fan_group)
    ordered_fan <- fan_size[labels$.restriction_fan_group] > 1L
    leader_paths <- lapply(seq_along(idx), function(member) {
      vector <- c(endpoint_x[member] - tips$x[member],
        endpoint_y[member] - tips$y[member])
      radial <- c(leader_frame$outward_x[member] * side_sign,
        leader_frame$outward_y[member] * side_sign)
      denominator <- sqrt(sum(vector^2) * sum(radial^2))
      centre_vector <- c(labels$.text_center_x[member] - tips$x[member],
        labels$.text_center_y[member] - tips$y[member])
      centre_length <- sqrt(sum(centre_vector^2))
      aligned <- is.finite(denominator) && denominator > 1e-12 &&
        is.finite(centre_length) && centre_length > 1e-12 &&
        sum(centre_vector * radial) > 0 &&
        abs(centre_vector[1] * radial[2] - centre_vector[2] * radial[1]) /
          centre_length < .02
      start <- c(tips$x[member], tips$y[member])
      endpoint <- c(endpoint_x[member], endpoint_y[member])
      if (identical(params$leader, "straight")) {
        return(data.frame(x = c(start[1], endpoint[1]),
          y = c(start[2], endpoint[2])))
      }
      if (identical(params$leader, "elbow") && side_sign < 0) {
        # Interior callouts must not enter the exterior-circle clearance
        # search: their tick roots are already inside that forbidden radius.
        shoulder <- min(.085, sqrt(sum(vector^2)) * .30)
        bend <- start + radial * shoulder
        return(data.frame(x = c(start[1], bend[1], endpoint[1]),
          y = c(start[2], bend[2], endpoint[2])))
      }
      # The leader may pass inside its own tick tip, but never inside the
      # visible sequence boundary where that tick begins.  Treating the tip
      # radius as the forbidden circle forces harmless chords into elbows;
      # those elbows can then cross an adjacent, correctly ordered leader.
      base <- unlist(bases[member, c("x", "y")])
      safe_radius <- sqrt(sum(base^2)) + 1e-6
      direct_is_safe <- point_segment_radius(start, endpoint) >= safe_radius
      # The shortest exterior segment is also the only order-preserving route
      # shared by circular fans and Cartesian rails.  Adding independently
      # sized radial stubs before an otherwise safe chord can make adjacent,
      # correctly ordered leaders cross one another.  Keep a bend only when a
      # direct chord would enter the plasmid interior.
      if (!identical(params$leader, "elbow") &&
          direct_is_safe && !ordered_fan[member]) {
        return(data.frame(x = c(start[1], endpoint[1]),
          y = c(start[2], endpoint[2])))
      }
      if (!identical(params$leader, "elbow") && ordered_fan[member]) {
        endpoint_angle <- atan2(endpoint[2], endpoint[1])
        main_radius <- labels$.restriction_base_contour_radius[member]
        if (!is.finite(main_radius)) main_radius <- sqrt(sum(endpoint^2))
        elbow <- main_radius * c(cos(endpoint_angle), sin(endpoint_angle))
        if (point_segment_radius(start, elbow) >= safe_radius &&
            point_segment_radius(elbow, endpoint) >= safe_radius) {
          return(data.frame(x = c(start[1], elbow[1], endpoint[1]),
            y = c(start[2], elbow[2], endpoint[2])))
        }
      }
      bend <- c(labels$.radial_bend_x[member],
        labels$.radial_bend_y[member])
      if (any(!is.finite(bend)) ||
          sqrt(sum((bend - unlist(tips[member, ]))^2)) < 1e-8) {
        distance <- sqrt(sum(vector^2))
        shoulder <- min(.085, max(.045, distance * .30))
        bend <- c(tips$x[member] + radial[1] * shoulder,
          tips$y[member] + radial[2] * shoulder)
      }
      if (point_segment_radius(start, bend) >= safe_radius &&
          point_segment_radius(bend, endpoint) >= safe_radius) {
        return(data.frame(x = c(start[1], bend[1], endpoint[1]),
          y = c(start[2], bend[2], endpoint[2])))
      }
      # Never turn an exterior site leader into a sampled circular arc.  When
      # the direct chord would enter the sequence circle, move one elbow out
      # along the root radius until the second straight segment clears it.
      # This preserves the reference convention: a straight leader whenever
      # possible, otherwise a simple two-segment polyline.
      start_radius <- sqrt(sum(start^2))
      endpoint_radius <- sqrt(sum(endpoint^2))
      root_direction <- if (start_radius > 1e-12) start / start_radius else radial
      search_limit <- max(endpoint_radius,
        labels$.restriction_base_contour_radius[member], safe_radius) + .60
      search_radii <- seq(max(start_radius + .012, safe_radius + .012),
        search_limit, length.out = 80L)
      candidates <- lapply(search_radii, function(radius) root_direction * radius)
      valid <- vapply(candidates, function(candidate) {
        point_segment_radius(start, candidate) >= safe_radius &&
          point_segment_radius(candidate, endpoint) >= safe_radius
      }, logical(1))
      if (any(valid)) {
        elbow <- candidates[[which(valid)[1L]]]
        return(data.frame(x = c(start[1], elbow[1], endpoint[1]),
          y = c(start[2], elbow[2], endpoint[2])))
      }

      # Numerical last resort: retain a single, visibly angular elbow.  The
      # larger radius keeps both segments outside the sequence boundary even
      # for an unusually displaced user-supplied label.
      start_angle <- atan2(start[2], start[1])
      end_angle <- atan2(endpoint[2], endpoint[1])
      delta_angle <- atan2(sin(end_angle - start_angle),
        cos(end_angle - start_angle))
      elbow_angle <- start_angle + delta_angle / 2
      elbow_radius <- search_limit + .60
      elbow <- elbow_radius * c(cos(elbow_angle), sin(elbow_angle))
      data.frame(x = c(start[1], elbow[1], endpoint[1]),
        y = c(start[2], elbow[2], endpoint[2]))
    })

    # Resolve the rare intersections that remain between straight chords near
    # a dense root cluster.  Search only angular one-elbow polylines; sampled
    # arcs are intentionally excluded.  Earlier accepted paths are fixed, so
    # this pass is deterministic and preserves the label solver's order.
    path_segments <- function(path) {
      if (!is.data.frame(path) || nrow(path) < 2L) return(data.frame())
      data.frame(
        x1 = path$x[-nrow(path)], y1 = path$y[-nrow(path)],
        x2 = path$x[-1L], y2 = path$y[-1L]
      )
    }
    segments_cross <- function(a, b, tolerance = 1e-9) {
      cross2 <- function(x, y) x[1L] * y[2L] - x[2L] * y[1L]
      p <- c(a$x1, a$y1)
      r <- c(a$x2 - a$x1, a$y2 - a$y1)
      q <- c(b$x1, b$y1)
      s <- c(b$x2 - b$x1, b$y2 - b$y1)
      denominator <- cross2(r, s)
      if (!is.finite(denominator) || abs(denominator) <= tolerance) return(FALSE)
      t <- cross2(q - p, s) / denominator
      u <- cross2(q - p, r) / denominator
      t > tolerance && t < 1 - tolerance &&
        u > tolerance && u < 1 - tolerance
    }
    path_crosses <- function(path, accepted) {
      candidate_segments <- path_segments(path)
      if (!nrow(candidate_segments) || !length(accepted)) return(FALSE)
      any(vapply(accepted, function(other) {
        other_segments <- path_segments(other)
        any(vapply(seq_len(nrow(candidate_segments)), function(i) {
          any(vapply(seq_len(nrow(other_segments)), function(j) {
            segments_cross(candidate_segments[i, ], other_segments[j, ])
          }, logical(1L)))
        }, logical(1L)))
      }, logical(1L)))
    }
    accepted <- list()
    for (member in seq_along(leader_paths)) {
      current <- leader_paths[[member]]
      if (!identical(params$leader, "straight") &&
          !(identical(params$leader, "elbow") && side_sign < 0) &&
          path_crosses(current, accepted)) {
        start <- c(tips$x[member], tips$y[member])
        endpoint <- c(endpoint_x[member], endpoint_y[member])
        base <- unlist(bases[member, c("x", "y")])
        safe_radius <- sqrt(sum(base^2)) + 1e-6
        start_angle <- atan2(start[2], start[1])
        endpoint_angle <- atan2(endpoint[2], endpoint[1])
        delta_angle <- atan2(sin(endpoint_angle - start_angle),
          cos(endpoint_angle - start_angle))
        main_radius <- labels$.restriction_base_contour_radius[member]
        if (!is.finite(main_radius)) main_radius <- sqrt(sum(endpoint^2))
        outer_radius <- max(main_radius, sqrt(sum(endpoint^2)))
        candidate_angles <- start_angle + delta_angle *
          c(1, .875, .75, .625, .5, .375, .25, .125, 0)
        candidate_radii <- outer_radius + c(0, .025, .05, .085, .13,
          .19, .27, .38, .52)
        alternatives <- list()
        for (radius in candidate_radii) for (angle in candidate_angles) {
          elbow <- radius * c(cos(angle), sin(angle))
          if (point_segment_radius(start, elbow) < safe_radius ||
              point_segment_radius(elbow, endpoint) < safe_radius) next
          alternative <- data.frame(
            x = c(start[1], elbow[1], endpoint[1]),
            y = c(start[2], elbow[2], endpoint[2])
          )
          if (!path_crosses(alternative, accepted)) {
            alternatives[[length(alternatives) + 1L]] <- alternative
          }
        }
        if (length(alternatives)) {
          lengths <- vapply(alternatives, function(path) sum(sqrt(
            diff(path$x)^2 + diff(path$y)^2
          )), numeric(1L))
          current <- alternatives[[which.min(lengths)]]
          leader_paths[[member]] <- current
        }
      }
      accepted[[length(accepted) + 1L]] <- current
    }
    for (member in seq_along(leader_paths)) {
      row <- idx[member]
      append_path(
        leader_paths[[member]],
        "leader", row, data$position[row], slot_positions[member],
        site_members[[member]],
        cluster_names[member], junction_names[member], directions[member],
        order_values[member], connection_side[member]
      )
    }

    for (member in seq_along(idx)) {
      row <- idx[member]
      direction <- directions[member]
      order_value <- order_values[member]
      gid <- gid + 1L
      output[[length(output) + 1L]] <- data.frame(
        x = labels$text_x[member], y = labels$text_y[member],
        label = label_text[member], plotmath_label = plotmath[member],
        .component = "label", restriction_component = "label",
        group = gid, source_row = row,
        anchor_position = data$position[row],
        slot_position = slot_positions[member],
        source_rows = I(list(as.integer(site_members[[member]]))),
        cluster_id = cluster_names[member],
        junction_id = junction_names[member],
        label_direction = direction, label_order = order_value,
        label_connection_side = connection_side[member],
        label_layout = if (rail[member]) "perimeter_rail" else if (
          column[member]) "contour_fan" else "contour",
        label_boundary = if (column[member])
          labels$.restriction_column_edge[member] else NA_real_,
        label_attachment_x = labels$.restriction_contour_x[member],
        label_attachment_y = labels$.restriction_contour_y[member],
        label_contour_radius = labels$.restriction_contour_radius[member],
        outer_annotation_region = "restriction",
        outer_track = labels$.restriction_outer_track[member],
        outer_slot = 0L,
        enzyme_label = enzyme_label[member],
        enzyme_fontface = site_enzyme_fontface[member],
        coordinate_label = coordinate_label[member],
        .ggchord_composite_label = isTRUE(params$composite_labels %||% TRUE),
        family = params$label_family %||% '', fontface = params$label_fontface %||% 'plain',
        lineheight = params$label_lineheight %||% 1.2,
        size = params$label_size %||% 2.9, angle = 0,
        hjust = labels$hjust[member], vjust = labels$vjust[member],
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(output)) {
    return(data.frame(
      x = numeric(), y = numeric(), label = character(),
      .component = character(), restriction_component = character(),
      group = integer(), source_row = integer(),
      anchor_position = numeric(), slot_position = numeric(),
      source_rows = I(list()), cluster_id = character(),
      junction_id = character(), label_direction = character(),
      label_order = character(), label_connection_side = character(),
      label_layout = character(), label_boundary = numeric(),
      label_attachment_x = numeric(), label_attachment_y = numeric(),
      label_contour_radius = numeric(),
      outer_annotation_region = character(), outer_track = integer(),
      outer_slot = integer(),
      stringsAsFactors = FALSE
    ))
  }
  rendered <- ggchord_rbind_fill(output)
  rendered$.leader_style <- params$leader %||% "radial"
  if ("display_warning" %in% names(data) && nrow(rendered)) {
    blocked <- vapply(rendered$source_rows, function(rows) {
      rows <- as.integer(rows)
      rows <- rows[is.finite(rows) & rows >= 1L & rows <= nrow(data)]
      length(rows) && any(data$display_warning[rows] %in%
        "methylation_blocked", na.rm = TRUE)
    }, logical(1L))
    rendered$colour <- ifelse(blocked, "#858585", "#202020")
  }
  rendered
}
