# Coord-owned perimeter annotation coordination for circular maps.

ggchord_external_candidates <- function(max_band = 8L, max_slot = 6L,
                                         class = "restriction") {
  candidates <- expand.grid(
    band = seq_len(max_band),
    slot = c(0L, as.vector(rbind(seq_len(max_slot), -seq_len(max_slot))))
  )
  # A shared dominant band is the first objective. Small local tangential
  # movements are cheaper than opening a new band, but remain bounded so an
  # annotation cannot migrate around the circle.
  radial_cost <- if (identical(class, "restriction")) 3.4 else 5.2
  candidates$score <- (candidates$band - 1L) * radial_cost +
    abs(candidates$slot) * .72
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
                                    class = "restriction") {
  x <- row$text_x %||% row$x
  y <- row$text_y %||% row$y
  radius <- sqrt(x^2 + y^2)
  if (!is.finite(radius) || radius <= 1e-8) {
    return(list(row = row, box = data.frame(), band = NA_integer_,
      slot = NA_integer_, collision = FALSE))
  }
  direction <- c(x, y) / radius
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
  radial_half <- abs(direction[1L]) * measured$w[1L] / 2 +
    abs(direction[2L]) * measured$h[1L] / 2
  if (is.null(base_edge_radius) || !is.finite(base_edge_radius)) {
    base_edge_radius <- max(0, radius - radial_half)
  }
  step <- max(.052, measured$h[1L] * 1.02)
  slot_step <- max(.038, measured$h[1L] * 1.05)
  candidates <- ggchord_external_candidates(class = class)
  chosen <- NULL
  for (i in seq_len(nrow(candidates))) {
    band <- candidates$band[i]
    slot <- candidates$slot[i]
    candidate_radius <- base_edge_radius + radial_half + (band - 1L) * step
    candidate <- row
    candidate$text_x <- direction[1L] * candidate_radius +
      tangent[1L] * slot * slot_step
    candidate$text_y <- direction[2L] * candidate_radius +
      tangent[2L] * slot * slot_step
    candidate$x <- candidate$text_x
    candidate$y <- candidate$text_y
    box <- ggchord_text_boxes(
      candidate, units_per_inch = units_per_inch, box_padding = .045
    )
    collision <- nrow(occupied) &&
      any(ggchord_oriented_box_overlaps(box, occupied))
    chosen <- list(row = candidate, box = box, band = band, slot = slot,
      collision = collision, radial_half = radial_half)
    if (!collision) break
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
  geometry$outer_spill_reason[row] <- if (candidate$band > 1L) {
    "bbox_collision"
  } else NA_character_
  geometry$leader_crossing_count[row] <- 0L
  segment_rows <- which(geometry$.component %in% "segment" &
    geometry$source_row == geometry$source_row[row])
  if (length(segment_rows)) {
    geometry$xend[segment_rows] <- geometry$xend[segment_rows] + dx
    geometry$yend[segment_rows] <- geometry$yend[segment_rows] + dy
    if ("x1" %in% names(geometry)) {
      geometry$x1[segment_rows] <- geometry$x1[segment_rows] + dx
      geometry$y1[segment_rows] <- geometry$y1[segment_rows] + dy
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
      rows <- which(geometry$.component %in% "segment" &
        is.finite(geometry$source_row))
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
        radius <- sqrt(sites$x[index]^2 + sites$y[index]^2)
        index <- index[order(radius, seq_along(index))]
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
# Restriction sites already own a four-sided, order-preserving contour solver
# in restriction-sites.R.  That geometry is the fixed perimeter scaffold:
# flattening it back onto polar bands here recreates the very "flower" layout
# the contour solver was designed to avoid.  Feature and primer callouts are
# therefore placed against the measured restriction boxes and against one
# another, while the restriction rail itself remains unchanged.
ggchord_share_external_annotations <- function(registry, layout) {
  if (!isTRUE(layout$circular) || !length(registry)) return(registry)
  units_per_inch <- layout$text_units_per_inch %||% .30
  occupied <- data.frame()

  # Register the already solved restriction contours first.  Their exact
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
    registry[[id]]$restriction_site <- sites
  }

  # Place feature and primer labels in stable circular order.  Candidate bands
  # are local to the callout rail, but every candidate sees the fixed
  # restriction contour through `occupied`.
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
    angle <- atan2(geometry$text_y[rows], geometry$text_x[rows])
    rows <- rows[ggchord_open_circular_order(
      angle, geometry$anchor_position[rows] %||% rows
    )]
    for (row in rows) {
      class <- if ("annotation_class" %in% names(geometry) &&
        identical(as.character(geometry$annotation_class[row]), "primer")) {
        "primer"
      } else "feature"
      candidate <- ggchord_external_place(
        geometry[row, , drop = FALSE], occupied, units_per_inch,
        class = class
      )
      geometry <- ggchord_shift_feature_callout(geometry, row, candidate)
      occupied <- if (!nrow(occupied)) candidate$box else
        ggchord_rbind_fill(list(occupied, candidate$box))
    }
    used_bands <- sort(unique(geometry$outer_track[rows]))
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
