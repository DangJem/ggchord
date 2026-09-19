# Perimeter annotation coordination for circular maps.

# Open a horizontal fan corridor only for genuinely dense restriction
# perimeters.  Sparse and medium maps keep their natural polar contour.  For a
# dense map, every leader retains its backbone root and genomic order while its
# outer points progressively follow the label into wider left/right space.
ggchord_expand_dense_restriction_fans <- function(registry) {
  label_count <- sum(vapply(registry, function(owner) {
    sites <- owner$restriction_site
    if (!is.data.frame(sites) || !nrow(sites)) return(0L)
    as.integer(sum(sites$.component %in% "label" & !is.na(sites$label)))
  }, integer(1L)))
  if (label_count <= 53L) return(registry)

  horizontal_gain <- min(1.25, .90 + .08 * (label_count - 54L))
  for (id in names(registry)) {
    sites <- registry[[id]]$restriction_site
    if (!is.data.frame(sites) || !nrow(sites)) next
    label_rows <- which(sites$.component %in% "label" & !is.na(sites$label))
    if (!length(label_rows)) next
    for (row in label_rows) {
      original_x <- sites$x[row]
      original_y <- sites$y[row]
      radius <- sqrt(original_x^2 + original_y^2)
      if (!is.finite(radius) || radius <= 0) next
      # Top/bottom labels receive less horizontal motion than side labels, so
      # the natural polar region remains recognizable.
      side_weight <- abs(original_x) / radius
      dx <- original_x * horizontal_gain * (.35 + .65 * side_weight)
      if (!is.finite(dx) || abs(dx) < 1e-8) next
      source <- sites$source_row[row]
      leader_rows <- which(
        !is.na(sites$source_row) & sites$source_row == source &
          sites$restriction_component %in% "leader"
      )
      if (length(leader_rows)) {
        leader_radius <- sqrt(
          sites$x[leader_rows]^2 + sites$y[leader_rows]^2
        )
        spread <- diff(range(leader_radius, na.rm = TRUE))
        progress <- if (is.finite(spread) && spread > 1e-8) {
          (leader_radius - min(leader_radius, na.rm = TRUE)) / spread
        } else seq(0, 1, length.out = length(leader_rows))
        sites$x[leader_rows] <- sites$x[leader_rows] + dx * progress
      }
      sites$x[row] <- original_x + dx
      if ("label_attachment_x" %in% names(sites) &&
          is.finite(sites$label_attachment_x[row])) {
        sites$label_attachment_x[row] <- sites$label_attachment_x[row] + dx
      }
    }
    registry[[id]]$restriction_site <- sites
  }
  registry
}

# Restriction-site labels already own a polar contour/fan solver that keeps
# sparse labels at their natural angles and opens dense clusters locally. Do
# not flatten those labels onto Cartesian side rails when feature callouts are
# also present: that destroys the genomic contour and produces very long
# leaders. Each layer-local solver retains ownership of its candidate geometry;
# this final coord-owned pass only registers occupied exterior boxes and moves
# colliding feature callouts to the nearest free outer track.
ggchord_share_external_annotations <- function(registry, layout) {
  if (!isTRUE(layout$circular) || !length(registry)) return(registry)
  registry <- ggchord_expand_dense_restriction_fans(registry)
  units_per_inch <- layout$text_units_per_inch %||% .30

  restriction_boxes <- data.frame()
  for (id in names(registry)) {
    sites <- registry[[id]]$restriction_site
    if (!is.data.frame(sites) || !nrow(sites)) next
    rows <- sites$.component %in% "label" & !is.na(sites$label)
    if (!any(rows)) next
    if (!"outer_annotation_region" %in% names(sites))
      sites$outer_annotation_region <- NA_character_
    if (!"outer_track" %in% names(sites)) sites$outer_track <- NA_integer_
    if (!"outer_slot" %in% names(sites)) sites$outer_slot <- NA_integer_
    sites$outer_annotation_region[rows] <- "restriction"
    missing_track <- rows & is.na(sites$outer_track)
    missing_slot <- rows & is.na(sites$outer_slot)
    sites$outer_track[missing_track] <- 1L
    sites$outer_slot[missing_slot] <- 0L
    registry[[id]]$restriction_site <- sites
    labels <- sites[rows, , drop = FALSE]
    labels$text <- labels$label
    labels$text_x <- labels$x
    labels$text_y <- labels$y
    labels$text_angle <- labels$angle %||% 0
    restriction_boxes <- rbind(restriction_boxes, ggchord_text_boxes(
      labels, units_per_inch = units_per_inch, box_padding = .045
    ))
  }

  occupied <- restriction_boxes
  for (id in names(registry)) {
    geometry <- registry[[id]]$gene_label_repel
    if (!is.data.frame(geometry) || !nrow(geometry) ||
        !"feature_label_mode" %in% names(geometry)) next
    rows <- which(geometry$.component %in% "text" &
      geometry$feature_label_mode %in% "external" &
      !is.na(geometry$label) & nzchar(geometry$label))
    if (!length(rows)) next
    if (!"outer_annotation_region" %in% names(geometry)) {
      geometry$outer_annotation_region <- NA_character_
      geometry$outer_track <- NA_integer_
      geometry$outer_slot <- NA_integer_
    }
    rows <- rows[order(geometry$anchor_position[rows], rows)]
    for (row in rows) {
      original_x <- geometry$text_x[row]
      original_y <- geometry$text_y[row]
      radius <- sqrt(original_x^2 + original_y^2)
      if (!is.finite(radius) || radius < 1e-8) next
      direction <- c(original_x, original_y) / radius
      tangent <- c(-direction[2], direction[1])
      measured <- ggchord_text_boxes(
        geometry[row, , drop = FALSE], units_per_inch = units_per_inch,
        box_padding = .045
      )
      step <- max(.055, measured$h * .90)
      chosen <- measured
      chosen_track <- 1L
      candidate <- geometry[row, , drop = FALSE]
      candidates <- expand.grid(
        track = 1:8,
        tangent = c(0L, as.vector(rbind(1:10, -1:-10)))
      )
      # Stay in the natural polar sector while space is available.  A new
      # radial band is more disruptive than a nearby fan slot, but the bounded
      # tangent search prevents a callout from migrating around the circle.
      candidates$score <- (candidates$track - 1L) * 1.65 +
        abs(candidates$tangent) * .58
      candidates <- candidates[order(candidates$score,
        abs(candidates$tangent), candidates$tangent < 0), , drop = FALSE]
      chosen_tangent <- 0L
      for (candidate_index in seq_len(nrow(candidates))) {
        track <- candidates$track[candidate_index]
        tangent_index <- candidates$tangent[candidate_index]
        radial_shift <- (track - 1L) * step
        tangent_shift <- tangent_index * max(step, measured$h * 1.08)
        candidate$text_x <- original_x + direction[1] * radial_shift +
          tangent[1] * tangent_shift
        candidate$text_y <- original_y + direction[2] * radial_shift +
          tangent[2] * tangent_shift
        candidate_radius <- sqrt(candidate$text_x^2 + candidate$text_y^2)
        if (candidate_radius < radius) {
          correction <- radius - candidate_radius
          candidate$text_x <- candidate$text_x + direction[1] * correction
          candidate$text_y <- candidate$text_y + direction[2] * correction
        }
        candidate$x <- candidate$text_x
        candidate$y <- candidate$text_y
        box <- ggchord_text_boxes(
          candidate, units_per_inch = units_per_inch, box_padding = .045
        )
        collision <- nrow(occupied) &&
          any(ggchord_oriented_box_overlaps(box, occupied))
        chosen <- box
        chosen_track <- track
        chosen_tangent <- tangent_index
        if (!collision) break
      }
      dx <- candidate$text_x - original_x
      dy <- candidate$text_y - original_y
      geometry$text_x[row] <- candidate$text_x
      geometry$text_y[row] <- candidate$text_y
      geometry$x[row] <- candidate$text_x
      geometry$y[row] <- candidate$text_y
      geometry$outer_annotation_region[row] <- "feature"
      geometry$outer_track[row] <- chosen_track
      geometry$outer_slot[row] <- chosen_tangent
      segment_rows <- which(geometry$.component %in% "segment" &
        geometry$source_row == geometry$source_row[row])
      if (length(segment_rows)) {
        geometry$xend[segment_rows] <- geometry$xend[segment_rows] + dx
        geometry$yend[segment_rows] <- geometry$yend[segment_rows] + dy
        if ("x1" %in% names(geometry)) {
          geometry$x1[segment_rows] <- geometry$x1[segment_rows] + dx
          geometry$y1[segment_rows] <- geometry$y1[segment_rows] + dy
        }
        geometry$outer_annotation_region[segment_rows] <- "feature_leader"
        geometry$outer_track[segment_rows] <- chosen_track
        geometry$outer_slot[segment_rows] <- chosen_tangent
      }
      occupied <- rbind(occupied, chosen)
    }
    registry[[id]]$gene_label_repel <- geometry
  }
  registry
}
