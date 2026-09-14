# Perimeter annotation coordination for circular maps.

# Restriction-site labels already own a polar contour/fan solver that keeps
# sparse labels at their natural angles and opens dense clusters locally. Do
# not flatten those labels onto Cartesian side rails when feature callouts are
# also present: that destroys the genomic contour and produces very long
# leaders. Each layer-local solver retains ownership of its candidate geometry;
# this final coord-owned pass only registers occupied exterior boxes and moves
# colliding feature callouts to the nearest free outer track.
ggchord_share_external_annotations <- function(registry, layout) {
  if (!isTRUE(layout$circular) || !length(registry)) return(registry)
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
        track = 1:5,
        tangent = c(0L, as.vector(rbind(1:8, -1:-8)))
      )
      candidates$score <- (candidates$track - 1L) * 1.15 +
        abs(candidates$tangent) * .72
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
