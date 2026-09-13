# Shared perimeter layout for circular-map annotations.

# Feature callouts and restriction-site labels remain separate public layers,
# but they cannot own separate pieces of perimeter space. This final pass runs
# after every layer has produced its preferred positions, then packs all
# horizontal external labels on common left/right rails. It preserves the
# vertical order of the biological anchors and only rewrites label positions
# and the final leader endpoint.
ggchord_share_external_annotations <- function(registry, layout) {
  if (!isTRUE(layout$circular) || !length(registry)) return(registry)

  entries <- list()
  add_entries <- function(layer_id, type, data, rows, anchor_x, anchor_y) {
    if (!length(rows)) return(invisible(NULL))
    for (j in seq_along(rows)) {
      row <- rows[j]
      face <- data$fontface[row] %||% 1
      if (is.numeric(face)) {
        face <- c("plain", "bold", "italic", "bold.italic")[
          pmax(1L, pmin(4L, as.integer(face[1L])))
        ]
      } else {
        face <- as.character(face[1L])
      }
      if (is.na(face) || !face %in% c(
          "plain", "bold", "italic", "bold.italic")) face <- "plain"
      entries[[length(entries) + 1L]] <<- data.frame(
        layer_id = layer_id, type = type, row = row,
        old_x = data$x[row], old_y = data$y[row],
        anchor_x = anchor_x[j], anchor_y = anchor_y[j],
        text = as.character(data$label[row]),
        size = as.numeric(data$size[row] %||% 2.9),
        family = as.character(data$family[row] %||% ""),
        fontface = face,
        lineheight = as.numeric(data$lineheight[row] %||% 1.2),
        stringsAsFactors = FALSE
      )
    }
    invisible(NULL)
  }

  for (layer_id in names(registry)) {
    feature <- registry[[layer_id]]$gene_label_repel
    if (is.data.frame(feature) && nrow(feature) &&
        all(c(".component", "feature_label_mode") %in% names(feature))) {
      rows <- which(feature$.component == "text" &
        feature$feature_label_mode == "external" &
        !is.na(feature$label) & nzchar(feature$label))
      if (length(rows)) {
        ax <- feature$anchor_x[rows] %||% feature$x[rows]
        ay <- feature$anchor_y[rows] %||% feature$y[rows]
        add_entries(layer_id, "feature", feature, rows, ax, ay)
      }
    }
    restriction <- registry[[layer_id]]$restriction_site
    if (is.data.frame(restriction) && nrow(restriction) &&
        ".component" %in% names(restriction)) {
      rows <- which(restriction$.component == "label" &
        !is.na(restriction$label) & nzchar(restriction$label))
      if (length(rows)) {
        # Use the biological leader root, not the restriction layer's former
        # label attachment. A local fan is allowed to swing a label across a
        # cardinal axis; using that displaced point here would assign the
        # annotation to the wrong global side and send its leader through the
        # plasmid.
        ax <- ay <- rep(NA_real_, length(rows))
        for (j in seq_along(rows)) {
          junction <- restriction$junction_id[rows[j]]
          tick_rows <- which(restriction$.component == "path" &
            restriction$restriction_component == "tick" &
            restriction$junction_id == junction)
          if (length(tick_rows)) {
            radius <- restriction$x[tick_rows]^2 + restriction$y[tick_rows]^2
            root <- tick_rows[which.max(radius)]
            ax[j] <- restriction$x[root]
            ay[j] <- restriction$y[root]
          }
        }
        invalid <- !is.finite(ax) | !is.finite(ay)
        ax[invalid] <- restriction$label_attachment_x[rows[invalid]] %||%
          restriction$x[rows[invalid]]
        ay[invalid] <- restriction$label_attachment_y[rows[invalid]] %||%
          restriction$y[rows[invalid]]
        add_entries(layer_id, "restriction", restriction, rows, ax, ay)
      }
    }
  }
  if (!length(entries)) return(registry)
  candidates <- do.call(rbind, entries)
  # A shared pass is material only when the two annotation systems coexist.
  if (!all(c("feature", "restriction") %in% candidates$type)) return(registry)

  measure <- data.frame(
    text_x = candidates$old_x, text_y = candidates$old_y,
    text = candidates$text, text_angle = 0, hjust = .5, vjust = .5,
    size = candidates$size, family = candidates$family,
    fontface = candidates$fontface, lineheight = candidates$lineheight
  )
  units_per_inch <- layout$text_units_per_inch %||% .35
  boxes <- ggchord_text_boxes(
    measure, units_per_inch = units_per_inch, box_padding = 0
  )
  callout_pad <- ifelse(candidates$type == "feature", .030, .014)
  half_width <- boxes$w / 2 + callout_pad
  half_height <- boxes$h / 2 + callout_pad * .62
  side <- ifelse(candidates$anchor_x < 0, "left", "right")
  near_axis <- abs(candidates$anchor_x) < .02
  side[near_axis] <- ifelse(candidates$old_x[near_axis] < 0,
    "left", "right")
  new_x <- candidates$old_x
  new_y <- candidates$old_y
  shared_order <- integer(nrow(candidates))

  for (side_name in c("left", "right")) {
    rows <- which(side == side_name)
    if (!length(rows)) next
    biological_y <- candidates$anchor_y[rows]
    # Ties are deterministic and retain the existing visual order.
    order_value <- rank(biological_y, ties.method = "first")
    new_y[rows] <- ggchord_pack_label_axis(
      biological_y, order_value,
      half_height[rows], half_height[rows], gap = .014
    )
    if (side_name == "right") {
      inner_edge <- max(candidates$old_x[rows] - half_width[rows],
        candidates$anchor_x[rows] + .055, na.rm = TRUE)
      new_x[rows] <- inner_edge + half_width[rows]
    } else {
      inner_edge <- min(candidates$old_x[rows] + half_width[rows],
        candidates$anchor_x[rows] - .055, na.rm = TRUE)
      new_x[rows] <- inner_edge - half_width[rows]
    }
    shared_order[rows[order(biological_y, seq_along(rows))]] <-
      seq_along(rows)
  }

  for (i in seq_len(nrow(candidates))) {
    layer_id <- candidates$layer_id[i]
    type <- candidates$type[i]
    component <- if (type == "feature") "gene_label_repel" else
      "restriction_site"
    data <- registry[[layer_id]][[component]]
    if (!"external_annotation_type" %in% names(data)) {
      data$external_annotation_type <- rep(NA_character_, nrow(data))
      data$shared_external_side <- rep(NA_character_, nrow(data))
      data$shared_external_order <- rep(NA_integer_, nrow(data))
    }
    row <- candidates$row[i]
    old_center <- c(data$x[row], data$y[row])
    data$x[row] <- new_x[i]
    data$y[row] <- new_y[i]
    if ("text_x" %in% names(data)) data$text_x[row] <- new_x[i]
    if ("text_y" %in% names(data)) data$text_y[row] <- new_y[i]
    data$hjust[row] <- .5
    data$vjust[row] <- .5
    data$external_annotation_type[row] <- type
    data$shared_external_side[row] <- side[i]
    data$shared_external_order[row] <- shared_order[i]
    endpoint <- c(
      new_x[i] + if (side[i] == "right") -half_width[i] else half_width[i],
      new_y[i]
    )

    if (type == "feature") {
      group_value <- data$group[row]
      segments <- which(data$.component == "segment" &
        data$group == group_value)
      if (length(segments)) {
        distance <- (data$xend[segments] - old_center[1L])^2 +
          (data$yend[segments] - old_center[2L])^2
        terminal <- segments[which.min(distance)]
        data$xend[terminal] <- endpoint[1L]
        data$yend[terminal] <- pmax(new_y[i] - half_height[i],
          pmin(data$yend[terminal], new_y[i] + half_height[i]))
      }
    } else {
      junction <- data$junction_id[row]
      leaders <- which(data$.component == "path" &
        data$restriction_component == "leader" &
        data$junction_id == junction)
      if (length(leaders)) {
        endpoint[2L] <- pmax(new_y[i] - half_height[i],
          pmin(endpoint[2L], new_y[i] + half_height[i]))
        # Once a label enters the shared perimeter, its old layer-local fan is
        # no longer a valid route. Rebuild the two-piece leader from its real
        # radial root through a short independent stub to the newly assigned
        # text edge. Keeping both endpoint orders aligned with anchor_y is what
        # prevents the long fan segments from crossing one another.
        path_groups <- unique(data$group[leaders])
        path_groups <- path_groups[order(path_groups)]
        root_rows <- leaders[data$group[leaders] == path_groups[1L]]
        root <- c(candidates$anchor_x[i], candidates$anchor_y[i])
        radial_length <- sqrt(sum(root^2))
        radial <- if (is.finite(radial_length) && radial_length > 0) {
          root / radial_length
        } else c(if (side[i] == "right") 1 else -1, 0)
        stub <- root + radial * .022
        write_segment <- function(rows, from, to) {
          fraction <- seq(0, 1, length.out = length(rows))
          data$x[rows] <<- from[1L] + fraction * (to[1L] - from[1L])
          data$y[rows] <<- from[2L] + fraction * (to[2L] - from[2L])
        }
        write_segment(root_rows, root, stub)
        if (length(path_groups) > 1L) {
          fan_rows <- leaders[data$group[leaders] == path_groups[2L]]
          write_segment(fan_rows, stub, endpoint)
        } else {
          write_segment(root_rows, root, endpoint)
        }
      }
      order_value <- if (side[i] == "right") {
        "enzyme_position"
      } else "position_enzyme"
      data$label_order[row] <- order_value
      data$label_connection_side[row] <- if (side[i] == "right")
        "left" else "right"
      if (all(c("enzyme_label", "coordinate_label") %in% names(data))) {
        data$label[row] <- if (order_value == "enzyme_position") {
          paste(data$enzyme_label[row], data$coordinate_label[row])
        } else paste(data$coordinate_label[row], data$enzyme_label[row])
      }
    }
    registry[[layer_id]][[component]] <- data
  }
  registry
}
