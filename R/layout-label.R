# Automatic and fixed gene-label layout.
ggchord_layout_label_step <- quote({
  # ====================================================================
  # Step 8b: wrap gene labels; optionally arrange them automatically
  # ====================================================================
  gene_label_segments <- data.frame(x0 = numeric(0), y0 = numeric(0),
                                    x1 = numeric(0), y1 = numeric(0),
                                    group = integer(0),
                                    alpha = numeric(0),
                                    occluded = logical(0),
                                    linetype = character(0),
                                    stringsAsFactors = FALSE)
  gene_label_clip_units <- NA_real_
  if (nrow(gene_labels) > 0) {
    if (!is.null(gene_label_wrap)) {
      gene_labels$text <- ggchord_label_wrap_text(gene_labels$text,
                                                  gene_label_wrap)
    }
    units_per_inch <- text_units_per_inch
    if (isTRUE(gene_label_repel_layer)) {
      # Layout modes own these internal clearances. Keeping them out of the
      # public API prevents combinations that violate the geometry invariants.
      layout_box_padding <- if (gene_label_layout %in% c("auto", "callout")) {
        # Physical inches on each side of a text box. The previous 0.18-inch
        # value made an eight-label vertical rail more than twice as tall as
        # necessary and forced top/bottom labels onto extra rows. About 1 mm
        # keeps labels visually separate while letting the cardinal rails use
        # the available horizontal and vertical perimeter efficiently.
        0.05
      } else {
        0.04
      }
      layout_point_padding <- if (gene_label_layout %in% c("auto", "callout")) {
        0.08
      } else {
        0.05
      }
      layout_min_segment <- 0.02
      layout_units <- text_units_per_inch
      base_gene_labels <- gene_labels
      layout_result <- NULL

      if (!identical(gene_label_layout, "feature") &&
          is.null(gene_label_wrap) &&
          !identical(gene_label_fit, "none")) {
        initial_obstacles <- ggchord_text_obstacle_boxes(
          seq_labels_df, axis_ticks, show_axis,
          units_per_inch = layout_units
        )
        fit_labels <- base_gene_labels
        # auto/radial ultimately draw horizontal text. Measure that final
        # orientation here; measuring the temporary tangent angle could make a
        # long label look artificially narrow and skip adaptive wrapping.
        if (!identical(gene_label_layout, "arc")) {
          fit_labels$text_angle <- 0
        }
        fit_labels <- ggchord_fit_label_text(
          fit_labels,
          fit = gene_label_fit,
          max_lines = gene_label_max_lines,
          units_per_inch = layout_units,
          box_padding = layout_box_padding,
          repel_boxes = initial_obstacles
        )
        base_gene_labels$text <- fit_labels$text
        gene_labels <- base_gene_labels
      }

      if (gene_label_layout %in% c("radial", "auto")) {
        # Reserve physical space for text on both sides before converting to
        # data units. This estimate depends on the device and fitted font
        # metrics, never on the positions produced by the layout solver.
        physical <- ggchord_text_boxes(base_gene_labels, units_per_inch = 1)
        device <- if (grDevices::dev.cur() == 1L) c(8, 8) else
          grDevices::dev.size("in")
        reserve <- c(1.5 + 2 * max(physical$w),
                     0.75 + 2 * max(physical$h) + 0.8)
        usable <- pmax(device - reserve, device * 0.35)
        layout_units <- max(diff(range(compact_x)) / usable[1],
                            diff(range(compact_y)) / usable[2])
        text_units_per_inch <- layout_units
      }

      # A second deterministic pass lets fixed obstacles and label boxes use
      # the same device-derived physical scale without feeding the expanded
      # label limits back into the estimate (which would create excess blank
      # space around already-distant labels).
      layout_passes <- if (identical(gene_label_layout, "arc")) 2L else 1L
      for (layout_pass in seq_len(layout_passes)) {
        layout_obstacles <- ggchord_text_obstacle_boxes(
          seq_labels_df, axis_ticks, show_axis,
          units_per_inch = layout_units
        )
        if (identical(gene_label_layout, "feature")) {
          layout_result <- ggchord_feature_label_lanes(
            base_gene_labels, gene_polys, seq_arcs,
            units_per_inch = layout_units,
            max_overlaps = gene_label_repel_max_overlaps,
            allow_external = feature_label_external
          )
        } else if (gene_label_layout %in% c("auto", "callout")) {
          layout_result <- ggchord_auto_label_lanes(
            base_gene_labels, seq_arcs, side = gene_label_side,
            units_per_inch = layout_units, box_padding = layout_box_padding,
            point_padding = layout_point_padding, repel_boxes = layout_obstacles)
        } else if (identical(gene_label_layout, "radial")) {
          layout_result <- ggchord_radial_label_lanes(
            base_gene_labels, seq_arcs, side = gene_label_side,
            units_per_inch = layout_units, box_padding = layout_box_padding,
            point_padding = layout_point_padding, repel_boxes = layout_obstacles)
        } else {
          layout_result <- ggchord_offset_label_tracks(
            base_gene_labels, seq_arcs,
            side = gene_label_side,
            orientation = if (identical(gene_label_layout, "arc")) {
              "arc"
            } else {
              "horizontal"
            },
            units_per_inch = layout_units,
            box_padding = layout_box_padding,
            point_padding = layout_point_padding,
            repel_boxes = layout_obstacles
          )
        }
        gene_labels <- layout_result$labels

      }

      label_lanes <- layout_result$lanes
      label_directions <- layout_result$directions
      draw_segment <- layout_result$draw_segment
      gene_labels$label_layout <- gene_label_layout
      gene_labels$label_track <- layout_result$tracks
      final_obstacles <- ggchord_text_obstacle_boxes(
        seq_labels_df, axis_ticks, show_axis,
        units_per_inch = layout_units
      )
      if (!identical(gene_label_layout, "feature")) {
        gene_labels <- ggchord_hide_conflicted_labels(
          gene_labels,
          max_overlaps = gene_label_repel_max_overlaps,
          units_per_inch = layout_units,
          repel_boxes = final_obstacles
        )
      }
      draw_segment <- draw_segment & !is.na(gene_labels$text) &
        nzchar(gene_labels$text)

      if (gene_label_layout %in% c("radial", "auto")) {
        rows <- which(draw_segment)
        gene_label_segments <- data.frame(
          x0 = c(gene_labels$anchor_x[rows], gene_labels$.radial_bend_x[rows]),
          y0 = c(gene_labels$anchor_y[rows], gene_labels$.radial_bend_y[rows]),
          x1 = c(gene_labels$.radial_bend_x[rows], gene_labels$text_x[rows]),
          y1 = c(gene_labels$.radial_bend_y[rows], gene_labels$text_y[rows]),
          group = rep(rows, 2L))
      } else if (identical(gene_label_layout, "arc")) {
        rows <- which(draw_segment)
        gene_label_segments <- data.frame(
          x0 = gene_labels$anchor_x[rows],
          y0 = gene_labels$anchor_y[rows],
          x1 = gene_labels$text_x[rows],
          y1 = gene_labels$text_y[rows],
          group = rows,
          stringsAsFactors = FALSE
        )
      } else if (identical(gene_label_layout, "feature")) {
        rows <- which(draw_segment)
        gene_label_segments <- data.frame(
          x0 = gene_labels$anchor_x[rows],
          y0 = gene_labels$anchor_y[rows],
          x1 = gene_labels$text_x[rows],
          y1 = gene_labels$text_y[rows],
          group = rows,
          stringsAsFactors = FALSE
        )
      } else {
        gene_label_segments <- ggchord_repel_segments(
          gene_labels, min_segment_length = layout_min_segment
        )
        gene_label_segments <- gene_label_segments[
          draw_segment[gene_label_segments$group], , drop = FALSE
        ]
        if (nrow(gene_label_segments) > 0) {
          seg <- gene_label_segments
          text_boxes <- ggchord_text_boxes(
            gene_labels, units_per_inch = text_units_per_inch
          )
          bends <- ggchord_elbow_bends(
            gene_labels, text_boxes$w, text_boxes$h,
            directions = label_directions
          )
          bx <- bends$x[seg$group]
          by <- bends$y[seg$group]
          elbow <- data.frame(
            x0 = c(seg$x0, bx),
            y0 = c(seg$y0, by),
            x1 = c(bx, seg$x1),
            y1 = c(by, seg$y1),
            group = c(seg$group, seg$group),
            stringsAsFactors = FALSE
          )
          direction <- label_directions[seg$group]
          local <- direction %in% c("top", "bottom")
          local_segments <- seg[local, , drop = FALSE]
          elbow_segments <- elbow[!c(local, local), , drop = FALSE]
          elbow_segments <- ggchord_collapse_crossed_elbows(
            elbow_segments, lanes = label_lanes
          )
          gene_label_segments <- rbind(elbow_segments, local_segments)
        }
      }

      if (nrow(gene_label_segments) > 0) {
        gene_label_segments <- gene_label_segments[
          (gene_label_segments$x1 - gene_label_segments$x0)^2 +
          (gene_label_segments$y1 - gene_label_segments$y0)^2 > 1e-16, , drop = FALSE]
      }
      if (nrow(gene_label_segments) > 0) {
        gene_label_clip_units <- layout_units
        gene_label_segments <- ggchord_clip_segments_to_labels(
          gene_label_segments, gene_labels,
          units_per_inch = gene_label_clip_units,
          include_own = gene_label_layout %in% c("arc", "feature"),
          overlap = gene_label_segment_overlap,
          overlap_alpha = gene_label_segment_overlap_alpha
        )
      }
      # "auto" is solid unless the requested side differs from the gene's
      # strand-based side, in which case the established dashed convention is
      # retained for every visible piece of that label's leader.
      if (nrow(gene_label_segments) > 0) {
        if (identical(gene_label_segment_linetype, "auto")) {
          flipped <- gene_labels$side_flipped[
            match(gene_label_segments$group, seq_len(nrow(gene_labels)))
          ]
          gene_label_segments$linetype <- ifelse(flipped, "dashed", "solid")
        } else {
          gene_label_segments$linetype <- rep(
            gene_label_segment_linetype,
            length.out = nrow(gene_label_segments)
          )
        }
      }
    } else if (identical(gene_label_overlap, "nudge")) {
      gene_labels <- ggchord_label_deoverlap(
        gene_labels, units_per_inch = units_per_inch
      )
    } else if (identical(gene_label_overlap, "hide")) {
      fixed_obstacles <- ggchord_text_obstacle_boxes(
        seq_labels_df, axis_ticks, show_axis,
        units_per_inch = units_per_inch, box_padding = 0.01
      )
      gene_labels <- ggchord_label_prune_overlaps(
        gene_labels, units_per_inch = units_per_inch,
        repel_boxes = fixed_obstacles
      )
    }
    # Drop hidden labels and remap segment group IDs so every segment keeps a
    # valid reference after max_overlaps removes an interior label row.
    visible_labels <- !is.na(gene_labels$text) & nzchar(gene_labels$text)
    if (nrow(gene_label_segments) > 0) {
      group_map <- match(seq_len(nrow(gene_labels)), which(visible_labels))
      keep_segments <- !is.na(group_map[gene_label_segments$group])
      gene_label_segments <- gene_label_segments[keep_segments, , drop = FALSE]
      gene_label_segments$group <- group_map[gene_label_segments$group]
    }
    gene_labels <- gene_labels[visible_labels, , drop = FALSE]
  }

})
