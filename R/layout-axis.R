# Sequence-axis lines, ticks and labels.
ggchord_layout_axis_step <- quote({
  # Step 5: generate axes (lines, ticks, labels)
  # ====================================================================
  axis_lines <- data.frame(x = numeric(0), y = numeric(0),
                            accver = character(0),
                            stringsAsFactors = FALSE)
  axis_ticks <- data.frame(x0 = numeric(0), y0 = numeric(0),
                           x1 = numeric(0), y1 = numeric(0),
                           label = character(0), label_x = numeric(0),
                           label_y = numeric(0), size = numeric(0),
                           label_angle = numeric(0),
                           label_angle_relative = logical(0),
                           accver = character(0),
                           stringsAsFactors = FALSE)

  if (show_axis) {
    # Axis lines
    axis_lines <- do.call(rbind, lapply(seqs, function(id) {
      ref <- seq_refs[[id]]
      r0 <- ref$r0 - axisGap[id]
      angles <- seq(starts[id], ends[id], length.out = nSeg)
      pts <- map_to_curve_many(angles, r0, ref)
      data.frame(x = pts[, 1], y = pts[, 2], accver = id, stringsAsFactors = FALSE)
    }))

    # Tick marks
    axis_ticks <- do.call(rbind, lapply(seqs, function(id) {
      ref <- seq_refs[[id]]
      r0 <- ref$r0 - axisGap[id]

      majors <- axis_breaks[[id]] %||% breakPointsFunc(lens[id], axisMaj[id])
      major_labels <- axis_labels[[id]]
      if (isTRUE(circular) && any(abs(majors) < sqrt(.Machine$double.eps))) {
        # Zero and sequence length are the same circular seam. Keep only the
        # zero-side tick even when a default ggplot2 scale supplied both.
        keep_major <- abs(majors - lens[id]) >= sqrt(.Machine$double.eps)
        majors <- majors[keep_major]
        if (!is.null(major_labels)) major_labels <- major_labels[keep_major]
      }
      minors <- axis_minor_breaks[[id]]
      if (is.null(minors)) {
        minors <- unlist(lapply(seq_len(length(majors) - 1), function(i) {
          seq(majors[i], majors[i + 1], length.out = axisMin[id] + 2)[-c(1, axisMin[id] + 2)]
        }))
      }
      major_labels <- major_labels %||% as.character(majors)
      pts <- data.frame(
        pos = c(majors, minors),
        is_major = c(rep(TRUE, length(majors)), rep(FALSE, length(minors))),
        display_label = c(as.character(major_labels), rep(NA_character_, length(minors)))
      )
      origin_tick <- isTRUE(circular) &
        pts$is_major & abs(pts$pos) < sqrt(.Machine$double.eps)
      # Circular maps use the seam itself as the origin cue. A longer radial
      # mark is clearer than printing a redundant zero over the 12-o'clock
      # feature stack.
      pts$display_label[origin_tick] <- NA_character_

      # Label orientation for this sequence. "horizontal" keeps the text
      # horizontal in the rendered plot; "parallel" aligns the text with the
      # axis (tangent direction); "perpendicular" aligns it with the radial
      # direction; numeric values are absolute angles in degrees (ggplot2
      # convention: counter-clockwise from horizontal).
      orient_val <- axisLabelOrientation[[id]]
      relative_angle <- is.character(orient_val) &&
        tolower(orient_val) %in% c("parallel", "perpendicular")

      frac <- if (orientation[id] == 1) pts$pos / lens[id] else 1 - pts$pos / lens[id]
      angle <- starts[id] + frac * (ends[id] - starts[id])

      # Tangent direction at each tick position (used for "parallel" and
      # "perpendicular" orientations).
      fi <- findInterval(angle, ref$angles)
      fi <- pmax(1, pmin(fi, length(ref$angles) - 1))
      idx <- ifelse(abs(ref$angles[fi] - angle) <=
                      abs(ref$angles[fi + 1] - angle), fi, fi + 1)
      idx_next <- pmin(idx + 1, nrow(ref$path))
      idx_prev <- pmax(idx - 1, 1)
      dx_t <- ref$path$x[idx_next] - ref$path$x[idx]
      dy_t <- ref$path$y[idx_next] - ref$path$y[idx]
      last <- idx == nrow(ref$path)
      if (any(last)) {
        dx_t[last] <- ref$path$x[idx[last]] - ref$path$x[idx_prev[last]]
        dy_t[last] <- ref$path$y[idx[last]] - ref$path$y[idx_prev[last]]
      }
      base_angle <- atan2(dy_t, dx_t) * 180 / pi

      if (is.character(orient_val) && tolower(orient_val) == "horizontal") {
        label_angle <- rep(0, nrow(pts))
      } else if (is.character(orient_val) &&
                 tolower(orient_val) == "parallel") {
        label_angle <- base_angle
      } else if (is.character(orient_val) &&
                 tolower(orient_val) == "perpendicular") {
        label_angle <- base_angle + 90
      } else {
        label_angle <- suppressWarnings(as.numeric(orient_val))
        if (is.na(label_angle)) label_angle <- 0
        label_angle <- rep(label_angle, nrow(pts))
      }

      dir <- if (axisGap[id] >= 0) -1 else 1
      len <- ifelse(pts$is_major, axisMajLen[id], axisMinLen[id])
      len[origin_tick] <- len[origin_tick] * 2.2
      base <- map_to_curve_many(angle, r0, ref)
      tip <- map_to_curve_many(angle, r0 + len * dir, ref)
      lbl <- map_to_curve_many(angle, r0 + (len + labelOffset[id]) * dir, ref)

      data.frame(
        x0 = base[, 1], y0 = base[, 2],
        x1 = tip[, 1], y1 = tip[, 2],
        label = pts$display_label,
        label_x = lbl[, 1], label_y = lbl[, 2],
        size = labelSize[[id]],
        label_angle = label_angle,
        label_angle_relative = relative_angle,
        is_major = pts$is_major,
        is_origin = origin_tick,
        accver = id,
        stringsAsFactors = FALSE
      )
    }))
  }

})
