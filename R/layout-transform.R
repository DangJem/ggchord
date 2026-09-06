# Common geometry rotation and text measurement scale.
ggchord_layout_transform_step <- quote({
  # ====================================================================
  # Step 8: rotate all elements uniformly
  # ====================================================================
  rotate_df <- function(df) {
    if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) return(df)

    if (all(c("x", "y") %in% names(df))) {
      x0 <- df$x; y0 <- df$y
      df$x <- x0 * cos(rot_rad) - y0 * sin(rot_rad)
      df$y <- x0 * sin(rot_rad) + y0 * cos(rot_rad)
    }
    if (all(c("x0", "y0", "x1", "y1") %in% names(df))) {
      X <- df$x0; Y <- df$y0
      df$x0 <- X * cos(rot_rad) - Y * sin(rot_rad)
      df$y0 <- X * sin(rot_rad) + Y * cos(rot_rad)
      X1 <- df$x1; Y1 <- df$y1
      df$x1 <- X1 * cos(rot_rad) - Y1 * sin(rot_rad)
      df$y1 <- X1 * sin(rot_rad) + Y1 * cos(rot_rad)
    }
    if (all(c("label_x", "label_y") %in% names(df))) {
      LX <- df$label_x; LY <- df$label_y
      df$label_x <- LX * cos(rot_rad) - LY * sin(rot_rad)
      df$label_y <- LX * sin(rot_rad) + LY * cos(rot_rad)
    }
    # Geometry-relative label angles ("parallel"/"perpendicular") follow the
    # rotated geometry; absolute angles ("horizontal"/numeric) stay fixed.
    if (all(c("label_angle", "label_angle_relative") %in% names(df))) {
      rel <- !is.na(df$label_angle_relative) & df$label_angle_relative
      df$label_angle[rel] <- df$label_angle[rel] + rotation
    }
    if (all(c("text_x", "text_y") %in% names(df))) {
      TX <- df$text_x; TY <- df$text_y
      df$text_x <- TX * cos(rot_rad) - TY * sin(rot_rad)
      df$text_y <- TX * sin(rot_rad) + TY * cos(rot_rad)
      df$text_angle <- df$text_angle + rotation
      # Keep the text readable after the global rotation: labels whose angle
      # ends up in (90, 270) would be upside down, so flip them by 180 degrees
      # (and flip their justification so the text box stays anchored).
      flip <- df$text_angle > 90 & df$text_angle < 270
      if (any(flip)) {
        df$text_angle[flip] <- df$text_angle[flip] + 180
        if ("hjust" %in% names(df)) df$hjust[flip] <- 1 - df$hjust[flip]
      }
      df$text_angle <- df$text_angle %% 360
    }
    if (all(c("anchor_x", "anchor_y") %in% names(df))) {
      AX <- df$anchor_x; AY <- df$anchor_y
      df$anchor_x <- AX * cos(rot_rad) - AY * sin(rot_rad)
      df$anchor_y <- AX * sin(rot_rad) + AY * cos(rot_rad)
    }
    df
  }

  # Apply to all geometric elements (avoid deep copies: modify the original object's columns in place)
  seq_arcs <- lapply(seq_arcs, rotate_df)
  if (nrow(axis_lines) > 0) axis_lines <- rotate_df(axis_lines)
  if (nrow(axis_ticks) > 0) {
    axis_ticks <- rotate_df(axis_ticks)
    # Flip geometry-relative labels that would read upside-down so that the
    # text always reads outward (text whose reading direction points left is
    # rotated by 180 degrees).
    rel <- !is.na(axis_ticks$label_angle_relative) &
      axis_ticks$label_angle_relative
    if (any(rel)) {
      ta <- (axis_ticks$label_angle + 360) %% 360
      flip <- rel & ta > 90 & ta < 270
      axis_ticks$label_angle[flip] <- ta[flip] + 180
    }
    # Align labels outward (away from the chord) so they do not overlap the
    # axis lines / ticks. Horizontal labels use simple quadrant-based
    # justification; rotated labels are justified in the text's local frame so
    # that the text still extends away from the chord center.
    eps <- 1e-3
    horiz <- (axis_ticks$label_angle %% 360) < 0.5
    axis_ticks$label_hjust <- ifelse(
      horiz,
      ifelse(axis_ticks$label_x > eps, 0,
             ifelse(axis_ticks$label_x < -eps, 1, 0.5)),
      0.5
    )
    axis_ticks$label_vjust <- ifelse(
      horiz,
      ifelse(axis_ticks$label_y > eps, 1,
             ifelse(axis_ticks$label_y < -eps, 0, 0.5)),
      0.5
    )
    if (any(!horiz)) {
      phi <- atan2(axis_ticks$label_y, axis_ticks$label_x) * 180 / pi
      alpha <- (phi - axis_ticks$label_angle) * pi / 180
      ca <- cos(alpha)
      sa <- sin(alpha)
      idx <- !horiz
      axis_ticks$label_hjust[idx] <- ifelse(ca[idx] > 0.05, 0,
                                            ifelse(ca[idx] < -0.05, 1, 0.5))
      axis_ticks$label_vjust[idx] <- ifelse(sa[idx] > 0.05, 0,
                                            ifelse(sa[idx] < -0.05, 1, 0.5))
    }
  }
  if (!is.null(ribbon_polys)) ribbon_polys <- rotate_df(ribbon_polys)
  if (nrow(region_polys) > 0) region_polys <- rotate_df(region_polys)
  if (nrow(gene_labels) > 0) gene_labels <- rotate_df(gene_labels)
  if (nrow(seq_labels_df) > 0) seq_labels_df <- rotate_df(seq_labels_df)
  # Horizontal sequence labels: keep every label horizontal (independent of
  # the global rotation) and let the text extend away from the chord center
  # unless the user supplied an explicit justification.
  if (identical(seq_label_orientation, "horizontal") &&
      nrow(seq_labels_df) > 0) {
    seq_labels_df$text_angle <- 0
    if (is.null(seq_label_hjust)) {
      seq_labels_df$hjust <- ifelse(seq_labels_df$text_x >= 0, 0, 1)
    }
  }
  if (nrow(gene_polys) > 0) {
    gene_polys <- rotate_df(gene_polys)
    gene_polys <- gene_polys[with(gene_polys, order(group, ord)), ]
  }

  # Use the undecorated chord geometry as the single physical scale for all
  # text measurement. Labels must not enlarge their own scale estimate: that
  # feedback produced excessive margins on large devices and underestimated
  # boxes on small devices when the old six-inch fallback was triggered.
  compact_x <- c(
    unlist(lapply(seq_arcs, `[[`, "x"), use.names = FALSE),
    gene_polys$x
  )
  compact_y <- c(
    unlist(lapply(seq_arcs, `[[`, "y"), use.names = FALSE),
    gene_polys$y
  )
  text_units_per_inch <- ggchord_device_units_per_inch(compact_x, compact_y)


})
