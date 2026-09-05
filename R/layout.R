# layout.R - core chord layout computation
# Pre-computes Cartesian (x, y) coordinates for all geometric elements from standardized parameters,
# for direct use by the geom_* layers

#' Compute the chord layout
#'
#' Pre-computes the coordinates of all geometric elements (sequence arcs, ribbons, gene arrows, axes, etc.)
#' into Cartesian (x, y) coordinates and stores them in a layout list.
#'
#' @param seqs Vector of sequence IDs (order already processed)
#' @param lens Named vector of sequence lengths (names = seq_id)
#' @param seq_labels Named vector of sequence labels
#' @param seqRadius Named vector of sequence radii
#' @param seq_curvature Named vector of sequence curvatures
#' @param orientation Named vector of sequence orientations (1 or -1)
#' @param seq_gap Named vector of sequence gap proportions
#' @param ribbonGap Named vector of ribbon gaps
#' @param ribbon_gap_auto Whether ribbon endpoints may move closer to sequence
#'   arcs when no gene or feature geometry occupies the local interval.
#' @param ribbon_obstacles Optional normalized gene/feature obstacle table.
#' @param ribbon_data Alignment data (already validated)
#' @param gene_data Gene data (already validated)
#' @param draw_gene_geometry Whether gene/feature polygons should be generated.
#' @param gene_label_layout Character, default "radial". Deterministic
#'   automatic layout: "radial", "auto", or "arc".
#' @param gene_label_side Character, default "auto". Which side of the arc the
#'   labels sit on: "auto" (strand-based placement), "inside" (toward the chord
#'   center) or "outside" (away from the center, avoiding ribbon overlap).
#' @param gene_label_segment_linetype Character or numeric, default "auto".
#'   Leader-line linetype; "auto" uses solid lines except for labels moved to
#'   the other side of their arc, which use dashed lines.
#' @param gene_label_segment_overlap Character: "fade", "clip", or "show".
#' @param gene_label_segment_overlap_alpha Relative opacity for covered leader
#'   portions when `gene_label_segment_overlap = "fade"`.
#' @param gene_label_orientation Fixed-label text orientation: "radial",
#'   "tangent", or "horizontal". Automatic repel layouts manage their own
#'   text orientation.
#' @param gene_label_overlap Fixed-label collision policy: "hide", "nudge",
#'   or "allow".
#' @param seq_label_orientation Character, default "arc". Sequence label text
#'   orientation: "arc" (rotated along the arc, kept readable) or "horizontal"
#'   (all labels horizontal, extending away from the chord center).
#' @param seq_label_hjust Optional named vector or NULL, default NULL. Per-seq
#'   horizontal justification; NULL uses 0.5 (arc mode) or a side-based value
#'   (horizontal mode).
#' @param seq_label_vjust Optional named vector or NULL, default NULL. Per-seq
#'   vertical justification; NULL uses 0.5.
#' @param rotation Global rotation angle (degrees)
#' @param debug Whether to output debug information
#'
#' @return A chord layout list
#' @keywords internal
compute_chord_layout <- function(
    seqs, lens, seq_labels, seq_colors,
    seqRadius, seq_curvature, orientation, seq_gap,
    # Ribbon parameters
    ribbon_data = NULL, ribbonGap,
    ribbon_gap_auto = FALSE, ribbon_obstacles = NULL,
    ribbon_color_scheme, ribbon_colors, ribbon_alpha,
    ribbon_color_by = NULL,
    ribbon_color_limits = NULL,
    ribbon_color_breaks = NULL,
    ribbon_color_name = NULL,
    ribbon_alpha_by = NULL,
    ribbon_alpha_range = c(0.15, 0.9),
    ribbon_outline_by = NULL,
    ribbon_outline_colors = NULL,
    ribbon_linetype_by = NULL,
    ribbon_linetypes = NULL,
    ribbon_direction = "none",
    ribbon_direction_colors = c(same = "black", reverse = "grey50"),
    ribbon_direction_linetypes = c(same = "solid", reverse = "dashed"),
    ribbon_direction_alpha = c(same = 1, reverse = 0.45),
    ribbon_ctrl_point,
    region_data = NULL,
    region_fill = "#F59E0B",
    region_color = "#B45309",
    region_alpha = 0.25,
    region_width = 0.08,
    region_offset = 0,
    region_side = "inside",
    ribbon_highlight_rows = integer(0),
    # Gene parameters
    gene_data = NULL, draw_gene_geometry = TRUE,
    geneGap, geneWidth,
    geneLabelRadialOffset, geneLabelCircumOffset,
    geneLabelCircumLimit, geneLabelRotation,
    gene_label_show, gene_label_size,
    gene_label_wrap = NULL,
    gene_label_fit = "wrap",
    gene_label_max_lines = 2L,
    gene_label_orientation = "horizontal",
    gene_label_overlap = "hide",
    gene_label_repel_layer = FALSE,
    gene_label_repel_max_overlaps = Inf,
    gene_label_layout = "radial",
    gene_label_side = "auto",
    gene_label_segment_overlap = "fade",
    gene_label_segment_overlap_alpha = 0.18,
    gene_label_segment_linetype = "auto",
    gene_color_scheme, gene_colors, gene_order,
    feature_shape_pal = NULL, feature_shape_order = NULL,
    # Sequence label parameters
    seq_label_text = NULL, seq_label_radius = NULL,
    seq_label_rotation = NULL, seq_label_size = NULL,
    seq_label_orientation = "arc",
    seq_label_hjust = NULL, seq_label_vjust = NULL,
    # Axis parameters
    axisGap, axisMaj, axisMajLen, axisMin, axisMinLen,
    labelSize, labelOffset, axisLabelOrientation,
    axis_breaks = NULL, axis_minor_breaks = NULL, axis_labels = NULL,
    axis_label_hide_overlaps = FALSE,
    show_axis,
    # Global parameters
    rotation, debug = FALSE,
    geometry_cache = NULL
) {
  n <- length(seqs)

  # ====================================================================
  # Step 1: compute angle allocation
  # ====================================================================
  total_circ <- 2 * pi

  total_gap_prop <- sum(seq_gap)

  if (total_gap_prop >= 1) {
    ggchord_stop("The sum of seq_gap cannot exceed 1 (no space left for sequences)")
  }

  seq_total_prop <- 1 - total_gap_prop
  sum_lens <- sum(lens)
  theta <- (lens / sum_lens) * total_circ * seq_total_prop
  gap_rads <- total_circ * seq_gap

  # Compute the start and end angles of each sequence
  starts <- numeric(n)
  names(starts) <- seqs
  starts[1] <- 0

  if (n > 1) {
    for (i in 2:n) {
      starts[i] <- starts[i - 1] + theta[i - 1] + gap_rads[i - 1]
    }
  }
  ends <- starts + theta
  names(ends) <- seqs

  # ====================================================================
  # Step 2: convert to radians
  # ====================================================================
  rot_rad <- rotation * pi / 180

  # ====================================================================
  # Step 3: generate reference paths and map_to_curve
  # ====================================================================
  nSeg <- 500
  nRef <- 2000

  # High-resolution reference path for each sequence
  reference_key <- list(
    seqs = seqs, starts = starts, ends = ends,
    radius = seqRadius, curvature = seq_curvature
  )
  if (!is.null(geometry_cache) &&
      identical(geometry_cache$reference_key, reference_key)) {
    seq_refs <- geometry_cache$seq_refs
  } else {
    seq_refs <- lapply(seqs, function(id) {
      path <- generate_curvature_path(
        starts[id], ends[id], seqRadius[id], seq_curvature[id], n_points = nRef
      )
      angles <- seq(starts[id], ends[id], length.out = nRef)
      list(path = path, angles = angles, r0 = seqRadius[id])
    })
    names(seq_refs) <- seqs
    if (!is.null(geometry_cache)) {
      geometry_cache$reference_key <- reference_key
      geometry_cache$seq_refs <- seq_refs
    }
  }

  # Curve coordinate mapping function
  map_to_curve_many <- function(angle, radius, ref) {
    n <- length(angle)
    nref <- length(ref$angles)

    # Vectorised nearest-reference lookup (O(n log n) overall).
    fi <- findInterval(angle, ref$angles)
    fi <- pmax(1, pmin(fi, nref - 1))
    idx <- ifelse(abs(ref$angles[fi] - angle) <=
                    abs(ref$angles[fi + 1] - angle), fi, fi + 1)

    base <- ref$path[idx, , drop = FALSE]
    idx_next <- pmin(idx + 1, nref)
    idx_prev <- pmax(idx - 1, 1)
    dx <- ref$path$x[idx_next] - base$x
    dy <- ref$path$y[idx_next] - base$y
    last <- idx == nref
    if (any(last)) {
      dx[last] <- base$x[last] - ref$path$x[idx_prev[last]]
      dy[last] <- base$y[last] - ref$path$y[idx_prev[last]]
    }

    norm_x <- -dy
    norm_y <- dx
    nl <- sqrt(norm_x^2 + norm_y^2)
    ok <- nl > 0
    norm_x[ok] <- norm_x[ok] / nl[ok]
    norm_y[ok] <- norm_y[ok] / nl[ok]

    offset <- radius - ref$r0
    cbind(x = base$x + norm_x * offset,
          y = base$y + norm_y * offset)
  }

  map_to_curve <- function(angle, radius, ref) {
    # ref$angles is sorted, so use findInterval (O(log n)) instead of a
    # full linear scan to locate the nearest reference angle.
    fi <- findInterval(angle, ref$angles)
    if (fi < 1) fi <- 1
    if (fi >= length(ref$angles)) fi <- length(ref$angles) - 1
    idx <- if (abs(ref$angles[fi] - angle) <= abs(ref$angles[fi + 1] - angle)) fi else fi + 1
    base <- ref$path[idx, ]
    if (idx < nrow(ref$path)) {
      dx <- ref$path$x[idx + 1] - base$x
      dy <- ref$path$y[idx + 1] - base$y
    } else {
      dx <- base$x - ref$path$x[idx - 1]
      dy <- base$y - ref$path$y[idx - 1]
    }
    norm <- c(-dy, dx)
    nl <- sqrt(sum(norm^2))
    if (nl > 0) norm <- norm / nl
    offset <- radius - ref$r0
    c(x = base$x + norm[1] * offset,
      y = base$y + norm[2] * offset)
  }

  gene_track_radius <- function(gene, sid, strand) {
    side <- if (".feature_stack_side" %in% names(gene)) {
      as.character(gene[[".feature_stack_side"]])
    } else if (strand == "+") {
      "inside"
    } else {
      "outside"
    }
    lane <- if (".feature_stack_lane" %in% names(gene)) {
      as.numeric(gene[[".feature_stack_lane"]])
    } else 0
    spacing <- if (".feature_stack_spacing" %in% names(gene)) {
      as.numeric(gene[[".feature_stack_spacing"]])
    } else 0
    direction <- if (identical(side, "inside")) -1 else 1
    seqRadius[sid] + direction *
      (geneGap[[sid]][strand] + lane * spacing)
  }

  # ====================================================================
  # Step 4: generate sequence arcs (outer layer)
  # ====================================================================
  seq_arcs <- stats::setNames(lapply(seqs, function(id) {
    path_data <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id], seq_curvature[id], nSeg
    )
    path_data$seq_id <- id
    if (orientation[id] == -1) {
      path_data <- path_data[nrow(path_data):1, ]
    }
    attr(path_data, "ggchord_path_direction") <- unname(orientation[id])
    path_data
  }), seqs)

  axis_reference_x <- unlist(lapply(seq_arcs, `[[`, "x"), use.names = FALSE)
  axis_reference_y <- unlist(lapply(seq_arcs, `[[`, "y"), use.names = FALSE)
  axis_units_per_inch <- ggchord_device_units_per_inch(
    axis_reference_x, axis_reference_y
  )
  axisGap <- axisGap * axis_units_per_inch
  axisMajLen <- axisMajLen * axis_units_per_inch
  axisMinLen <- axisMinLen * axis_units_per_inch
  labelOffset <- labelOffset * axis_units_per_inch

  # ====================================================================
  # Step 5: generate axes (lines, ticks, labels)
  # ====================================================================
  axis_lines <- data.frame(x = numeric(0), y = numeric(0),
                            seq_id = character(0),
                            stringsAsFactors = FALSE)
  axis_ticks <- data.frame(x0 = numeric(0), y0 = numeric(0),
                           x1 = numeric(0), y1 = numeric(0),
                           label = character(0), label_x = numeric(0),
                           label_y = numeric(0), size = numeric(0),
                           label_angle = numeric(0),
                           label_angle_relative = logical(0),
                           seq_id = character(0),
                           stringsAsFactors = FALSE)

  if (show_axis) {
    # Axis lines
    axis_lines <- do.call(rbind, lapply(seqs, function(id) {
      ref <- seq_refs[[id]]
      r0 <- ref$r0 - axisGap[id]
      angles <- seq(starts[id], ends[id], length.out = nSeg)
      pts <- map_to_curve_many(angles, r0, ref)
      data.frame(x = pts[, 1], y = pts[, 2], seq_id = id, stringsAsFactors = FALSE)
    }))

    # Tick marks
    axis_ticks <- do.call(rbind, lapply(seqs, function(id) {
      ref <- seq_refs[[id]]
      r0 <- ref$r0 - axisGap[id]

      majors <- axis_breaks[[id]] %||% breakPointsFunc(lens[id], axisMaj[id])
      minors <- axis_minor_breaks[[id]]
      if (is.null(minors)) {
        minors <- unlist(lapply(seq_len(length(majors) - 1), function(i) {
          seq(majors[i], majors[i + 1], length.out = axisMin[id] + 2)[-c(1, axisMin[id] + 2)]
        }))
      }
      major_labels <- axis_labels[[id]] %||% as.character(majors)
      pts <- data.frame(
        pos = c(majors, minors),
        is_major = c(rep(TRUE, length(majors)), rep(FALSE, length(minors))),
        display_label = c(as.character(major_labels), rep(NA_character_, length(minors)))
      )

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
        seq_id = id,
        stringsAsFactors = FALSE
      )
    }))
  }

  # ====================================================================
  # Step 6: generate ribbon polygons
  # ====================================================================
  ribbon_polys <- NULL
  ribbon_color_info <- list(scheme = ribbon_color_scheme)
  ribbon_use_outline <- !is.null(ribbon_outline_by) || identical(ribbon_direction, "outline")
  ribbon_use_linetype <- !is.null(ribbon_linetype_by) || identical(ribbon_direction, "linetype")

  resolve_discrete_ribbon_values <- function(vals, supplied, default_fun, arg_name) {
    uniq <- unique(vals)
    if (is.null(supplied)) {
      return(setNames(default_fun(length(uniq)), uniq))
    }
    if (is.null(names(supplied))) {
      if (length(supplied) == 1) {
        return(setNames(rep(supplied, length(uniq)), uniq))
      }
      if (length(supplied) == length(uniq)) {
        return(setNames(as.character(supplied), uniq))
      }
      ggchord_stop(arg_name, " must be length 1 or match the number of unique values")
    }
    unknown <- setdiff(names(supplied), uniq)
    if (length(unknown) > 0) {
      ggchord_stop(arg_name, " contains unknown value(s): ", paste(unknown, collapse = ", "))
    }
    out <- setNames(default_fun(length(uniq)), uniq)
    out[names(supplied)] <- as.character(supplied)
    out
  }

  if (!is.null(ribbon_data) && nrow(ribbon_data) > 0) {
    ribbons <- list()
    cntValid <- 0
    cntInvalid <- 0

    # Pre-process ribbon colors
    singleCol <- NULL
    queryCols <- NULL

    if (ribbon_color_scheme == "single") {
      singleCol <- if (length(ribbon_colors) > 1) ribbon_colors[[1]] else ribbon_colors
    } else if (ribbon_color_scheme %in% c("query", "subject")) {
      queryCols <- ribbon_colors
    }

    # Pre-extract columns for faster iteration and validation.
    n_ribbons <- nrow(ribbon_data)
    rib_q <- ribbon_data$qaccver
    rib_s <- ribbon_data$saccver
    rib_qstart <- ribbon_data$qstart
    rib_qend <- ribbon_data$qend
    rib_sstart <- ribbon_data$sstart
    rib_send <- ribbon_data$send
    rib_pident <- ribbon_data$pident

    valid <- rib_q != rib_s & rib_q %in% seqs & rib_s %in% seqs
    valid_idx <- which(valid)
    cntInvalid <- n_ribbons - length(valid_idx)

    ribbon_polys_list <- vector("list", length(valid_idx))
    ribbon_group <- integer(length(valid_idx))
    ribbon_fill <- character(length(valid_idx))
    ribbon_pident <- numeric(length(valid_idx))
    ribbon_value <- numeric(length(valid_idx))
    ribbon_alpha_vec <- numeric(length(valid_idx))
    ribbon_outline_vec <- character(length(valid_idx))
    ribbon_linetype_vec <- character(length(valid_idx))
    ribbon_dir_vec <- character(length(valid_idx))
    ribbon_q_gap <- numeric(length(valid_idx))
    ribbon_s_gap <- numeric(length(valid_idx))

    # With the default NULL ribbon_gap, determine spacing independently at
    # each ribbon endpoint. Only polygons that occupy the ribbon-facing side
    # of the sequence and overlap that endpoint's genomic interval count as
    # obstacles. Text and leader segments never enter this table.
    endpoint_ribbon_gap <- function(seq_id, start, end) {
      configured <- unname(ribbonGap[[seq_id]])
      if (!isTRUE(ribbon_gap_auto)) return(configured)

      close_gap <- min(configured, 0.035)
      if (is.null(ribbon_obstacles) || nrow(ribbon_obstacles) == 0L) {
        return(close_gap)
      }
      lo <- min(start, end)
      hi <- max(start, end)
      hit <- ribbon_obstacles$seq_id == seq_id &
        ribbon_obstacles$start <= hi & ribbon_obstacles$end >= lo &
        ribbon_obstacles$outer_offset > 0
      if (!any(hit)) return(close_gap)

      # The standard gene offset (0.10) plus half-width (0.025) and this
      # clearance reproduces the historical safe 0.15 gap where an obstacle
      # actually exists, while wider/custom-offset features remain protected.
      max(close_gap, max(ribbon_obstacles$outer_offset[hit]) + 0.025)
    }

    # Continuous / discrete value preparation on the valid rows only.
    alpha_norm <- rep(1, length(valid_idx))
    if (!is.null(ribbon_alpha_by)) {
      av <- as.numeric(ribbon_data[[ribbon_alpha_by]][valid_idx])
      if (length(av) > 1 && diff(range(av, na.rm = TRUE)) > 0) {
        rng <- range(av, na.rm = TRUE)
        alpha_norm <- (av - rng[1]) / (rng[2] - rng[1])
      } else {
        alpha_norm <- rep(0.5, length(av))
      }
      alpha_norm <- ribbon_alpha_range[1] +
        alpha_norm * (ribbon_alpha_range[2] - ribbon_alpha_range[1])
    }

    if (!is.null(ribbon_outline_by)) {
      ov <- as.character(ribbon_data[[ribbon_outline_by]][valid_idx])
      outline_map <- resolve_discrete_ribbon_values(
        ov, ribbon_outline_colors, chord_default_palette, "ribbon_outline_colors")
      ribbon_outline_vec <- unname(outline_map[ov])
    } else if (identical(ribbon_direction, "outline")) {
      ribbon_outline_vec <- rep("black", length(valid_idx))
    } else {
      ribbon_outline_vec <- rep("black", length(valid_idx))
    }

    if (!is.null(ribbon_linetype_by)) {
      lv <- as.character(ribbon_data[[ribbon_linetype_by]][valid_idx])
      linetype_map <- resolve_discrete_ribbon_values(
        lv, ribbon_linetypes,
        function(n) rep_len(c("solid", "dashed", "dotted", "dotdash",
                              "longdash", "twodash"), n),
        "ribbon_linetypes")
      ribbon_linetype_vec <- unname(linetype_map[lv])
    } else if (identical(ribbon_direction, "linetype")) {
      ribbon_linetype_vec <- rep("solid", length(valid_idx))
    } else {
      ribbon_linetype_vec <- rep("solid", length(valid_idx))
    }

    for (j in seq_along(valid_idx)) {
      i <- valid_idx[j]
      q <- rib_q[i]
      s <- rib_s[i]

      q_ref <- seq_refs[[q]]
      s_ref <- seq_refs[[s]]
      q_gap <- endpoint_ribbon_gap(q, rib_qstart[i], rib_qend[i])
      s_gap <- endpoint_ribbon_gap(s, rib_sstart[i], rib_send[i])
      rq <- seqRadius[q] + q_gap
      rs <- seqRadius[s] + s_gap
      ribbon_q_gap[j] <- q_gap
      ribbon_s_gap[j] <- s_gap

      q_frac_start <- if (orientation[q] == 1) (rib_qstart[i] - 1) / lens[q] else 1 - (rib_qstart[i] - 1) / lens[q]
      q_angle_start <- starts[q] + q_frac_start * (ends[q] - starts[q])
      q_frac_end <- if (orientation[q] == 1) (rib_qend[i] - 1) / lens[q] else 1 - (rib_qend[i] - 1) / lens[q]
      q_angle_end <- starts[q] + q_frac_end * (ends[q] - starts[q])

      s_frac_start <- if (orientation[s] == 1) (rib_sstart[i] - 1) / lens[s] else 1 - (rib_sstart[i] - 1) / lens[s]
      s_angle_start <- starts[s] + s_frac_start * (ends[s] - starts[s])
      s_frac_end <- if (orientation[s] == 1) (rib_send[i] - 1) / lens[s] else 1 - (rib_send[i] - 1) / lens[s]
      s_angle_end <- starts[s] + s_frac_end * (ends[s] - starts[s])

      q_angles <- seq(q_angle_start, q_angle_end, length.out = 50)
      s_angles <- seq(s_angle_start, s_angle_end, length.out = 50)
      q_coords <- map_to_curve_many(q_angles, rq, q_ref)
      s_coords <- map_to_curve_many(s_angles, rs, s_ref)

      # Bezier control points
      if (!is.null(ribbon_ctrl_point)) {
        if (is.list(ribbon_ctrl_point)) {
          cp_idx <- ifelse(i > length(ribbon_ctrl_point), length(ribbon_ctrl_point), i)
          cp <- ribbon_ctrl_point[[cp_idx]]
          if (length(cp) >= 2) {
            c1 <- cp[[1]]
            c2 <- cp[[2]]
          } else {
            c1 <- c2 <- if (length(cp) == 1) cp[[1]] else c(0, 0)
          }
        } else {
          if (length(ribbon_ctrl_point) == 2) {
            c1 <- c2 <- ribbon_ctrl_point
          } else if (length(ribbon_ctrl_point) == 4) {
            c1 <- ribbon_ctrl_point[1:2]
            c2 <- ribbon_ctrl_point[3:4]
          } else {
            warning("ribbon_ctrl_point vector must have length 2 or 4; using default values")
            c1 <- c2 <- c(0, 0)
          }
        }
      } else {
        mid_angle_q <- (q_angle_start + q_angle_end) / 2
        mid_angle_s <- (s_angle_start + s_angle_end) / 2
        mid_point_q <- map_to_curve(
          mid_angle_q, seqRadius[q] + q_gap * 0.5, q_ref
        )
        mid_point_s <- map_to_curve(
          mid_angle_s, seqRadius[s] + s_gap * 0.5, s_ref
        )
        c1 <- (mid_point_q + mid_point_s) / 2
        c2 <- c1
      }

      b1 <- bezier_pts(q_coords[1, ], s_coords[1, ], c1, c1, n = 50)
      b2 <- bezier_pts(q_coords[50, ], s_coords[50, ], c2, c2, n = 50)

      ribbon_polys_list[[j]] <- cbind(
        x = c(q_coords[, 1], b2[, 1], rev(s_coords[, 1]), rev(b1[, 1])),
        y = c(q_coords[, 2], b2[, 2], rev(s_coords[, 2]), rev(b1[, 2]))
      )
      ribbon_group[j] <- j

      dir_val <- if ("direction" %in% names(ribbon_data) &&
                     !is.na(ribbon_data$direction[i]) &&
                     ribbon_data$direction[i] %in% c("same", "reverse")) {
        as.character(ribbon_data$direction[i])
      } else {
        dir_sign <- (rib_qend[i] - rib_qstart[i]) *
          (rib_send[i] - rib_sstart[i])
        ifelse(dir_sign >= 0, "same", "reverse")
      }
      ribbon_dir_vec[j] <- dir_val

      if (ribbon_color_scheme == "pident") {
        ribbon_pident[j] <- rib_pident[i]
      } else if (ribbon_color_scheme == "value") {
        ribbon_value[j] <- as.numeric(ribbon_data[[ribbon_color_by]][i])
      } else {
        ribbon_fill[j] <- switch(ribbon_color_scheme,
                                 single = singleCol,
                                 query = queryCols[q],
                                 subject = queryCols[s])
      }

      dir_alpha_factor <- if (identical(ribbon_direction, "alpha")) {
        as.numeric(ribbon_direction_alpha[dir_val])
      } else {
        1
      }
      ribbon_alpha_vec[j] <- ribbon_alpha * alpha_norm[j] * dir_alpha_factor

      if (!is.null(ribbon_outline_by)) {
        ribbon_outline_vec[j] <- unname(outline_map[as.character(ribbon_data[[ribbon_outline_by]][i])])
      } else if (identical(ribbon_direction, "outline")) {
        ribbon_outline_vec[j] <- as.character(ribbon_direction_colors[dir_val])
      }

      if (!is.null(ribbon_linetype_by)) {
        ribbon_linetype_vec[j] <- unname(linetype_map[as.character(ribbon_data[[ribbon_linetype_by]][i])])
      } else if (identical(ribbon_direction, "linetype")) {
        ribbon_linetype_vec[j] <- as.character(ribbon_direction_linetypes[dir_val])
      }
    }
    cntValid <- length(valid_idx)

    if (debug) {
      cat("Valid ribbons:", cntValid, "invalid ribbons:", cntInvalid, "
")
    }

    if (cntValid > 0) {
      m <- do.call(rbind, ribbon_polys_list)
      group_vals <- rep(ribbon_group, each = 200)
      source_row_vals <- rep(valid_idx, each = 200)
      alpha_vals <- rep(ribbon_alpha_vec, each = 200)
      outline_vals <- rep(ribbon_outline_vec, each = 200)
      linetype_vals <- rep(ribbon_linetype_vec, each = 200)
      dir_vals <- rep(ribbon_dir_vec, each = 200)
      q_gap_vals <- rep(ribbon_q_gap, each = 200)
      s_gap_vals <- rep(ribbon_s_gap, each = 200)

      if (ribbon_color_scheme == "pident") {
        ribbon_polys <- data.frame(
          x = m[, 1], y = m[, 2],
          pident = rep(ribbon_pident, each = 200),
          group = group_vals,
          source_row = source_row_vals,
          alpha = alpha_vals,
          outline_col = outline_vals,
          linetype_val = linetype_vals,
          direction = dir_vals,
          q_gap = q_gap_vals,
          s_gap = s_gap_vals,
          stringsAsFactors = FALSE
        )
      } else if (ribbon_color_scheme == "value") {
        ribbon_polys <- data.frame(
          x = m[, 1], y = m[, 2],
          value = rep(ribbon_value, each = 200),
          group = group_vals,
          source_row = source_row_vals,
          alpha = alpha_vals,
          outline_col = outline_vals,
          linetype_val = linetype_vals,
          direction = dir_vals,
          q_gap = q_gap_vals,
          s_gap = s_gap_vals,
          stringsAsFactors = FALSE
        )
      } else {
        ribbon_polys <- data.frame(
          x = m[, 1], y = m[, 2],
          fill = rep(ribbon_fill, each = 200),
          group = group_vals,
          source_row = source_row_vals,
          alpha = alpha_vals,
          outline_col = outline_vals,
          linetype_val = linetype_vals,
          direction = dir_vals,
          q_gap = q_gap_vals,
          s_gap = s_gap_vals,
          stringsAsFactors = FALSE
        )
      }
    } else {
      warning("No valid alignment data available for plotting")
    }
  }

  # ====================================================================
  # Step 6b: generate sequence-region bands and ribbon highlights (v0.9.0)
  # ====================================================================
  region_polys <- data.frame()
  if (!is.null(region_data) && nrow(region_data) > 0) {
    req <- c("seq_id", "start", "end")
    if (!all(req %in% colnames(region_data))) {
      ggchord_stop("regions must contain seq_id, start and end columns")
    }
    region_data$start <- as.numeric(region_data$start)
    region_data$end <- as.numeric(region_data$end)
    ok_rows <- is.finite(region_data$start) & is.finite(region_data$end) &
      region_data$seq_id %in% seqs
    if (any(!ok_rows)) region_data <- region_data[ok_rows, , drop = FALSE]
    if (nrow(region_data) > 0) {
      region_poly_list <- lapply(seq_len(nrow(region_data)), function(i) {
        row <- region_data[i, ]
        sid <- as.character(row$seq_id)
        len <- lens[sid]
        sp <- min(row$start, row$end)
        ep <- max(row$start, row$end)
        if (sp < 1) sp <- 1
        if (ep > len) ep <- len
        if (ep <= sp) return(NULL)

        frac_sp <- if (orientation[sid] == 1) sp / len else 1 - sp / len
        frac_ep <- if (orientation[sid] == 1) ep / len else 1 - ep / len
        a_start <- starts[sid] + frac_sp * (ends[sid] - starts[sid])
        a_end <- starts[sid] + frac_ep * (ends[sid] - starts[sid])

        ref <- seq_refs[[sid]]
        # Determine which side of the local curve normal points toward the
        # chord centre. This remains correct for non-circular curvature.
        mid_angle <- (a_start + a_end) / 2
        mid_base <- map_to_curve(mid_angle, seqRadius[sid], ref)
        mid_plus <- map_to_curve(mid_angle, seqRadius[sid] + 1e-4, ref)
        normal <- mid_plus - mid_base
        inward_sign <- if (sum(normal * -mid_base) >= 0) 1 else -1
        side <- if (identical(region_side, "auto")) "inside" else region_side
        side_sign <- if (identical(side, "inside")) inward_sign else -inward_sign
        base_radius <- seqRadius[sid] + side_sign * region_offset

        n <- 30
        angs <- seq(a_start, a_end, length.out = n)
        outer_r <- base_radius + region_width / 2
        inner_r <- base_radius - region_width / 2
        orig_ang <- c(angs, rev(angs))
        orig_rad <- c(rep(outer_r, n), rep(inner_r, n))
        mapped <- map_to_curve_many(orig_ang, orig_rad, ref)

        if (all(c("colour", "color") %in% colnames(region_data))) {
          ggchord_stop("geom_seq_region(): data may contain only one of colour and color")
        }
        colour_column <- intersect(c("colour", "color"), colnames(region_data))
        fill_col <- if (length(colour_column) && !is.na(row[[colour_column]])) {
          as.character(row[[colour_column]])
        } else {
          region_fill
        }
        data.frame(
          x = mapped[, 1], y = mapped[, 2],
          group = i,
          zregionfill = fill_col,
          colour = region_color,
          alpha = region_alpha,
          label = if ("label" %in% colnames(region_data)) as.character(row$label) else NA_character_,
          category = if ("category" %in% colnames(region_data)) as.character(row$category) else NA_character_,
          source_row = i,
          stringsAsFactors = FALSE
        )
      })
      region_poly_list <- Filter(Negate(is.null), region_poly_list)
      if (length(region_poly_list) > 0) {
        region_polys <- do.call(rbind, region_poly_list)
      }
    }
  }

  ribbon_highlight_polys <- data.frame()
  if (length(ribbon_highlight_rows) > 0 && !is.null(ribbon_polys) &&
      nrow(ribbon_polys) > 0) {
    ribbon_highlight_polys <- ribbon_polys[
      ribbon_polys$source_row %in% ribbon_highlight_rows, , drop = FALSE
    ]
  }

  # ====================================================================
  # Step 7: generate gene arrow polygons
  # ====================================================================
  gene_polys <- data.frame()
  gene_labels <- data.frame()

  if (!is.null(gene_data) && nrow(gene_data) > 0) {
    valid_gene_rows <- which(gene_data$seq_id %in% seqs)
    valid_genes <- gene_data[valid_gene_rows, , drop = FALSE]
    valid_genes$.source_row <- valid_gene_rows

    # Process gene colors
    gene_pal <- NULL
    final_gene_order <- NULL
    if (nrow(valid_genes) > 0) {
      unique_anno <- unique(valid_genes$anno)

      if (!is.null(gene_order)) {
        final_gene_order <- c(gene_order, setdiff(unique_anno, gene_order))
      } else {
        final_gene_order <- unique_anno
      }

      if (gene_color_scheme == "strand") {
        gene_pal <- process_strand_colors(gene_colors)
      } else if (gene_color_scheme == "manual") {
        gene_pal <- process_manual_colors(gene_colors, unique_anno, gene_order)
      }
    } else {
      gene_pal <- character(0)
      final_gene_order <- character(0)
    }

    # Generate feature polygons only when a geom_gene()/geom_feature() layer
    # is present. Label-only plots still use valid_genes below, without paying
    # for polygons that will never be drawn.
    gene_poly_list <- list()
    gene_rows_to_draw <- if (isTRUE(draw_gene_geometry)) {
      seq_len(nrow(valid_genes))
    } else {
      integer(0)
    }
    for (i in gene_rows_to_draw) {
      gene <- valid_genes[i, ]
      sid <- gene$seq_id
      strand <- gene$strand
      anno <- gene$anno

      width <- geneWidth[[sid]][strand]
      if (!is.numeric(width) || width <= 0) width <- 0.1

      seq_len <- lens[sid]
      sp <- min(gene$start, gene$end)
      ep <- max(gene$start, gene$end)
      if (ep <= sp) next

      frac_sp <- if (orientation[sid] == 1) sp / seq_len else 1 - sp / seq_len
      frac_ep <- if (orientation[sid] == 1) ep / seq_len else 1 - ep / seq_len
      a_start <- starts[sid] + frac_sp * (ends[sid] - starts[sid])
      a_end <- starts[sid] + frac_ep * (ends[sid] - starts[sid])
      if (strand == "-") { tmp <- a_start; a_start <- a_end; a_end <- tmp }

      r0 <- gene_track_radius(gene, sid, strand)

      ref <- seq_refs[[sid]]
      feature_shape <- if (".feature_shape" %in% names(gene)) {
        as.character(gene[[".feature_shape"]])
      } else {
        "arrow"
      }
      span <- a_end - a_start

      shape_parts <- switch(
        feature_shape,
        block = {
          ang <- seq(a_start, a_end, length.out = 60)
          list(list(
            angle = c(ang, rev(ang)),
            radius = c(rep(r0 + width / 2, length(ang)),
                       rep(r0 - width / 2, length(ang)))
          ))
        },
        chevron = {
          shoulder <- seq(a_start, a_start + 0.68 * span, length.out = 30)
          list(list(
            angle = c(shoulder, a_end, rev(shoulder),
                      a_start + 0.28 * span),
            radius = c(rep(r0 + width / 2, length(shoulder)), r0,
                       rep(r0 - width / 2, length(shoulder)), r0)
          ))
        },
        lollipop = {
          mid <- (a_start + a_end) / 2
          angle_half <- min(
            abs(span) * 0.08,
            width * 0.10 / max(abs(r0), 0.1)
          )
          stem_start <- seqRadius[sid]
          stem_end <- r0
          stem <- list(
            angle = c(mid - angle_half, mid + angle_half,
                      mid + angle_half, mid - angle_half),
            radius = c(stem_start, stem_start, stem_end, stem_end)
          )
          theta <- seq(0, 2 * pi, length.out = 48)
          head_radius <- width * 0.58
          center <- as.numeric(map_to_curve_many(mid, r0, ref)[1, ])
          delta <- max(abs(span) * 1e-4, 1e-7)
          tangent_pts <- map_to_curve_many(
            c(mid - delta, mid + delta), rep(r0, 2), ref
          )
          tangent <- as.numeric(tangent_pts[2, ] - tangent_pts[1, ])
          tangent_norm <- sqrt(sum(tangent^2))
          if (!is.finite(tangent_norm) || tangent_norm <= 1e-12) {
            tangent <- c(1, 0)
          } else {
            tangent <- tangent / tangent_norm
          }
          normal <- c(-tangent[2], tangent[1])
          head <- list(
            xy = cbind(
              center[1] + head_radius *
                (cos(theta) * tangent[1] + sin(theta) * normal[1]),
              center[2] + head_radius *
                (cos(theta) * tangent[2] + sin(theta) * normal[2])
            )
          )
          list(stem, head)
        },
        {
          n_body <- 30
          n_head <- 15
          body_ang <- seq(
            a_start, a_start + 0.6 * span, length.out = n_body
          )
          head_ang <- seq(utils::tail(body_ang, 1), a_end,
                          length.out = n_head)
          ang <- c(body_ang, head_ang)
          width_factor <- c(rep(1, n_body), seq(1, 0, length.out = n_head))
          list(list(
            angle = c(ang, rev(ang)),
            radius = c(r0 + (width / 2) * width_factor,
                       rev(r0 - (width / 2) * width_factor))
          ))
        }
      )

      for (part in seq_along(shape_parts)) {
        mapped <- if (!is.null(shape_parts[[part]]$xy)) {
          shape_parts[[part]]$xy
        } else {
          map_to_curve_many(
            shape_parts[[part]]$angle, shape_parts[[part]]$radius, ref
          )
        }
        gene_poly_list[[length(gene_poly_list) + 1]] <- data.frame(
          x = mapped[, 1],
          y = mapped[, 2],
          group = i * 10L + part,
          anno = anno,
          strand = strand,
          feature_shape = feature_shape,
          source_row = gene$.source_row,
          ord = seq_len(nrow(mapped)),
          stringsAsFactors = FALSE
        )
      }
    }
    gene_polys <- if (length(gene_poly_list)) do.call(rbind, gene_poly_list) else data.frame()

    # Generate gene labels
    if (gene_label_show && nrow(valid_genes) > 0) {
      gene_labels <- do.call(rbind, lapply(seq_len(nrow(valid_genes)), function(i) {
        gene <- valid_genes[i, ]
        sid <- gene$seq_id
        strand <- gene$strand
        seq_len <- lens[sid]
        ref <- seq_refs[[sid]]
        orient <- orientation[sid]

        sp <- min(gene$start, gene$end)
        ep <- max(gene$start, gene$end)
        frac_mid <- (sp + ep) / (2 * seq_len)

        circum_ratio <- geneLabelCircumOffset[[sid]][strand]
        if (geneLabelCircumLimit[[sid]][strand]) {
          gene_length_ratio <- (ep - sp) / seq_len
          max_offset_ratio <- gene_length_ratio * 0.5
          circum_ratio <- pmin(max_offset_ratio, pmax(-max_offset_ratio, circum_ratio))
        }
        frac_mid <- frac_mid + circum_ratio
        frac_mid <- pmin(1, pmax(0, frac_mid))

        if (orient != 1) frac_mid <- 1 - frac_mid

        ref_n <- length(ref$angles)
        idx <- round(frac_mid * (ref_n - 1)) + 1
        idx <- pmin(ref_n, pmax(1, idx))

        if (idx < ref_n) {
          dx <- ref$path$x[idx + 1] - ref$path$x[idx]
          dy <- ref$path$y[idx + 1] - ref$path$y[idx]
        } else {
          dx <- ref$path$x[idx] - ref$path$x[idx - 1]
          dy <- ref$path$y[idx] - ref$path$y[idx - 1]
        }
        dx <- dx * orient
        dy <- dy * orient

        width <- geneWidth[[sid]][strand]

        r0 <- gene_track_radius(gene, sid, strand)

        center_r <- r0
        center_pt <- map_to_curve(angle = ref$angles[idx], radius = center_r, ref = ref)

        normal_x <- -dy
        normal_y <- dx
        nl <- sqrt(normal_x^2 + normal_y^2)
        if (nl > 0) {
          normal_x <- normal_x / nl
          normal_y <- normal_y / nl
        }

        direction_factor <- ifelse(strand == "+", 1, -1) * orient
        normal_x <- normal_x * direction_factor
        normal_y <- normal_y * direction_factor

        text_x <- center_pt[1] - normal_x * geneLabelRadialOffset[[sid]][strand]
        text_y <- center_pt[2] - normal_y * geneLabelRadialOffset[[sid]][strand]
        # Leader-line origin: the fixed label position next to the gene. When
        # the label is moved to the other side of its arc, the line still
        # starts here (at the gene) and only the repelled text position moves.
        anchor_x <- text_x
        anchor_y <- text_y

        # Optional side flip: mirror labels across their sequence arc (e.g.
        # inner labels to the outside so they do not overlap the ribbons).
        # map_to_curve() uses a radius parameter R where R < seqRadius places
        # points outside the chord and R > seqRadius inside; reflecting across
        # the arc means R' = 2 * seqRadius - R, which preserves the label's
        # distance from the arc.
        side_flipped <- FALSE
        R_label <- r0 - direction_factor * geneLabelRadialOffset[[sid]][strand]
        if ((identical(gene_label_side, "outside") && R_label > seqRadius[sid]) ||
            (identical(gene_label_side, "inside") && R_label < seqRadius[sid])) {
          dR <- 2 * seqRadius[sid] - 2 * R_label
          # text = base + norm * (R - r0) with norm the direction-adjusted
          # normal; the unadjusted normal is normal / direction_factor.
          text_x <- text_x + (normal_x / direction_factor) * dR
          text_y <- text_y + (normal_y / direction_factor) * dR
          side_flipped <- TRUE
        }

        base_angle <- atan2(dy, dx) * 180 / pi
        text_angle <- switch(
          gene_label_orientation,
          radial = base_angle + 90,
          tangent = base_angle,
          # coord_chord() subsequently rotates every grob by `rotation`;
          # compensate here so "horizontal" means horizontal on the device.
          horizontal = -rotation
        ) + geneLabelRotation[[sid]][strand]

        if (strand == "+" && orient == 1) {
          hjust <- 1
        } else if (strand == "+" && orient != 1) {
          hjust <- 0
        } else if (strand == "-" && orient == 1) {
          hjust <- 0
        } else {
          hjust <- 1
        }

        text_angle <- (text_angle + 360) %% 360
        if (text_angle > 90 && text_angle < 270) {
          text_angle <- text_angle + 180
          hjust <- 1 - hjust
        }
        text_angle <- text_angle %% 360
        vjust <- 0.5

        if (identical(gene_label_orientation, "tangent")) {
          hjust <- 0.5
        } else if (identical(gene_label_orientation, "horizontal")) {
          rotation_rad <- rotation * pi / 180
          # Anchor text by the direction in which it was actually displaced
          # from its own sequence curve. Using its position relative to the
          # global origin made an inside label extend back across the arc in
          # some quadrants, especially with gene_label_side = "auto".
          label_dx <- text_x - center_pt[1]
          label_dy <- text_y - center_pt[2]
          device_x <- cos(rotation_rad) * label_dx -
            sin(rotation_rad) * label_dy
          device_y <- sin(rotation_rad) * label_dx +
            cos(rotation_rad) * label_dy
          # Prefer a left/right anchor in diagonal quadrants so horizontal
          # text extends away from the chord rather than half back across it.
          if (abs(device_x) >= 0.75 * abs(device_y)) {
            hjust <- if (device_x >= 0) 0 else 1
            vjust <- 0.5
          } else {
            hjust <- 0.5
            vjust <- if (device_y >= 0) 0 else 1
          }
        }

        data.frame(
          text = gene$anno,
          text_x = text_x,
          text_y = text_y,
          text_angle = text_angle,
          hjust = hjust,
          vjust = vjust,
          size = gene_label_size,
          seq_id = sid,
          group = i,
          source_row = gene$.source_row,
          anchor_x = anchor_x,
          anchor_y = anchor_y,
          side_flipped = side_flipped,
          stringsAsFactors = FALSE
        )
      }))
    }
  } else {
    gene_pal <- character(0)
    final_gene_order <- character(0)
  }

  # ====================================================================
  # Step 7b: generate sequence labels (if requested)
  # ====================================================================
  seq_labels_df <- data.frame()
  if (!is.null(seq_label_text)) {
    seq_labels_df <- do.call(rbind, lapply(seqs, function(id) {
      ref <- seq_refs[[id]]
      mid_angle <- (starts[id] + ends[id]) / 2
      # seq_label_radius is a multiplier of the arc radius: 1 = on the arc,
      # > 1 = outside (away from the chord center), < 1 = inside. map_to_curve()
      # measures its radius parameter along the inward normal, so the
      # multiplier must be mirrored: R = seqRadius * (2 - multiplier).
      r <- seqRadius[id] * (2 - seq_label_radius[id])
      pt <- map_to_curve(mid_angle, r, ref)
      # Tangent angle at the midpoint, used to orient the label along the arc.
      idx <- which.min(abs(ref$angles - mid_angle))
      if (idx < length(ref$angles)) {
        dx <- ref$path$x[idx + 1] - ref$path$x[idx]
        dy <- ref$path$y[idx + 1] - ref$path$y[idx]
      } else {
        dx <- ref$path$x[idx] - ref$path$x[idx - 1]
        dy <- ref$path$y[idx] - ref$path$y[idx - 1]
      }
      text_angle <- atan2(dy, dx) * 180 / pi + 90 + seq_label_rotation[id]
      hjust <- if (is.null(seq_label_hjust)) 0.5 else seq_label_hjust[[id]]
      vjust <- if (is.null(seq_label_vjust)) 0.5 else seq_label_vjust[[id]]
      text_angle <- (text_angle + 360) %% 360
      if (text_angle > 90 && text_angle < 270) {
        text_angle <- text_angle + 180
        # keep the text box anchored when a user-supplied hjust is flipped
        if (!is.null(seq_label_hjust)) hjust <- 1 - hjust
      }
      text_angle <- text_angle %% 360
      data.frame(
        text_x = pt[1], text_y = pt[2],
        label = seq_label_text[id],
        text_angle = text_angle,
        size = seq_label_size[id],
        hjust = hjust, vjust = vjust,
        seq_id = id,
        stringsAsFactors = FALSE
      )
    }))
  }

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
  if (nrow(ribbon_highlight_polys) > 0) ribbon_highlight_polys <- rotate_df(ribbon_highlight_polys)
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
      layout_box_padding <- if (identical(gene_label_layout, "auto")) {
        # Physical inches on each side of a text box. The previous 0.18-inch
        # value made an eight-label vertical rail more than twice as tall as
        # necessary and forced top/bottom labels onto extra rows. About 1 mm
        # keeps labels visually separate while letting the cardinal rails use
        # the available horizontal and vertical perimeter efficiently.
        0.05
      } else {
        0.04
      }
      layout_point_padding <- if (identical(gene_label_layout, "auto")) {
        0.08
      } else {
        0.05
      }
      layout_min_segment <- 0.02
      layout_units <- text_units_per_inch
      base_gene_labels <- gene_labels
      layout_result <- NULL

      if (is.null(gene_label_wrap) &&
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
        if (identical(gene_label_layout, "auto")) {
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
      gene_labels <- ggchord_hide_conflicted_labels(
        gene_labels,
        max_overlaps = gene_label_repel_max_overlaps,
        units_per_inch = layout_units,
        repel_boxes = final_obstacles
      )
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
          include_own = identical(gene_label_layout, "arc"),
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

  # ====================================================================
  # Step 8c: optionally hide axis labels that overlap other elements
  # ====================================================================
  if (isTRUE(axis_label_hide_overlaps) && nrow(axis_ticks) > 0 &&
      any(!is.na(axis_ticks$label))) {
    content <- ggchord_repel_points(seq_arcs, gene_polys,
                                    data.frame(x = numeric(0), y = numeric(0)),
                                    data.frame(), show_axis = FALSE)
    axis_ticks <- ggchord_hide_text_overlaps(
      axis_ticks, content,
      units_per_inch = text_units_per_inch
    )
  }

  # ====================================================================
  # Step 9: compute plot extremes
  # ====================================================================
  extremes <- get_plot_extremes(
    allRibbon = ribbon_polys,
    seqArcs = seq_arcs,
    axisLines = axis_lines,
    axisTicks = axis_ticks,
    gene_polys = gene_polys,
    gene_arrows = gene_labels,
    seq_labels = seq_labels_df,
    show_axis = show_axis
  )

  # ====================================================================
  # Step 10: assemble and return the layout object
  # ====================================================================
  layout <- list(
    # Geometric data
    seq_arcs       = seq_arcs,
    ribbon_polys   = ribbon_polys,
    region_polys   = region_polys,
    ribbon_highlight_polys = ribbon_highlight_polys,
    gene_polys     = gene_polys,
    gene_labels    = gene_labels,
    gene_label_segments = gene_label_segments,
    gene_label_clip_units = gene_label_clip_units,
    text_units_per_inch = text_units_per_inch,
    gene_label_layout = gene_label_layout,
    feature_shape_pal = feature_shape_pal,
    feature_shape_order = feature_shape_order,
    seq_labels_df  = seq_labels_df,
    axis_lines     = axis_lines,
    axis_ticks     = axis_ticks,

    # Extremes
    extremes       = extremes,

    # Colors and labels
    seq_colors     = seq_colors,
    seq_labels     = seq_labels,
    seqs           = seqs,
    seqRadius      = seqRadius,

    # Ribbon-related
    ribbon_color_scheme = ribbon_color_scheme,
    ribbon_colors  = ribbon_colors,
    ribbon_alpha   = ribbon_alpha,
    ribbon_color_by = ribbon_color_by,
    ribbon_color_limits = ribbon_color_limits,
    ribbon_color_breaks = ribbon_color_breaks,
    ribbon_color_name = ribbon_color_name,
    ribbon_use_outline = ribbon_use_outline,
    ribbon_use_linetype = ribbon_use_linetype,

    # Gene-related
    gene_pal           = gene_pal,
    gene_color_scheme  = gene_color_scheme,
    final_gene_order   = final_gene_order,
    gene_label_show    = gene_label_show,
    gene_label_size    = gene_label_size,

    # Axis-related
    show_axis           = show_axis,
    axisLabelOrientation = axisLabelOrientation,

    # Metadata
    rotation        = rotation,
    n_sequences     = n
  )

  class(layout) <- "chord_layout"
  layout
}
