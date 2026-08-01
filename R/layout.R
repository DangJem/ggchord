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
#' @param ribbon_data Alignment data (already validated)
#' @param gene_data Gene data (already validated)
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
    ribbon_color_scheme, ribbon_colors, ribbon_alpha,
    ribbon_ctrl_point,
    # Gene parameters
    gene_data = NULL,
    geneGap, geneWidth,
    geneLabelRadialOffset, geneLabelCircumOffset,
    geneLabelCircumLimit, geneLabelRotation,
    gene_label_show, gene_label_size,
    gene_color_scheme, gene_colors, gene_order,
    # Sequence label parameters
    seq_label_text = NULL, seq_label_radius = NULL,
    seq_label_rotation = NULL, seq_label_size = NULL,
    # Axis parameters
    axisGap, axisMaj, axisMajLen, axisMin, axisMinLen,
    labelSize, labelOffset, axisLabelOrientation,
    show_axis,
    # Global parameters
    rotation, debug = FALSE
) {
  n <- length(seqs)

  # ====================================================================
  # Step 1: compute angle allocation
  # ====================================================================
  total_circ <- 2 * pi
  total_gap_prop <- sum(seq_gap)

  if (total_gap_prop >= 1) {
    stop("The sum of seq_gap cannot exceed 1 (no space left for sequences)")
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
  seq_refs <- lapply(seqs, function(id) {
    path <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id], seq_curvature[id], n_points = nRef
    )
    angles <- seq(starts[id], ends[id], length.out = nRef)
    list(path = path, angles = angles, r0 = seqRadius[id])
  })
  names(seq_refs) <- seqs

  # Curve coordinate mapping function
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

  # ====================================================================
  # Step 4: generate sequence arcs (outer layer)
  # ====================================================================
  seq_arcs <- lapply(seqs, function(id) {
    path_data <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id], seq_curvature[id], nSeg
    )
    path_data$seq_id <- id
    if (orientation[id] == -1) {
      path_data <- path_data[nrow(path_data):1, ]
    }
    path_data
  })

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
                           seq_id = character(0),
                           stringsAsFactors = FALSE)

  if (show_axis) {
    # Axis lines
    axis_lines <- do.call(rbind, lapply(seqs, function(id) {
      ref <- seq_refs[[id]]
      r0 <- ref$r0 - axisGap[id]
      angles <- seq(starts[id], ends[id], length.out = nSeg)
      pts <- t(sapply(angles, function(angle) {
        map_to_curve(angle, r0, ref)
      }))
      data.frame(x = pts[, 1], y = pts[, 2], seq_id = id, stringsAsFactors = FALSE)
    }))

    # Tick marks
    axis_ticks <- do.call(rbind, lapply(seqs, function(id) {
      ref <- seq_refs[[id]]
      r0 <- ref$r0 - axisGap[id]

      majors <- breakPointsFunc(lens[id], axisMaj[id])
      minors <- unlist(lapply(seq_len(length(majors) - 1), function(i) {
        seq(majors[i], majors[i + 1], length.out = axisMin[id] + 2)[-c(1, axisMin[id] + 2)]
      }))
      pts <- data.frame(
        pos = c(majors, minors),
        is_major = c(rep(TRUE, length(majors)), rep(FALSE, length(minors)))
      )

      do.call(rbind, lapply(seq_len(nrow(pts)), function(j) {
        p <- pts[j, ]
        frac <- if (orientation[id] == 1) p$pos / lens[id] else 1 - p$pos / lens[id]
        angle <- starts[id] + frac * (ends[id] - starts[id])

        base <- map_to_curve(angle, r0, ref)
        dir <- if (axisGap[id] >= 0) -1 else 1
        len <- if (p$is_major) axisMajLen[id] else axisMinLen[id]
        tip <- map_to_curve(angle, r0 + len * dir, ref)
        lbl <- map_to_curve(angle, r0 + len * (1.5 + labelOffset[id]) * dir, ref)

        data.frame(
          x0 = base[1], y0 = base[2],
          x1 = tip[1], y1 = tip[2],
          label = if (p$is_major) as.character(p$pos) else NA,
          label_x = lbl[1], label_y = lbl[2],
          size = labelSize[id],
          seq_id = id,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }

  # ====================================================================
  # Step 6: generate ribbon polygons
  # ====================================================================
  ribbon_polys <- NULL
  ribbon_color_info <- list(scheme = ribbon_color_scheme)

  if (!is.null(ribbon_data) && nrow(ribbon_data) > 0) {
    ribbons <- list()
    cntValid <- 0
    cntInvalid <- 0

    # Pre-process ribbon colors
    singleCol <- NULL
    queryCols <- NULL
    rampFunc <- NULL

    if (ribbon_color_scheme == "single") {
      singleCol <- if (length(ribbon_colors) > 1) ribbon_colors[[1]] else ribbon_colors
    } else if (ribbon_color_scheme %in% c("query", "subject")) {
      queryCols <- ribbon_colors
    } else if (ribbon_color_scheme == "pident") {
      rampFunc <- colorRampPalette(ribbon_colors)
    }

    for (i in seq_len(nrow(ribbon_data))) {
      row <- ribbon_data[i, ]
      q <- row$qaccver
      s <- row$saccver
      if (q == s || !q %in% seqs || !s %in% seqs) {
        cntInvalid <- cntInvalid + 1
        next
      }

      # Query sequence coordinates
      q_ref <- seq_refs[[q]]
      q_frac_start <- if (orientation[q] == 1) (row$qstart - 1) / lens[q] else 1 - (row$qstart - 1) / lens[q]
      q_angle_start <- starts[q] + q_frac_start * (ends[q] - starts[q])
      q_start_coord <- map_to_curve(q_angle_start, seqRadius[q] + ribbonGap[q], q_ref)

      q_frac_end <- if (orientation[q] == 1) (row$qend - 1) / lens[q] else 1 - (row$qend - 1) / lens[q]
      q_angle_end <- starts[q] + q_frac_end * (ends[q] - starts[q])
      q_end_coord <- map_to_curve(q_angle_end, seqRadius[q] + ribbonGap[q], q_ref)

      q_angles <- seq(q_angle_start, q_angle_end, length.out = 50)
      q_coords <- do.call(rbind, lapply(q_angles, function(angle) {
        map_to_curve(angle, seqRadius[q] + ribbonGap[q], q_ref)
      }))
      segQ <- data.frame(x = q_coords[, 1], y = q_coords[, 2])

      # Subject sequence coordinates
      s_ref <- seq_refs[[s]]
      s_frac_start <- if (orientation[s] == 1) (row$sstart - 1) / lens[s] else 1 - (row$sstart - 1) / lens[s]
      s_angle_start <- starts[s] + s_frac_start * (ends[s] - starts[s])
      s_start_coord <- map_to_curve(s_angle_start, seqRadius[s] + ribbonGap[s], s_ref)

      s_frac_end <- if (orientation[s] == 1) (row$send - 1) / lens[s] else 1 - (row$send - 1) / lens[s]
      s_angle_end <- starts[s] + s_frac_end * (ends[s] - starts[s])
      s_end_coord <- map_to_curve(s_angle_end, seqRadius[s] + ribbonGap[s], s_ref)

      s_angles <- seq(s_angle_start, s_angle_end, length.out = 50)
      s_coords <- do.call(rbind, lapply(s_angles, function(angle) {
        map_to_curve(angle, seqRadius[s] + ribbonGap[s], s_ref)
      }))
      segS <- data.frame(x = s_coords[, 1], y = s_coords[, 2])

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
        mid_point_q <- map_to_curve(mid_angle_q, seqRadius[q] + ribbonGap[q] * 0.5, q_ref)
        mid_point_s <- map_to_curve(mid_angle_s, seqRadius[s] + ribbonGap[s] * 0.5, s_ref)
        c1 <- colMeans(rbind(mid_point_q, mid_point_s))
        c2 <- c1
      }

      # Generate Bezier curves
      b1 <- bezier_pts(
        as.numeric(segQ[1, ]),
        as.numeric(segS[1, ]),
        c1, c1,
        n = 50
      )
      b2 <- bezier_pts(
        as.numeric(segQ[nrow(segQ), ]),
        as.numeric(segS[nrow(segS), ]),
        c2, c2,
        n = 50
      )

      # Close the polygon
      poly <- rbind(
        segQ,
        b2,
        segS[nrow(segS):1, ],
        b1[nrow(b1):1, ]
      )

      if (ribbon_color_scheme == "pident") {
        poly$pident <- row$pident
      } else {
        poly$fill <- switch(ribbon_color_scheme,
                            single = singleCol,
                            query = queryCols[q],
                            subject = queryCols[s])
      }
      poly$group <- cntValid + 1
      poly$alpha <- ribbon_alpha
      ribbons[[length(ribbons) + 1]] <- poly
      cntValid <- cntValid + 1
    }

    if (debug) {
      cat("Valid ribbons:", cntValid, "invalid ribbons:", cntInvalid, "\n")
    }

    if (cntValid > 0) {
      ribbon_polys <- do.call(rbind, ribbons)
    } else {
      warning("No valid alignment data available for plotting")
    }
  }

  # ====================================================================
  # Step 7: generate gene arrow polygons
  # ====================================================================
  gene_polys <- data.frame()
  gene_labels <- data.frame()

  if (!is.null(gene_data) && nrow(gene_data) > 0) {
    valid_genes <- gene_data[gene_data$seq_id %in% seqs, ]

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

    # Generate arrow polygons
    for (i in seq_len(nrow(valid_genes))) {
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

      n_body <- 30
      n_head <- 15
      body_ang <- seq(a_start, a_start + 0.6 * (a_end - a_start), length.out = n_body)
      head_ang <- seq(tail(body_ang, 1), a_end, length.out = n_head)
      angs <- c(body_ang, head_ang)
      total_pt <- length(angs)
      widths <- c(rep(1, n_body), seq(1, 0, length.out = n_head))

      if (strand == "+") {
        r0 <- seqRadius[sid] - geneGap[[sid]][strand]
      } else {
        r0 <- seqRadius[sid] + geneGap[[sid]][strand]
      }

      outer_r <- r0 + (width / 2) * widths
      inner_r <- r0 - (width / 2) * widths

      orig_ang <- c(angs, rev(angs))
      orig_rad <- c(outer_r, rev(inner_r))

      ref <- seq_refs[[sid]]
      pts_list <- mapply(function(ang, rad) {
        map_to_curve(ang, rad, ref)
      }, orig_ang, orig_rad, SIMPLIFY = FALSE)

      mapped <- do.call(rbind, pts_list)

      gene_poly <- data.frame(
        x = mapped[, 1],
        y = mapped[, 2],
        group = i,
        anno = anno,
        strand = strand,
        ord = seq_len(2 * total_pt),
        stringsAsFactors = FALSE
      )
      gene_polys <- rbind(gene_polys, gene_poly)
    }

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

        if (strand == "+") {
          r0 <- seqRadius[sid] - geneGap[[sid]][strand]
        } else {
          r0 <- seqRadius[sid] + geneGap[[sid]][strand]
        }

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

        base_angle <- atan2(dy, dx) * 180 / pi
        text_angle <- base_angle + 90 + geneLabelRotation[[sid]][strand]

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

        data.frame(
          text = gene$anno,
          text_x = text_x,
          text_y = text_y,
          text_angle = text_angle,
          hjust = hjust,
          vjust = 0.5,
          seq_id = sid,
          group = i,
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
      r <- seqRadius[id] * seq_label_radius[id]
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
      text_angle <- (text_angle + 360) %% 360
      if (text_angle > 90 && text_angle < 270) {
        text_angle <- text_angle + 180
      }
      text_angle <- text_angle %% 360
      data.frame(
        text_x = pt[1], text_y = pt[2],
        label = seq_label_text[id],
        text_angle = text_angle,
        size = seq_label_size[id],
        hjust = 0.5, vjust = 0.5,
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
    if (all(c("text_x", "text_y") %in% names(df))) {
      TX <- df$text_x; TY <- df$text_y
      df$text_x <- TX * cos(rot_rad) - TY * sin(rot_rad)
      df$text_y <- TX * sin(rot_rad) + TY * cos(rot_rad)
      df$text_angle <- df$text_angle + rotation
    }
    df
  }

  # Apply to all geometric elements (avoid deep copies: modify the original object's columns in place)
  seq_arcs <- lapply(seq_arcs, rotate_df)
  if (nrow(axis_lines) > 0) axis_lines <- rotate_df(axis_lines)
  if (nrow(axis_ticks) > 0) axis_ticks <- rotate_df(axis_ticks)
  if (!is.null(ribbon_polys)) ribbon_polys <- rotate_df(ribbon_polys)
  if (nrow(gene_labels) > 0) gene_labels <- rotate_df(gene_labels)
  if (nrow(seq_labels_df) > 0) seq_labels_df <- rotate_df(seq_labels_df)
  if (nrow(gene_polys) > 0) {
    gene_polys <- rotate_df(gene_polys)
    gene_polys <- gene_polys[with(gene_polys, order(group, ord)), ]
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
    gene_polys     = gene_polys,
    gene_labels    = gene_labels,
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
