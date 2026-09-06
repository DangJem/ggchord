# Alignment ribbon geometry and mapped ribbon metadata.
ggchord_layout_ribbon_step <- quote({
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
    rib_pident <- ribbon_data$pident %||% rep(NA_real_, nrow(ribbon_data))

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
    ribbon_q_n <- ribbon_s_n <- integer(length(valid_idx))

    # With the default NULL ribbon_gap, determine spacing independently at
    # each ribbon endpoint. Only polygons that occupy the ribbon-facing side
    # of the sequence and overlap that endpoint's genomic interval count as
    # obstacles. Text and leader segments never enter this table.
    endpoint_ribbon_gap <- function(accver, start, end) {
      configured <- unname(ribbonGap[[accver]])
      if (!isTRUE(ribbon_gap_auto)) return(configured)
      max(ggchord_gap_profile(c(start, end), accver, ribbon_obstacles,
        if (link_avoid == "none") "none" else "uniform", min(configured, .035)))
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

    front_fallbacks <- 0L
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

      make_front <- function(id, start, end, ref, a0, a1, gap) {
        positions <- seq(start, end, length.out = 50)
        if (ribbon_gap_auto && link_avoid == "smooth" && nrow(ribbon_obstacles)) {
          obs <- ribbon_obstacles[ribbon_obstacles$accver == id,,drop=FALSE]
          span <- max(abs(end-start)*.08, 1)
          knots <- c(obs$start-span, obs$start, obs$end, obs$end+span)
          knots <- knots[knots > min(start,end) & knots < max(start,end)]
          positions <- sort(unique(c(positions, knots)), decreasing = end < start)
        }
        front <- function(pos) {
          angles <- if (end == start) rep(a0, length(pos)) else a0 + (pos-start)/(end-start)*(a1-a0)
          gaps <- if (ribbon_gap_auto) ggchord_gap_profile(pos, id, ribbon_obstacles,
            link_avoid, min(unname(ribbonGap[[id]]), .035)) else rep(gap,length(pos))
          xy <- map_to_curve_many(angles, seqRadius[id] + gaps, ref)
          list(xy=xy, baseline=map_to_curve_many(angles,seqRadius[id],ref), gaps=gaps)
        }
        result <- front(positions)
        if (ribbon_gap_auto && link_avoid == "smooth" && ggchord_front_invalid(result$xy,result$baseline)) {
          positions <- sort(unique(c(positions, (head(positions,-1)+tail(positions,-1))/2)), decreasing=end<start)
          result <- front(positions)
          if (ggchord_front_invalid(result$xy,result$baseline)) {
            angles <- if (end == start) rep(a0, length(positions)) else a0+(positions-start)/(end-start)*(a1-a0)
            result$xy <- map_to_curve_many(angles,seqRadius[id]+gap,ref)
            front_fallbacks <<- front_fallbacks + 1L
          }
        }
        result$xy
      }
      q_coords <- make_front(q, rib_qstart[i], rib_qend[i], q_ref, q_angle_start, q_angle_end, q_gap)
      s_coords <- make_front(s, rib_sstart[i], rib_send[i], s_ref, s_angle_start, s_angle_end, s_gap)

      ribbon_q_n[j] <- nrow(q_coords)
      ribbon_s_n[j] <- nrow(s_coords)
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
      b2 <- bezier_pts(q_coords[nrow(q_coords), ], s_coords[nrow(s_coords), ], c2, c2, n = 50)

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
      polygon_sizes <- vapply(ribbon_polys_list, nrow, integer(1))
      group_vals <- rep(ribbon_group, times = polygon_sizes)
      source_row_vals <- rep(valid_idx, times = polygon_sizes)
      alpha_vals <- rep(ribbon_alpha_vec, times = polygon_sizes)
      outline_vals <- rep(ribbon_outline_vec, times = polygon_sizes)
      linetype_vals <- rep(ribbon_linetype_vec, times = polygon_sizes)
      dir_vals <- rep(ribbon_dir_vec, times = polygon_sizes)
      q_gap_vals <- rep(ribbon_q_gap, times = polygon_sizes)
      s_gap_vals <- rep(ribbon_s_gap, times = polygon_sizes)

      if (ribbon_color_scheme == "pident") {
        if (front_fallbacks > 0L) warning(sprintf("Smooth avoidance used uniform fallback at %d endpoint(s)", front_fallbacks), call. = FALSE)
    ribbon_polys <- data.frame(
          x = m[, 1], y = m[, 2],
          pident = rep(ribbon_pident, times = polygon_sizes),
          group = group_vals,
          source_row = source_row_vals,
          alpha = alpha_vals,
          outline_col = outline_vals,
          linetype_val = linetype_vals,
          direction = dir_vals,
          .q_n = rep(ribbon_q_n, polygon_sizes), .s_n = rep(ribbon_s_n, polygon_sizes),
          q_gap = q_gap_vals,
          s_gap = s_gap_vals,
          stringsAsFactors = FALSE
        )
      } else if (ribbon_color_scheme == "value") {
        ribbon_polys <- data.frame(
          x = m[, 1], y = m[, 2],
          value = rep(ribbon_value, times = polygon_sizes),
          group = group_vals,
          source_row = source_row_vals,
          alpha = alpha_vals,
          outline_col = outline_vals,
          linetype_val = linetype_vals,
          direction = dir_vals,
          .q_n = rep(ribbon_q_n, polygon_sizes), .s_n = rep(ribbon_s_n, polygon_sizes),
          q_gap = q_gap_vals,
          s_gap = s_gap_vals,
          stringsAsFactors = FALSE
        )
      } else {
        ribbon_polys <- data.frame(
          x = m[, 1], y = m[, 2],
          fill = rep(ribbon_fill, times = polygon_sizes),
          group = group_vals,
          source_row = source_row_vals,
          alpha = alpha_vals,
          outline_col = outline_vals,
          linetype_val = linetype_vals,
          direction = dir_vals,
          .q_n = rep(ribbon_q_n, polygon_sizes), .s_n = rep(ribbon_s_n, polygon_sizes),
          q_gap = q_gap_vals,
          s_gap = s_gap_vals,
          stringsAsFactors = FALSE
        )
      }
    } else {
      warning("No valid alignment data available for plotting")
    }
  }

})
