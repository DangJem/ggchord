# Region, gene and fixed sequence-label geometry.
ggchord_layout_annotation_step <- quote({
  # ====================================================================
  # Step 6b: generate sequence-region bands
  # ====================================================================
  region_polys <- data.frame()
  if (!is.null(region_data) && nrow(region_data) > 0) {
    req <- c("accver", "start", "end")
    if (!all(req %in% colnames(region_data))) {
      ggchord_stop("regions must contain accver, start and end columns")
    }
    region_data$start <- as.numeric(region_data$start)
    region_data$end <- as.numeric(region_data$end)
    ok_rows <- is.finite(region_data$start) & is.finite(region_data$end) &
      region_data$accver %in% seqs
    if (any(!ok_rows)) region_data <- region_data[ok_rows, , drop = FALSE]
    if (nrow(region_data) > 0) {
      region_poly_list <- lapply(seq_len(nrow(region_data)), function(i) {
        row <- region_data[i, ]
        sid <- as.character(row$accver)
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

  # ====================================================================
  # Step 7: generate gene arrow polygons
  # ====================================================================
  gene_polys <- data.frame()
  gene_labels <- data.frame()

  if (!is.null(gene_data) && nrow(gene_data) > 0) {
    valid_gene_rows <- which(gene_data$accver %in% seqs)
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
      sid <- gene$accver
      strand <- gene$strand
      anno <- gene$anno

      width <- geneWidth[[sid]][strand]
      width <- width * as.numeric(gene$.feature_width_factor %||% 1)
      if (!is.numeric(width) || width <= 0) width <- 0.1

      sequence_length <- lens[sid]
      if (!is.finite(gene$start) || !is.finite(gene$end) ||
          gene$start == gene$end) next
      ref <- seq_refs[[sid]]
      feature_shape <- if (".feature_shape" %in% names(gene)) {
        as.character(gene[[".feature_shape"]])
      } else {
        "arrow"
      }
      pieces <- ggchord_feature_intervals(
        gene$start, gene$end, sequence_length, strand,
        circular = circular
      )
      part_index <- 0L
      for (piece in pieces) {
        interval <- ggchord_feature_angle_interval(
          piece, sequence_length, orientation[sid], starts[sid], ends[sid]
        )
        r0 <- gene_track_radius(
          gene, sid, strand, mean(interval)
        )
        shape_parts <- ggchord_feature_geometry(
          feature_shape, interval[1L], interval[2L], r0, width,
          seqRadius[sid], ref,
          arrow_head_length = arrow_head_length,
          arrow_head_width = arrow_head_width,
          arrow_head_style = arrow_head_style,
          short_feature = short_feature,
          draw_head = piece$draw_head,
          bidirectional = identical(
            as.character(gene$.feature_biological_strand %||% strand), "+/-"
          )
        )
        for (part in seq_along(shape_parts)) {
          part_index <- part_index + 1L
          mapped <- if (!is.null(shape_parts[[part]]$xy)) {
            shape_parts[[part]]$xy
          } else {
            map_to_curve_many(
              shape_parts[[part]]$angle, shape_parts[[part]]$radius, ref
            )
          }
          gene_poly_list[[length(gene_poly_list) + 1L]] <- data.frame(
          x = mapped[, 1],
          y = mapped[, 2],
          group = i * 100L + part_index,
          anno = anno,
          strand = strand,
          feature_shape = feature_shape,
          biological_strand = as.character(
            gene$.feature_biological_strand %||% strand
          ),
          feature_fill_explicit = if ("feature_color" %in% names(gene)) {
            as.character(gene$feature_color)
          } else NA_character_,
          position_name = as.character(gene$.position_name %||% "identity"),
          base_offset = as.numeric(gene$.position_base_offset %||% 0),
          lane = as.integer(gene$.feature_stack_lane %||% 0L),
          lane_offset = as.numeric(gene$.position_lane_offset %||% 0),
          normal_offset = as.numeric(gene$.normal_offset %||% 0),
          source_row = gene$.source_row,
          ord = seq_len(nrow(mapped)),
          .component = "polygon",
          stringsAsFactors = FALSE
        )
        }

        boundaries <- if (".feature_boundaries" %in% names(gene)) {
          gene$.feature_boundaries[[1L]]
        } else numeric()
        boundary_styles <- if (".feature_boundary_styles" %in% names(gene)) {
          as.character(gene$.feature_boundary_styles[[1L]])
        } else rep(NA_character_, length(boundaries))
        if (length(boundaries) && length(pieces) == 1L) {
          keep_boundaries <-
            boundaries > min(gene$start, gene$end) &
              boundaries < max(gene$start, gene$end)
          boundaries <- boundaries[keep_boundaries]
          boundary_styles <- boundary_styles[keep_boundaries]
          for (boundary_index in seq_along(boundaries)) {
            boundary <- boundaries[boundary_index]
            boundary_piece <- list(start = boundary, end = boundary)
            boundary_angle <- ggchord_feature_angle_interval(
              boundary_piece, sequence_length, orientation[sid],
              starts[sid], ends[sid]
            )[1L]
            boundary_xy <- map_to_curve_many(
              rep(boundary_angle, 2L),
              c(r0 - width / 2, r0 + width / 2), ref
            )
            part_index <- part_index + 1L
            gene_poly_list[[length(gene_poly_list) + 1L]] <- data.frame(
              x = boundary_xy[, 1L], y = boundary_xy[, 2L],
              group = i * 100L + part_index, anno = anno,
              strand = strand, feature_shape = feature_shape,
              position_name = as.character(gene$.position_name %||% "identity"),
              base_offset = as.numeric(gene$.position_base_offset %||% 0),
              lane = as.integer(gene$.feature_stack_lane %||% 0L),
              lane_offset = as.numeric(gene$.position_lane_offset %||% 0),
              normal_offset = as.numeric(gene$.normal_offset %||% 0),
              boundary_linetype = boundary_styles[boundary_index],
              source_row = gene$.source_row, ord = seq_len(nrow(boundary_xy)),
              .component = "boundary", stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    gene_polys <- if (length(gene_poly_list)) {
      ggchord_rbind_fill(gene_poly_list)
    } else data.frame()

    # Generate gene labels
    if (gene_label_show && nrow(valid_genes) > 0) {
      label_measure_units <- ggchord_device_units_per_inch(
        unlist(lapply(seq_arcs, `[[`, "x"), use.names = FALSE),
        unlist(lapply(seq_arcs, `[[`, "y"), use.names = FALSE)
      )
      gene_labels <- do.call(rbind, lapply(seq_len(nrow(valid_genes)), function(i) {
        gene <- valid_genes[i, ]
        sid <- gene$accver
        strand <- gene$strand
        seq_len <- lens[sid]
        ref <- seq_refs[[sid]]
        orient <- orientation[sid]

        sp <- min(gene$start, gene$end)
        ep <- max(gene$start, gene$end)
        frac_mid <- if (isTRUE(circular) && gene$start > gene$end) {
          ((gene$start + ((seq_len - gene$start) + gene$end) / 2) %% seq_len) /
            seq_len
        } else {
          (sp + ep) / (2 * seq_len)
        }

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
        width <- width * as.numeric(gene$.feature_width_factor %||% 1)

        r0 <- gene_track_radius(gene, sid, strand, ref$angles[idx])

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
        resolved_label_orientation <- gene_label_orientation
        feature_label_inside <- NA
        if (identical(resolved_label_orientation, "feature")) {
          feature_fraction <- if (isTRUE(circular) && gene$start > gene$end) {
            ((seq_len - gene$start) + gene$end) / seq_len
          } else (ep - sp) / seq_len
          feature_arc_length <- feature_fraction * sum(sqrt(
            diff(ref$path$x)^2 + diff(ref$path$y)^2
          ))
          measured_label <- ggchord_text_boxes(data.frame(
            text = as.character(gene$anno), text_x = 0, text_y = 0,
            text_angle = 0, hjust = .5, vjust = .5,
            size = gene_label_size, family = gene_label_family,
            fontface = gene_label_fontface,
            lineheight = gene_label_lineheight
          ), units_per_inch = label_measure_units)
          shape <- as.character(gene$.feature_shape %||% "arrow")
          head_reserve <- if (shape == "arrow") {
            arrow_head_length
          } else if (shape %in% c(
              "compact_arrow", "promoter_arrow", "primer_arrow")) {
            arrow_head_length * .65
          } else 0
          body_padding <- max(.008, measured_label$h * .24)
          available_length <- max(0,
            feature_arc_length - head_reserve - 2 * body_padding)
          feature_label_inside <- measured_label$w <= available_length
          # Plasmid feature text follows the interval direction whether it is
          # inside the polygon or immediately adjacent to it. Compact labels
          # move toward the map centre instead of changing to radial text.
          resolved_label_orientation <- "tangent"
          if (!isTRUE(feature_label_inside)) {
            centre_length <- sqrt(sum(center_pt^2))
            if (is.finite(centre_length) && centre_length > 1e-8) {
              compact_shape <- shape %in% c(
                "compact_arrow", "promoter_arrow", "primer_arrow", "marker"
              )
              adjacent_gap <- measured_label$h *
                if (compact_shape) .76 else .64
              adjacent_offset <- width / 2 + adjacent_gap + .014
              inward_x <- -center_pt[1L] / centre_length
              inward_y <- -center_pt[2L] / centre_length
              anchor_x <- center_pt[1L] + inward_x * width / 2
              anchor_y <- center_pt[2L] + inward_y * width / 2
              text_x <- text_x + inward_x * adjacent_offset
              text_y <- text_y + inward_y * adjacent_offset
            }
          }
        }
        text_angle <- switch(
          resolved_label_orientation,
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

        upright <- ggchord_normalize_text_orientation(text_angle, hjust)
        text_angle <- upright$angle
        hjust <- upright$hjust
        vjust <- 0.5

        if (identical(resolved_label_orientation, "tangent")) {
          hjust <- 0.5
        } else if (identical(resolved_label_orientation, "horizontal")) {
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

        feature_label_colour <- if (
            isTRUE(feature_label_inside) &&
            "feature_label_colour" %in% names(gene)) {
          as.character(gene$feature_label_colour)
        } else "#202020"

        data.frame(
          text = gene$anno,
          text_x = text_x,
          text_y = text_y,
          text_angle = text_angle,
          hjust = hjust,
          vjust = vjust,
          size = gene_label_size,
          family = gene_label_family,
          fontface = gene_label_fontface,
          lineheight = gene_label_lineheight,
          accver = sid,
          group = i,
          source_row = gene$.source_row,
          position_name = as.character(gene$.position_name %||% "identity"),
          base_offset = as.numeric(gene$.position_base_offset %||% 0),
          lane = as.integer(gene$.feature_stack_lane %||% 0L),
          lane_offset = as.numeric(gene$.position_lane_offset %||% 0),
          normal_offset = as.numeric(gene$.normal_offset %||% 0),
          anchor_x = anchor_x,
          anchor_y = anchor_y,
          side_flipped = side_flipped,
          feature_label_colour = feature_label_colour,
          feature_label_orientation = resolved_label_orientation,
          feature_label_inside = feature_label_inside,
          .feature_width = width,
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
        accver = id,
        stringsAsFactors = FALSE
      )
    }))
  }

})
