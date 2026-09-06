# Region, gene and fixed sequence-label geometry.
ggchord_layout_annotation_step <- function(context) evalq({
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
        sid <- gene$accver
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
          accver = sid,
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
        accver = id,
        stringsAsFactors = FALSE
      )
    }))
  }

}, envir = context)
