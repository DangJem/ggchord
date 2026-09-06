# ====================================================================
# ggplot_build.ggchord: compute the layout, inject data into cloned
# layers, add scales, and set the coordinate system.  Everything is
# driven by the plot object itself, so the plot can be printed, saved
# with ggsave(), or built with ggplot_build() any number of times and
# in any order, without cross-talk between plots.
# ====================================================================

compute_chord_geometry_single <- function(plot, geometry_cache = NULL) {
  # Step 1: collect data and parameters from the plot object
  chord <- plot$ggchord
  if (is.null(chord)) {
    ggchord_stop("Not a valid ggchord object: no data stored on the plot. ",
         "Please build the plot with ggchord().")
  }
  data_list <- chord$data
  global    <- chord$global

  seq_params    <- list()
  ribbon_params <- list()
  gene_params   <- list()
  seq_data_override <- NULL
  ribbon_data_override <- NULL
  gene_data_override <- NULL
  region_data_override <- NULL
  gene_label_params <- list()
  gene_repel_params <- list()
  axis_params   <- list()
  seq_label_params <- list()
  seq_region_params <- list()
  seq_layer_requested <- FALSE
  seq_ring_mapped <- FALSE
  ribbon_layer_requested <- FALSE
  gene_geometry_layer_requested <- FALSE
  gene_label_layer <- FALSE
  gene_repel_layer <- FALSE
  feature_stack_position <- NULL

  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    pp <- lyr$ggchord_params
    if (is.null(pp)) next
    switch(pp$type,
      seq               = {
        seq_params <- pp
        seq_layer_requested <- TRUE
        seq_ring_mapped <- "seq_ring" %in% names(lyr$ggchord_input_mapping)
        seq_data_override <- ggchord_resolve_layer_input(
          lyr, data_list$seq_data
        )
      },
      ribbon            = {
        ribbon_params <- pp
        ribbon_layer_requested <- TRUE
        ribbon_data_override <- ggchord_resolve_layer_input(
          lyr, data_list$ribbon_data
        )
      },
      gene              = {
        if (isTRUE(lyr$position$ggchord_feature_stack)) {
          pp$feature_stack_position <- lyr$position
          feature_stack_position <- lyr$position
        }
        gene_params <- pp
        gene_geometry_layer_requested <- TRUE
        gene_data_override <- pp$gene_data_override %||%
          ggchord_resolve_layer_input(lyr, data_list$gene_data)

      },
      gene_label        = {
        if (isTRUE(lyr$position$ggchord_feature_stack)) {
          feature_stack_position <- lyr$position
        }
        gene_label_params <- pp
        gene_label_layer <- TRUE
        label_data <- ggchord_resolve_layer_input(lyr, data_list$gene_data)
        if (!is.null(lyr$ggchord_input_data) ||
            length(intersect(names(lyr$ggchord_input_mapping),
                             lyr$ggchord_role_aes)) > 0) {
          gene_data_override <- label_data
        }
      },
      gene_label_repel  = {
        if (isTRUE(lyr$position$ggchord_feature_stack)) {
          feature_stack_position <- lyr$position
        }
        gene_repel_params <- pp
        gene_repel_layer <- TRUE
        label_data <- ggchord_resolve_layer_input(lyr, data_list$gene_data)
        if (!is.null(lyr$ggchord_input_data) ||
            length(intersect(names(lyr$ggchord_input_mapping),
                             lyr$ggchord_role_aes)) > 0) {
          gene_data_override <- label_data
        }
      },
      seq_label         = seq_label_params <- pp,
      seq_region        = {
        seq_region_params <- pp
        region_data_override <- ggchord_resolve_layer_input(lyr, pp$regions)
      }

    )
  }

  # --- Process sequences ---
  seq_data <- seq_data_override %||% data_list$seq_data
  seqs     <- seq_data$accver
  lens     <- setNames(seq_data$length, seqs)

  if (!is.null(seq_params$seq_order)) {
    if (!all(seq_params$seq_order %in% seqs)) {
      ggchord_stop("seq_order contains unknown sequence IDs")
    }
    seqs <- seq_params$seq_order
    lens <- lens[seqs]
  }
  n <- length(seqs)

  seq_labels    <- process_sequence_param(seq_params$seq_labels, seqs,
                                          "seq_labels", default_value = seqs)
  seqRadius     <- process_sequence_param(seq_params$seq_radius, seqs,
                                          "seq_radius", 1.0)
  orientation   <- process_sequence_param(seq_params$seq_orientation, seqs,
                                          "seq_orientation", 1)
  seq_gap       <- process_sequence_param(seq_params$seq_gap, seqs,
                                          "seq_gap", 0.03)
  seq_curvature <- process_sequence_param(seq_params$seq_curvature, seqs,
                                          "seq_curvature", 1.0)

  if (isTRUE(plot$coordinates$ggchord_genome)) {
    if (n != 1L) {
      ggchord_stop("coord_genome() requires exactly one sequence in seq_data")
    }
    incompatible <- vapply(plot$layers, function(layer) {
      (layer$ggchord_params$type %||% "") %in% c("ribbon", "link")
    }, logical(1))
    if (any(incompatible)) {
      ggchord_stop("coord_genome() does not support ribbon or link layers")
    }
    if (!is.null(seq_params$seq_gap)) {
      ggchord_stop("coord_genome() owns the opening; use its `gap` argument instead of geom_seq(seq_gap = ...)")
    }
    if (!is.null(seq_params$seq_orientation)) {
      ggchord_stop("coord_genome() owns genomic direction; use its `direction` argument instead of geom_seq(seq_orientation = ...)")
    }
    seq_gap[] <- plot$coordinates$genome_gap / 360
    orientation[] <- if (identical(plot$coordinates$genome_direction, "clockwise")) -1 else 1
  }

  # Rings are explicit input roles. Their scale values are radii, so no ring
  # count or spacing is guessed from the data. This deliberately keeps the
  # existing seq_radius interface unchanged for plots without a ring mapping.
  seq_ring <- NULL
  if (isTRUE(seq_ring_mapped)) {
    if (!is.null(seq_params$seq_radius)) {
      ggchord_stop(
        "geom_seq(): `seq_radius` cannot be combined with a `seq_ring` ",
        "mapping; set radii with scale_seq_ring_manual(values = ...)"
      )
    }
    ring_scale <- plot$scales$get_scales("seq_ring")
    if (is.null(ring_scale)) {
      ggchord_stop(
        "A mapped `seq_ring` requires scale_seq_ring_manual(values = ...)"
      )
    }
    raw_ring <- as.character(seq_data$seq_ring[match(seqs, seq_data$accver)])
    if (anyNA(raw_ring) || any(!nzchar(raw_ring))) {
      ggchord_stop("Mapped `seq_ring` values must be non-missing")
    }
    trained_ring_scale <- ring_scale$clone()
    trained_ring_scale$train(raw_ring)
    mapped_radius <- suppressWarnings(as.numeric(
      trained_ring_scale$map(raw_ring)
    ))
    if (anyNA(mapped_radius) || any(!is.finite(mapped_radius)) ||
        any(mapped_radius <= 0)) {
      ggchord_stop(
        "scale_seq_ring_manual(): every used ring must map to a finite ",
        "positive radius"
      )
    }
    seq_ring <- stats::setNames(raw_ring, seqs)
    seqRadius <- stats::setNames(mapped_radius, seqs)
  }

  if (!is.numeric(seqRadius) || any(!is.finite(seqRadius)) || any(seqRadius <= 0)) {
    ggchord_stop("seq_radius must contain finite positive numbers")
  }
  if (!is.numeric(orientation) || any(!is.finite(orientation)) ||
      any(!orientation %in% c(-1, 1))) {
    ggchord_stop("seq_orientation can only be 1 or -1")
  }
  if (!is.numeric(seq_gap) || any(!is.finite(seq_gap)) ||
      any(seq_gap < 0 | seq_gap >= 0.5)) {
    ggchord_stop("seq_gap must be in the [0, 0.5) range")
  }
  if (!is.numeric(seq_curvature) || any(!is.finite(seq_curvature))) {
    ggchord_stop("seq_curvature must contain finite numbers")
  }

  if (!is.null(seq_params$seq_colors)) {
    seq_colors <- process_sequence_param(seq_params$seq_colors, seqs, "seq_colors")
  } else {
    pal <- chord_default_palette(n)
    seq_colors <- setNames(pal, seqs)
  }

  # --- Process ribbons ---
  ribbon_gap_auto <- is.null(ribbon_params$ribbon_gap)
  ribbonGap  <- process_sequence_param(ribbon_params$ribbon_gap %||% 0.15,
                                       seqs, "ribbon_gap", 0.15)
  ribbon_color_scheme <- ribbon_params$ribbon_color_scheme %||% "pident"
  ribbon_color_by     <- ribbon_params$ribbon_color_by
  ribbon_color_limits <- ribbon_params$ribbon_color_limits
  ribbon_color_breaks <- ribbon_params$ribbon_color_breaks
  ribbon_color_name   <- ribbon_params$ribbon_color_name
  ribbon_alpha    <- ribbon_params$ribbon_alpha %||% 0.42
  ribbon_alpha_by <- ribbon_params$ribbon_alpha_by
  ribbon_alpha_range <- ribbon_params$ribbon_alpha_range %||% c(0.15, 0.9)
  ribbon_ctrl_pt  <- ribbon_params$ribbon_ctrl_point %||% c(0, 0)
  ribbon_outline_by   <- ribbon_params$ribbon_outline_by
  ribbon_outline_colors <- ribbon_params$ribbon_outline_colors
  ribbon_linetype_by  <- ribbon_params$ribbon_linetype_by
  ribbon_linetypes    <- ribbon_params$ribbon_linetypes
  ribbon_direction    <- ribbon_params$ribbon_direction %||% "none"
  ribbon_direction_colors <- ribbon_params$ribbon_direction_colors %||% c(same = "black", reverse = "grey50")
  ribbon_direction_linetypes <- ribbon_params$ribbon_direction_linetypes %||% c(same = "solid", reverse = "dashed")
  ribbon_direction_alpha <- ribbon_params$ribbon_direction_alpha %||% c(same = 1, reverse = 0.45)

  ribbon_colors <- ribbon_params$ribbon_colors
  if (!ribbon_color_scheme %in% c("pident", "query", "subject", "single")) {
    ggchord_stop("ribbon_color_scheme must be 'pident', 'query', 'subject', or 'single'")
  }
  if (!is.null(ribbon_color_by)) {
    ribbon_color_scheme <- "value"
    ribbon_color_name <- ribbon_color_name %||% ribbon_color_by
  }
  if (!is.numeric(ribbon_alpha) || length(ribbon_alpha) != 1 ||
      !is.finite(ribbon_alpha) || ribbon_alpha < 0 || ribbon_alpha > 1) {
    ggchord_stop("ribbon_alpha must be in the [0, 1] range")
  }
  if (!is.numeric(ribbonGap) || any(!is.finite(ribbonGap))) {
    ggchord_stop("ribbon_gap must contain finite numbers")
  }
  if (!is.null(ribbon_color_limits) &&
      (!is.numeric(ribbon_color_limits) || length(ribbon_color_limits) != 2 ||
       !is.finite(ribbon_color_limits[1]) || !is.finite(ribbon_color_limits[2]) ||
       ribbon_color_limits[1] >= ribbon_color_limits[2])) {
    ggchord_stop("ribbon_color_limits must be a length-2 increasing numeric vector")
  }
  if (!is.null(ribbon_color_breaks) &&
      (!is.numeric(ribbon_color_breaks) || any(!is.finite(ribbon_color_breaks)))) {
    ggchord_stop("ribbon_color_breaks must be a finite numeric vector")
  }
  if (!is.numeric(ribbon_alpha_range) || length(ribbon_alpha_range) != 2 ||
      any(!is.finite(ribbon_alpha_range)) || ribbon_alpha_range[1] < 0 ||
      ribbon_alpha_range[2] > 1 || ribbon_alpha_range[1] > ribbon_alpha_range[2]) {
    ggchord_stop("ribbon_alpha_range must be two increasing values within [0, 1]")
  }
  if (!ribbon_direction %in% c("none", "alpha", "outline", "linetype")) {
    ggchord_stop("ribbon_direction must be 'none', 'alpha', 'outline', or 'linetype'")
  }

  # ribbon_colors validation only runs when ribbon_data is actually present
  ribbon_data <- if (ribbon_layer_requested) {
    ribbon_data_override %||%
      data_list$ribbon_data
  } else {
    NULL
  }
  ribbon_stat_report <- NULL
  ribbon_stat_data <- NULL
  ribbon_stat <- ribbon_params$ribbon_stat
  if (!is.null(ribbon_stat) && !is.null(ribbon_data)) {
    if (identical(ribbon_stat$type, "bundle")) {
      computed <- bundle_ggchord_ribbons(
        ribbon_data = ribbon_data,
        seq_data = seq_data,
        bins = ribbon_stat$bins,
        min_bundle = ribbon_stat$min_bundle,
        weight = ribbon_stat$weight,
        group_by = ribbon_stat$group_by
      )
      computed$data$bundle_n <- computed$data$.bundle_n
      computed$data$bundle_weight <- computed$data$.bundle_weight
      computed$data$density <- computed$data$.bundle_density
    } else if (identical(ribbon_stat$type, "density")) {
      computed <- ggchord_ribbon_density(
        ribbon_data = ribbon_data,
        seq_data = seq_data,
        bins = ribbon_stat$bins,
        weight = ribbon_stat$weight,
        group_by = ribbon_stat$group_by,
        caller = "stat_ribbon_density()"
      )
    } else {
      ggchord_stop("Unknown ggchord ribbon stat: ", ribbon_stat$type)
    }
    ribbon_data <- computed$data
    ribbon_stat_data <- computed$data
    ribbon_stat_report <- computed$report
  }
  has_ribbon_data <- !is.null(ribbon_data) && nrow(ribbon_data) > 0

  if (has_ribbon_data) {
    if (is.null(ribbon_colors)) {
      ribbon_colors <- switch(ribbon_color_scheme,
        single = "steelblue",
        query  = {
          mix <- 0.5
          sapply(seq_colors, function(col) {
            cols <- grDevices::col2rgb(col)
            light_cols <- cols + (255 - cols) * mix
            grDevices::rgb(light_cols[1,], light_cols[2,], light_cols[3,],
                maxColorValue = 255)
          })
        },
        subject = {
          mix <- 0.5
          sapply(seq_colors, function(col) {
            cols <- grDevices::col2rgb(col)
            light_cols <- cols + (255 - cols) * mix
            grDevices::rgb(light_cols[1,], light_cols[2,], light_cols[3,],
                maxColorValue = 255)
          })
        },
        pident = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF",
                    "#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF",
                    "#B4DE2CFF","#FDE725FF"),
        value = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF",
                    "#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF",
                    "#B4DE2CFF","#FDE725FF"))
    }
    if (ribbon_color_scheme %in% c("query", "subject")) {
      ribbon_colors <- process_sequence_param(ribbon_colors, seqs,
                                              "ribbon_colors")
    } else if (ribbon_color_scheme %in% c("pident", "value") && length(ribbon_colors) < 2) {
      ggchord_stop("The 'pident'/'value' scheme requires at least two ribbon_colors")
    } else if (ribbon_color_scheme == "single" && length(ribbon_colors) < 1) {
      ggchord_stop("The 'single' scheme requires at least one ribbon_colors")
    }

    if (!is.null(ribbon_color_by)) {
      if (!ribbon_color_by %in% colnames(ribbon_data)) {
        ggchord_stop("ribbon_color_by column '", ribbon_color_by, "' not found in ribbon_data")
      }
      if (!is.numeric(ribbon_data[[ribbon_color_by]]) ||
          any(!is.finite(ribbon_data[[ribbon_color_by]]))) {
        ggchord_stop("ribbon_color_by column '", ribbon_color_by, "' must be numeric and finite")
      }
    }
    if (!is.null(ribbon_alpha_by)) {
      if (!ribbon_alpha_by %in% colnames(ribbon_data)) {
        ggchord_stop("ribbon_alpha_by column '", ribbon_alpha_by, "' not found in ribbon_data")
      }
      if (!is.numeric(ribbon_data[[ribbon_alpha_by]]) ||
          any(!is.finite(ribbon_data[[ribbon_alpha_by]]))) {
        ggchord_stop("ribbon_alpha_by column '", ribbon_alpha_by, "' must be numeric and finite")
      }
    }
    if (!is.null(ribbon_outline_by) &&
        !ribbon_outline_by %in% colnames(ribbon_data)) {
      ggchord_stop("ribbon_outline_by column '", ribbon_outline_by, "' not found in ribbon_data")
    }
    if (!is.null(ribbon_linetype_by) &&
        !ribbon_linetype_by %in% colnames(ribbon_data)) {
      ggchord_stop("ribbon_linetype_by column '", ribbon_linetype_by, "' not found in ribbon_data")
    }
  }

  # --- Process genes ---
  gene_off  <- gene_params$gene_offset %||% 0.1
  gene_w    <- gene_params$gene_width %||% 0.05
  gene_cs   <- gene_params$gene_color_scheme %||% "strand"
  gene_cols <- gene_params$gene_colors
  gene_ord  <- gene_params$gene_order
  # Gene label settings come from the dedicated geom_gene_label() layer, with
  # the legacy geom_gene() arguments as fallback.
  # The repel layer takes priority over the fixed label layer
  lbl <- if (gene_repel_layer) gene_repel_params else gene_label_params
  gene_ls   <- gene_label_layer || gene_repel_layer ||
    isTRUE(gene_params$show_label_override) ||
    isTRUE(gene_params$gene_label_show)
  gene_lsz  <- lbl$gene_label_size %||%
    gene_params$label_size_override %||%
    gene_params$gene_label_size %||%
    ggchord_theme_text_size(plot, "ggchord.gene.label", 2.5)
  # Repelled labels now use mode-owned deterministic positioning. Manual
  # rotation and offsets remain available through geom_gene_label(), but are
  # intentionally not inherited by geom_gene_label_repel().
  if (gene_repel_layer) {
    gene_lr <- gene_lro <- gene_lco <- 0
    gene_lcl <- TRUE
  } else {
    gene_lr <- lbl$gene_label_rotation %||%
      gene_params$gene_label_rotation %||% 0
    gene_lro <- lbl$gene_label_radial_offset %||%
      gene_params$gene_label_radial_offset %||% 0
    gene_lco <- lbl$gene_label_circum_offset %||%
      gene_params$gene_label_circum_offset %||% 0
    gene_lcl <- if (is.null(lbl$gene_label_circum_limit)) {
      if (is.null(gene_params$gene_label_circum_limit)) TRUE
      else gene_params$gene_label_circum_limit
    } else lbl$gene_label_circum_limit
  }
  gene_lwrap  <- lbl$gene_label_wrap %||% gene_params$gene_label_wrap
  gene_lorientation <- if (gene_repel_layer) {
    "radial"
  } else {
    lbl$gene_label_orientation %||% "horizontal"
  }
  gene_loverlap <- if (gene_repel_layer) {
    "allow"
  } else {
    lbl$gene_label_overlap %||% "hide"
  }
  gene_lrepel_layer <- gene_repel_layer
  gene_lrepel_maxov <- gene_repel_params$max_overlaps %||% Inf
  gene_lrepel_layout <- gene_repel_params$gene_label_layout %||% "radial"
  gene_lrepel_fit    <- gene_repel_params$gene_label_fit %||% "wrap"
  gene_lrepel_lines  <- gene_repel_params$gene_label_max_lines %||% 2L
  gene_lrepel_side   <- if (gene_repel_layer) {
    gene_repel_params$gene_label_side %||% "outside"
  } else {
    lbl$gene_label_side %||% "outside"
  }
  gene_lrepel_segment_overlap <-
    gene_repel_params$gene_label_segment_overlap %||% "fade"
  gene_lrepel_segment_overlap_alpha <-
    gene_repel_params$gene_label_segment_overlap_alpha %||% 0.18
  gene_lrepel_ltype  <- gene_repel_params$gene_label_segment_linetype %||% "auto"

  if (!gene_cs %in% c("strand", "manual")) {
    ggchord_stop("gene_color_scheme must be 'strand' or 'manual'")
  }
  if (!is.numeric(gene_lsz) || length(gene_lsz) != 1 || !is.finite(gene_lsz) ||
      gene_lsz <= 0) {
    ggchord_stop("gene_label_size must be a finite positive number")
  }
  if (!is.null(gene_lwrap) && (!is.numeric(gene_lwrap) ||
      length(gene_lwrap) != 1 || !is.finite(gene_lwrap) || gene_lwrap < 0)) {
    ggchord_stop("gene_label_wrap must be NULL or a finite non-negative number")
  }
  if (!is.numeric(gene_lrepel_maxov) || length(gene_lrepel_maxov) != 1 ||
      is.na(gene_lrepel_maxov) || gene_lrepel_maxov < 0) {
    ggchord_stop("max_overlaps must be a non-negative number or Inf")
  }

  geneGap    <- process_gene_param(gene_off, seqs, "gene_offset", 0.1, FALSE)
  geneWidth  <- process_gene_param(gene_w, seqs, "gene_width", 0.05, FALSE)

  # Feature geometry is resolved before coordinate generation because these
  # values change the actual polygon, not only its appearance. A user-supplied
  # feature-shape scale is cloned and trained here so the layout and legend use
  # exactly the same category-to-shape mapping.
  gene_data_layout <- if (gene_geometry_layer_requested || gene_label_layer ||
      gene_repel_layer) {
    gene_data_override %||% data_list$gene_data
  } else {
    NULL
  }
  ribbon_obstacles <- plot$ggchord$obstacles %||% ggchord_empty_obstacles()
  if (!is.null(gene_data_layout) && !"anno" %in% names(gene_data_layout)) {
    gene_data_layout$anno <- rep(NA_character_, nrow(gene_data_layout))
  }
  feature_shape_pal <- NULL
  feature_shape_order <- NULL
  if (isTRUE(gene_params$is_feature) && !is.null(gene_data_layout) &&
      nrow(gene_data_layout) > 0L) {
    raw_shape <- as.character(
      gene_data_layout$.feature_shape_raw %||%
        rep(gene_params$feature_shape %||% "arrow", nrow(gene_data_layout))
    )
    feature_shape_order <- unique(raw_shape)
    if (isTRUE(gene_params$feature_shape_mapped)) {
      shape_scale <- plot$scales$get_scales("feature_shape")
      if (!is.null(shape_scale)) {
        shape_scale <- shape_scale$clone()
        shape_scale$train(raw_shape)
        mapped_shape <- as.character(shape_scale$map(raw_shape))
      } else {
        allowed_shape <- c("arrow", "block", "chevron", "lollipop")
        if (all(feature_shape_order %in% allowed_shape)) {
          feature_shape_pal <- stats::setNames(
            feature_shape_order, feature_shape_order
          )
        } else {
          feature_shape_pal <- stats::setNames(
            rep(allowed_shape, length.out = length(feature_shape_order)),
            feature_shape_order
          )
        }
        mapped_shape <- unname(feature_shape_pal[raw_shape])
      }
    } else {
      mapped_shape <- raw_shape
    }
    allowed_shape <- c("arrow", "block", "chevron", "lollipop")
    mapped_shape[is.na(mapped_shape)] <- "arrow"
    if (any(!mapped_shape %in% allowed_shape)) {
      ggchord_stop(
        "feature_shape scale values must use 'arrow', 'block', ",
        "'chevron', or 'lollipop'"
      )
    }
    if (is.null(feature_shape_pal) &&
        isTRUE(gene_params$feature_shape_mapped)) {
      feature_shape_pal <- stats::setNames(
        mapped_shape[match(feature_shape_order, raw_shape)],
        feature_shape_order
      )
    }
    gene_data_layout <- as.data.frame(gene_data_layout, stringsAsFactors = FALSE)
    gene_data_layout$.feature_shape <- mapped_shape
  }
  if (!is.null(feature_stack_position) && !is.null(gene_data_layout)) {
    gene_data_layout <- ggchord_stack_feature_tracks(
      gene_data_layout, feature_stack_position
    )
  }
  geneLabelRadialOffset <- process_gene_param(gene_lro, seqs,
                                              "gene_label_radial_offset", 0, FALSE)
  geneLabelCircumOffset <- process_gene_param(gene_lco, seqs,
                                              "gene_label_circum_offset", 0, FALSE)
  geneLabelCircumLimit  <- process_gene_param(gene_lcl, seqs,
                                              "gene_label_circum_limit", TRUE, TRUE)
  geneLabelRotation     <- process_gene_param(gene_lr, seqs,
                                              "gene_label_rotation", 0, FALSE)

  # --- Process axes ---
  axis_theme <- ggchord_plot_settings(plot)$axis
  show_axis <- !isTRUE(axis_theme$hidden)
  # Pass physical distances as inches. compute_chord_layout() converts them
  # against the actual curved sequence span once that geometry is available.
  axis_unit_data <- ggchord_unit_inches
  axisGap    <- process_sequence_param(axis_unit_data(axis_theme$gap),
                                       seqs, "axis.gap", 0.04)
  axisMaj    <- process_sequence_param(axis_params$axis_tick_major_number %||% 3,
                                       seqs, "axis_tick_major_number", 3)
  axisMajLen <- process_sequence_param(axis_unit_data(axis_theme$ticks.length),
                                       seqs, "axis.ticks.length", 0.02)
  axisMin    <- process_sequence_param(axis_params$axis_tick_minor_number %||% 4,
                                       seqs, "axis_tick_minor_number", 4)
  axisMinLen <- process_sequence_param(
    axis_unit_data(axis_theme$minor.ticks.length),
    seqs, "axis.minor.ticks.length", 0.01
  )
  axis_theme_size <- ggchord_theme_text_size(
    plot, "ggchord.axis.text", 3
  )
  labelSize  <- process_sequence_param(
    axis_theme_size,
    seqs, "axis_label_size", axis_theme_size
  )
  labelOffset <- process_sequence_param(axis_unit_data(axis_theme$text.offset),
                                        seqs, "axis.text.offset", 0.02)
  axisLabelHide <- isTRUE(axis_theme$text.check.overlap)
  axisLabelOri <- process_axis_orientation(
    axis_theme$text.orientation, seqs
  )
  if (!is.logical(show_axis) || length(show_axis) != 1 || is.na(show_axis)) {
    ggchord_stop("show_axis must be TRUE or FALSE")
  }
  axis_numeric <- list(
    axis_gap = axisGap,
    axis_tick_major_length = axisMajLen,
    axis_tick_minor_length = axisMinLen,
    axis_label_size = labelSize,
    axis_text_offset = labelOffset
  )
  for (nm in names(axis_numeric)) {
    value <- axis_numeric[[nm]]
    if (!is.numeric(value) || any(!is.finite(value))) {
      ggchord_stop(nm, " must contain finite numbers")
    }
  }
  if (!is.numeric(axisMaj) || any(!is.finite(axisMaj)) ||
      any(axisMaj < 1 | axisMaj != as.integer(axisMaj))) {
    ggchord_stop("axis_tick_major_number must contain positive integers")
  }
  if (!is.numeric(axisMin) || any(!is.finite(axisMin)) ||
      any(axisMin < 0 | axisMin != as.integer(axisMin))) {
    ggchord_stop("axis_tick_minor_number must contain non-negative integers")
  }
  axis_breaks <- axis_minor_breaks <- axis_labels <- NULL
  position_scale <- plot$scales$get_scales("seq_position")
  if (!is.null(position_scale)) {
    axis_breaks <- axis_minor_breaks <- axis_labels <- setNames(
      vector("list", length(seqs)), seqs
    )
    for (id in seqs) {
      sc <- position_scale$clone()
      sc$train(c(0, lens[[id]]))
      br <- sc$get_breaks()
      br <- br[is.finite(br) & br >= 0 & br <= lens[[id]]]
      axis_breaks[[id]] <- br
      minor <- sc$get_breaks_minor()
      axis_minor_breaks[[id]] <- minor[
        is.finite(minor) & minor >= 0 & minor <= lens[[id]]
      ]
      axis_labels[[id]] <- sc$get_labels(br)
    }
  }

  # --- Process sequence labels ---
  seq_label_text <- NULL
  seq_label_radius <- NULL
  seq_label_rotation <- NULL
  seq_label_size <- NULL
  seq_label_orientation <- "arc"
  seq_label_hjust <- NULL
  seq_label_vjust <- NULL
  if (length(seq_label_params) > 0) {
    seq_label_text <- if (is.null(seq_label_params$seq_labels)) {
      seq_labels
    } else {
      # Process through the standard parameter helper so that unnamed vectors
      # are matched positionally to the sequences (named by accver).
      process_sequence_param(seq_label_params$seq_labels, seqs, "seq_labels",
                             default_value = seqs)
    }
    seq_label_radius <- process_sequence_param(
      seq_label_params$seq_label_radius, seqs, "seq_label_radius", 1)
    seq_label_rotation <- process_sequence_param(
      seq_label_params$seq_label_rotation, seqs, "seq_label_rotation", 0)
    seq_theme_size <- ggchord_theme_text_size(
      plot, "ggchord.seq.label", 3
    )
    seq_label_size <- process_sequence_param(
      seq_label_params$seq_label_size, seqs, "seq_label_size", seq_theme_size
    )
    seq_label_orientation <- seq_label_params$seq_label_orientation %||% "arc"
    seq_label_hjust <- if (is.null(seq_label_params$seq_label_hjust)) {
      if (identical(seq_label_orientation, "arc")) {
        process_sequence_param(-0.2, seqs, "seq_label_hjust", -0.2)
      } else {
        NULL
      }
    } else {
      process_sequence_param(seq_label_params$seq_label_hjust, seqs,
                             "seq_label_hjust", 0.5)
    }
    seq_label_vjust <- if (is.null(seq_label_params$seq_label_vjust)) {
      NULL
    } else {
      process_sequence_param(seq_label_params$seq_label_vjust, seqs,
                             "seq_label_vjust", 0.5)
    }
    seq_label_numeric <- list(
      seq_label_radius = seq_label_radius,
      seq_label_rotation = seq_label_rotation,
      seq_label_size = seq_label_size,
      seq_label_hjust = seq_label_hjust,
      seq_label_vjust = seq_label_vjust
    )
    for (nm in names(seq_label_numeric)) {
      value <- seq_label_numeric[[nm]]
      if (!is.null(value) && (!is.numeric(value) || any(!is.finite(value)))) {
        ggchord_stop(nm, " must contain finite numbers")
      }
    }
    if (any(seq_label_size <= 0)) {
      ggchord_stop("seq_label_size must contain positive numbers")
    }
  }

  # --- Process sequence-region highlight data ---
  region_data <- region_data_override %||% seq_region_params$regions
  region_fill   <- seq_region_params$region_fill %||% "#F59E0B"
  region_color  <- seq_region_params$region_color %||% "#B45309"
  region_alpha  <- seq_region_params$region_alpha %||% 0.25
  region_width  <- seq_region_params$region_width %||% 0.08
  region_offset <- seq_region_params$region_offset %||% 0
  region_side   <- seq_region_params$region_side %||% "inside"
  if (!is.null(region_data) && !is.data.frame(region_data)) {
    ggchord_stop("geom_seq_region(): regions must be a data.frame")
  }

  # ====================================================================
  # Step 2: compute the layout
  # ====================================================================
  coord_rotation <- if (isTRUE(plot$coordinates$ggchord_coord)) {
    plot$coordinates$rotation + if (isTRUE(plot$coordinates$ggchord_genome) &&
        identical(plot$coordinates$genome_direction, "clockwise")) {
      plot$coordinates$genome_gap
    } else 0
  } else {
    global$rotation
  }

  layout <- compute_chord_layout(
    seqs = seqs, lens = lens, seq_labels = seq_labels,
    seq_colors = seq_colors, seqRadius = seqRadius,
    seq_curvature = seq_curvature, orientation = orientation,
    seq_gap = seq_gap,
    ribbon_data = ribbon_data, ribbonGap = ribbonGap,
    ribbon_gap_auto = ribbon_gap_auto,
    link_avoid = ribbon_params$link_avoid %||% "none",
    ribbon_obstacles = ribbon_obstacles,
    ribbon_color_scheme = ribbon_color_scheme,
    ribbon_colors = ribbon_colors, ribbon_alpha = ribbon_alpha,
    ribbon_color_by = ribbon_color_by,
    ribbon_color_limits = ribbon_color_limits,
    ribbon_color_breaks = ribbon_color_breaks,
    ribbon_color_name = ribbon_color_name,
    ribbon_alpha_by = ribbon_alpha_by,
    ribbon_alpha_range = ribbon_alpha_range,
    ribbon_outline_by = ribbon_outline_by,
    ribbon_outline_colors = ribbon_outline_colors,
    ribbon_linetype_by = ribbon_linetype_by,
    ribbon_linetypes = ribbon_linetypes,
    ribbon_direction = ribbon_direction,
    ribbon_direction_colors = ribbon_direction_colors,
    ribbon_direction_linetypes = ribbon_direction_linetypes,
    ribbon_direction_alpha = ribbon_direction_alpha,
    ribbon_ctrl_point = ribbon_ctrl_pt,
    region_data = region_data,
    region_fill = region_fill,
    region_color = region_color,
    region_alpha = region_alpha,
    region_width = region_width,
    region_offset = region_offset,
    region_side = region_side,
    gene_data = gene_data_layout,
    draw_gene_geometry = gene_geometry_layer_requested,
    geneGap = geneGap, geneWidth = geneWidth,
    geneLabelRadialOffset = geneLabelRadialOffset,
    geneLabelCircumOffset = geneLabelCircumOffset,
    geneLabelCircumLimit = geneLabelCircumLimit,
    geneLabelRotation = geneLabelRotation,
    gene_label_show = gene_ls, gene_label_size = gene_lsz,
    gene_label_wrap = gene_lwrap,
    gene_label_fit = gene_lrepel_fit,
    gene_label_max_lines = gene_lrepel_lines,
    gene_label_orientation = gene_lorientation,
    gene_label_overlap = gene_loverlap,
    gene_label_repel_layer = gene_lrepel_layer,
    gene_label_repel_max_overlaps = gene_lrepel_maxov,
    gene_label_layout = gene_lrepel_layout,
    gene_label_side = gene_lrepel_side,
    gene_label_segment_overlap = gene_lrepel_segment_overlap,
    gene_label_segment_overlap_alpha = gene_lrepel_segment_overlap_alpha,
    gene_label_segment_linetype = gene_lrepel_ltype,
    gene_color_scheme = gene_cs, gene_colors = gene_cols,
    gene_order = gene_ord,
    feature_shape_pal = feature_shape_pal,
    feature_shape_order = feature_shape_order,
    seq_label_text = seq_label_text,
    seq_label_radius = seq_label_radius,
    seq_label_rotation = seq_label_rotation,
    seq_label_size = seq_label_size,
    seq_label_orientation = seq_label_orientation,
    seq_label_hjust = seq_label_hjust,
    seq_label_vjust = seq_label_vjust,
    axisGap = axisGap, axisMaj = axisMaj, axisMajLen = axisMajLen,
    axisMin = axisMin, axisMinLen = axisMinLen,
    labelSize = labelSize, labelOffset = labelOffset,
    axisLabelOrientation = axisLabelOri,
    axis_breaks = axis_breaks,
    axis_minor_breaks = axis_minor_breaks,
    axis_labels = axis_labels,
    axis_label_hide_overlaps = axisLabelHide,
    show_axis = show_axis,
    rotation = coord_rotation, debug = global$debug,
    geometry_cache = geometry_cache
  )

  layout$seq_ring <- seq_ring
  layout$seq_ring_radius <- if (is.null(seq_ring)) NULL else seqRadius
  layout$ribbon_stat_data <- ribbon_stat_data
  layout$ribbon_stat_report <- ribbon_stat_report
  if (!is.null(seq_ring) && length(layout$seq_arcs)) {
    layout$seq_arcs <- lapply(names(layout$seq_arcs), function(id) {
      arc <- layout$seq_arcs[[id]]
      arc$seq_ring <- unname(seq_ring[[id]])
      arc
    }) |>
      stats::setNames(names(layout$seq_arcs))
  }

  layout
}
#' Extract one drawable component from a computed layout
#' @noRd
ggchord_layout_component <- function(layout, type, fallback = data.frame()) {
  switch(type,
    seq = if (length(layout$seq_arcs) > 0) do.call(rbind, layout$seq_arcs) else fallback,
    ribbon = layout$ribbon_polys %||% fallback,
    link = layout$link_lines %||% fallback,
    gene_poly = layout$gene_polys %||% fallback,
    gene_text = layout$gene_labels %||% fallback,
    gene_text_repel = layout$gene_labels %||% fallback,
    gene_label_segment = layout$gene_label_segments %||% fallback,
    gene_label_repel = ggchord_repel_geometry(layout),
    seq_label = layout$seq_labels_df %||% fallback,
    seq_region = layout$region_polys %||% fallback,
    restriction_site = layout$restriction_sites %||% fallback,
    axis_line = layout$axis_lines %||% fallback,
    axis_seg = layout$axis_ticks %||% fallback,
    axis_text = {
      d <- layout$axis_ticks %||% fallback
      if (nrow(d) > 0 && "label" %in% names(d)) d[!is.na(d$label), , drop = FALSE] else d
    },
    axis = ggchord_axis_geometry(layout),
    fallback
  )
}

#' Compute a plot-owned, per-layer geometry registry
#' @noRd
compute_chord_geometry <- function(plot) {
  chord <- plot$ggchord
  if (is.null(chord)) {
    ggchord_stop("Not a valid ggchord object: no data stored on the plot")
  }
  # Layers added through ordinary ggplot2 mechanisms may not have passed the
  # list branch of +.ggchord. Assign deterministic IDs before grouping.
  next_id <- 1L
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    if (is.null(lyr$ggchord_type)) next
    if (is.null(lyr$ggchord_layer_id)) {
      lyr$ggchord_layer_id <- sprintf("layer-%04d", next_id)
      plot$layers[[i]] <- lyr
    }
    next_id <- next_id + 1L
  }

  geometry_cache <- new.env(parent = emptyenv())
  avoid_requested <- any(vapply(plot$layers, function(layer) {
    params <- layer$ggchord_params
    !is.null(params) && params$type %in% c("ribbon", "link") &&
      !identical(params$link_avoid %||% "none", "none") &&
      is.null(params$ribbon_gap %||% params$link_gap)
  }, logical(1)))
  plot$ggchord$obstacles <- if (avoid_requested) {
    ggchord_collect_obstacles(plot, geometry_cache)
  } else {
    ggchord_empty_obstacles()
  }
  primary <- compute_chord_geometry_single(plot, geometry_cache)
  primary$obstacles <- plot$ggchord$obstacles
  ids <- vapply(plot$layers, function(x) x$ggchord_layer_id %||% "", character(1))
  groups <- split(which(nzchar(ids)), ids[nzchar(ids)])

  group_type <- vapply(groups, function(idx) {
    types <- vapply(idx, function(i) plot$layers[[i]]$ggchord_params$type %||% "",
                    character(1))
    main <- types[types %in% c(
      "seq", "ribbon", "link", "gene", "gene_label", "gene_label_repel", "axis",
      "seq_label", "seq_region", "restriction_site"
    )]
    if (length(main)) main[length(main)] else ""
  }, character(1))
  type_counts <- table(group_type[nzchar(group_type)])

  first_group <- function(type) {
    hit <- names(group_type)[group_type == type]
    if (length(hit)) groups[[hit[1]]] else integer(0)
  }
  seq_dep <- first_group("seq")
  gene_dep <- first_group("gene")
  ribbon_dep <- first_group("ribbon")
  gene_geometry_dep <- unlist(
    groups[names(group_type)[group_type == "gene"]], use.names = FALSE
  )

  registry <- list()
  inputs <- list()
  layouts <- list()
  for (id in names(groups)) {
    idx <- groups[[id]]
    main_type <- group_type[[id]]
    needs_own <- nzchar(main_type) && type_counts[[main_type]] > 1
    sub_layout <- primary
    if (isTRUE(needs_own)) {
      deps <- seq_dep
      if (main_type %in% c("gene_label", "gene_label_repel")) {
        deps <- c(deps, gene_dep)
      }
      if (main_type == "ribbon") deps <- c(deps, gene_geometry_dep)
      sub_plot <- plot
      sub_plot$layers <- plot$layers[sort(unique(c(deps, idx)))]
      sub_layout <- compute_chord_geometry_single(sub_plot, geometry_cache)
    }
    if (main_type == "link") {
      link_layer <- plot$layers[[idx[1L]]]
      link_input <- ggchord_resolve_layer_input(link_layer)
      sub_layout$link_lines <- ggchord_attach_input_columns(
        ggchord_link_geometry(link_input, link_layer$ggchord_params, primary), link_input)
    }
    if (main_type == "restriction_site") {
      site_layer <- plot$layers[[idx[1L]]]
      site_input <- ggchord_resolve_layer_input(site_layer)
      sub_layout$restriction_sites <- ggchord_restriction_geometry(
        site_input, site_layer$ggchord_params, primary, chord$data$seq_data
      )
    }
    layouts[[id]] <- sub_layout
    registry[[id]] <- list()
    inputs[[id]] <- list()
    for (i in idx) {
      lyr <- plot$layers[[i]]
      component <- lyr$ggchord_type
      registry[[id]][[component]] <- ggchord_layout_component(
        sub_layout, component, lyr$ggchord_placeholder %||% data.frame()
      )
      fallback <- switch(component,
        seq = chord$data$seq_data,
        axis_line = chord$data$seq_data,
        axis_seg = chord$data$seq_data,
        axis_text = chord$data$seq_data,
        axis = chord$data$seq_data,
        seq_label = chord$data$seq_data,
        ribbon = chord$data$ribbon_data,
        gene_poly = chord$data$gene_data,
        gene_text = chord$data$gene_data,
        gene_text_repel = chord$data$gene_data,
        gene_label_segment = chord$data$gene_data,
        gene_label_repel = chord$data$gene_data,
        seq_region = lyr$ggchord_params$regions,
        restriction_site = lyr$ggchord_input_data,
        NULL
      )
      inputs[[id]][[component]] <- if (
          component == "ribbon" && !is.null(sub_layout$ribbon_stat_data)) {
        sub_layout$ribbon_stat_data
      } else {
        ggchord_resolve_layer_input(lyr, fallback)
      }
      registry[[id]][[component]] <- ggchord_attach_input_columns(
        registry[[id]][[component]], inputs[[id]][[component]]
      )
    }
  }
  primary$layer_geometry <- registry
  primary$layer_inputs <- inputs
  primary$layer_layouts <- layouts

  # Plot limits must see all independent layers, not only the compatibility
  # fields in the primary layout.
  collect <- function(component) {
    values <- lapply(registry, `[[`, component)
    values <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, values)
    if (length(values)) ggchord_rbind_fill(values) else data.frame()
  }
  for (pair in list(
    c("ribbon_polys", "ribbon"), c("link_lines", "link"), c("gene_polys", "gene_poly"),
    c("gene_labels", "gene_text"), c("gene_label_segments", "gene_label_segment"),
    c("seq_labels_df", "seq_label"), c("region_polys", "seq_region"),
    c("axis_lines", "axis_line"), c("axis_ticks", "axis_seg"),
    c("restriction_sites", "restriction_site")
  )) {
    combined <- collect(pair[2])
    if (nrow(combined) > 0) primary[[pair[1]]] <- combined
  }
  repel_labels <- collect("gene_text_repel")
  fixed_labels <- collect("gene_text")
  all_gene_labels <- Filter(function(x) nrow(x) > 0,
                            list(fixed_labels, repel_labels))
  if (length(all_gene_labels)) {
    primary$gene_labels <- ggchord_rbind_fill(all_gene_labels)
  }

  palettes <- lapply(layouts, function(x) x$gene_pal)
  palettes <- Filter(function(x) !is.null(x) && length(x) > 0, palettes)
  if (length(palettes)) {
    pal <- do.call(c, unname(palettes))
    primary$gene_pal <- pal[!duplicated(names(pal), fromLast = TRUE)]
    orders <- unlist(lapply(layouts, function(x) x$final_gene_order),
                     use.names = FALSE)
    primary$final_gene_order <- unique(orders)
  }
  shape_palettes <- lapply(layouts, function(x) x$feature_shape_pal)
  shape_palettes <- Filter(
    function(x) !is.null(x) && length(x) > 0L, shape_palettes
  )
  if (length(shape_palettes)) {
    shape_pal <- do.call(c, unname(shape_palettes))
    primary$feature_shape_pal <- shape_pal[
      !duplicated(names(shape_pal), fromLast = TRUE)
    ]
    primary$feature_shape_order <- unique(unlist(lapply(
      layouts, function(x) x$feature_shape_order
    ), use.names = FALSE))
  }
  primary$extremes <- get_plot_extremes(
    allRibbon = ggchord_rbind_fill(Filter(Negate(is.null), list(
      primary$ribbon_polys, primary$link_lines, primary$restriction_sites
    ))),
    seqArcs = primary$seq_arcs,
    axisLines = primary$axis_lines,
    axisTicks = primary$axis_ticks,
    gene_polys = primary$gene_polys,
    gene_arrows = primary$gene_labels,
    seq_labels = primary$seq_labels_df,
    show_axis = primary$show_axis
  )

  if (!is.null(plot$ggchord$ref)) plot$ggchord$ref$layout <- primary
  plot$ggchord$layout <- primary
  primary
}


# ====================================================================
# Shared helpers used by ggplot_build.ggchord() and layout preparation.
# ====================================================================

#' Reconstruct a layer with the given data (and optional remapped mapping).
#'
#' LayerInstance objects cannot be cloned with \code{ggproto(NULL, .)}, so the
#' layer is rebuilt through \code{layer()} with the same geom/stat/mapping/params.
#' @noRd
reconstruct_layer <- function(lyr, data, mapping = NULL) {
  params <- c(lyr$geom_params, lyr$stat_params, lyr$aes_params)
  params <- params[!duplicated(names(params))]
  new <- ggplot2::layer(
    geom = lyr$geom, stat = lyr$stat, data = data,
    mapping = mapping %||% ggchord_effective_mapping(lyr), position = lyr$position,
    params = params,
    inherit.aes = lyr$inherit.aes,
    show.legend = lyr$show.legend,
    check.aes = FALSE
  )
  # Preserve the ggchord custom fields on the reconstructed layer
  for (fld in c(
    "ggchord_type", "ggchord_params", "ggchord_placeholder",
    "ggchord_layer_id", "ggchord_input_data", "ggchord_input_mapping",
    "ggchord_role_aes", "ggchord_resolved_input", "ggchord_theme_element",
    "ggchord_theme_components", "ggchord_input_transform", "ggchord_obstacle_provider"
  )) {
    if (!is.null(lyr[[fld]])) new[[fld]] <- lyr[[fld]]
  }
  new
}

#' Classify the ggchord layers of a plot by their ggchord_type marker
#' @noRd
classify_ggchord_layers <- function(plot) {
  idx <- list(seq = integer(0), ribbon = integer(0), gene_poly = integer(0),
              gene_text = integer(0), gene_text_repel = integer(0),
              gene_label_segment = integer(0),
              axis_line = integer(0), axis_seg = integer(0),
              axis_text = integer(0), axis = integer(0),
              gene_label_repel = integer(0), seq_label = integer(0),
              seq_region = integer(0), restriction_site = integer(0))
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    type <- lyr$ggchord_type %||% ""
    if (type %in% names(idx)) idx[[type]] <- c(idx[[type]], i)
  }
  idx
}
