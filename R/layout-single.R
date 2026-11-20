# Compute geometry for one sequence layout and its drawing layers.
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
  ribbon_data_override <- NULL
  gene_data_override <- NULL
  region_data_override <- NULL
  gene_label_params <- list()
  gene_repel_params <- list()
  axis_params   <- list()
  seq_label_params <- list()
  seq_region_params <- list()
  seq_layer_requested <- FALSE
  ribbon_layer_requested <- FALSE
  gene_geometry_layer_requested <- FALSE
  gene_label_layer <- FALSE
  gene_repel_layer <- FALSE

  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    pp <- lyr$ggchord_params
    if (is.null(pp)) next
    switch(pp$type,
      seq               = {
        seq_params <- pp
        seq_layer_requested <- TRUE
      },
      ribbon            = {
        ribbon_params <- pp
        ribbon_layer_requested <- TRUE
        ribbon_data_override <- ggchord_resolve_layer_input(
          lyr, data_list$ribbon_data
        )
      },
      gene              = {
        gene_params <- pp
        gene_geometry_layer_requested <- TRUE
        gene_data_override <- pp$gene_data_override %||%
          ggchord_resolve_layer_input(lyr, ggchord_annotation_fallback(data_list))

      },
      gene_label        = {
        gene_label_params <- pp
        gene_label_layer <- TRUE
        label_data <- ggchord_resolve_layer_input(lyr, ggchord_annotation_fallback(data_list))
        gene_label_params$label_text <- label_data$anno
        if (!isTRUE(gene_geometry_layer_requested)) {
          gene_data_override <- label_data
        }
      },
      gene_label_repel  = {
        gene_repel_params <- pp
        gene_repel_layer <- TRUE
        label_data <- ggchord_resolve_layer_input(lyr, ggchord_annotation_fallback(data_list))
        gene_repel_params$label_text <- label_data$anno
        if (!isTRUE(gene_geometry_layer_requested)) {
          gene_data_override <- label_data
        }
      },
      seq_label         = seq_label_params <- pp,
      seq_region        = {
        seq_region_params <- pp
        region_data_override <- ggchord_resolve_layer_input(lyr, ggchord_layer_fallback(lyr, data_list))
      }

    )
  }

  arrangement <- chord$sequence_layout
  # --- Process sequences ---
  seq_data <- data_list$seq_data
  sequence_unit <- ggchord_sequence_unit(seq_data, "ggchord layout")
  if (identical(sequence_unit, "aa")) {
    if (any(seq_data$length != floor(seq_data$length))) ggchord_stop(
      "ggchord(): protein lengths must be whole amino-acid counts")
    for (layer in plot$layers) {
      params <- layer$ggchord_params
      if (is.null(params)) next
      gene <- identical(params$type, "gene") && !isTRUE(params$is_feature)
      gene_label <- params$type %in% c("gene_label", "gene_label_repel") &&
        !isTRUE(params$is_feature_label)
      if (gene || gene_label || isTRUE(params$is_primer) ||
          isTRUE(params$is_primer_label) || identical(params$type, "restriction_site")) {
        ggchord_stop("ggchord(): gene, primer and restriction layers require bp coordinates; use generic feature layers for protein domains")
      }
    }
  }
  seqs     <- seq_data$accver
  lens     <- setNames(seq_data$length, seqs)

  if (!is.null(arrangement$order)) {
    if (!is.character(arrangement$order) || !length(arrangement$order) ||
        anyNA(arrangement$order) || anyDuplicated(arrangement$order) ||
        !all(arrangement$order %in% seqs)) {
      ggchord_stop("layout_sequences(): order must contain known sequence IDs from the displayed reference")
    }
    seqs <- arrangement$order
  }
  lens <- lens[seqs]
  n <- length(seqs)

  seq_labels    <- process_sequence_param(seq_params$seq_labels, seqs,
                                          "seq_labels", default_value = seqs)
  seqRadius     <- process_sequence_param(arrangement$radius, seqs,
                                          "layout_sequences(radius)", 1)
  orientation   <- process_sequence_param(arrangement$orientation, seqs,
                                          "layout_sequences(orientation)", 1)
  collinear <- isTRUE(plot$coordinates$ggchord_collinear)
  clockwise <- identical(plot$coordinates$reference_direction, "clockwise")
  if (clockwise) orientation <- -orientation
  seq_gap       <- if (collinear) ggchord_row_gap(plot, arrangement, seqs) else
    process_sequence_param(arrangement$gap, seqs, "angular_gap", 0)
  seq_curvature <- process_sequence_param(arrangement$curvature, seqs,
                                          "layout_sequences(curvature)", 1)

  if (isTRUE(plot$coordinates$ggchord_circular)) {
    source <- chord$display_instances$accver %||% seqs
    if (length(unique(source)) != 1L)
      ggchord_stop("coord_circular() requires one reference per panel")
  }

  seq_ring <- chord$display_instances$ring_id
  if (!is.null(seq_ring))
    seq_ring <- stats::setNames(seq_ring[match(seqs, chord$display_instances$key)], seqs)

  if (!is.numeric(seqRadius) || any(!is.finite(seqRadius)) || any(seqRadius <= 0)) {
    ggchord_stop("seq_radius must contain finite positive numbers")
  }
  if (!is.numeric(orientation) || any(!is.finite(orientation)) ||
      any(!orientation %in% c(-1, 1))) {
    ggchord_stop("layout_sequences(): orientation can only be 1 or -1")
  }
  if (!is.numeric(seq_gap) || any(!is.finite(seq_gap) | seq_gap < 0))
    ggchord_stop("Reference gaps must be finite and nonnegative")
  if (!is.numeric(seq_curvature) || any(!is.finite(seq_curvature))) {
    ggchord_stop("seq_curvature must contain finite numbers")
  }

  if (!is.null(seq_params$seq_colors)) {
    seq_colors <- process_sequence_param(seq_params$seq_colors, seqs, "seq_colors")
  } else {
    pal <- if (isTRUE(plot$coordinates$ggchord_circular)) "#252525" else
      chord_default_palette(n)
    seq_colors <- setNames(pal, seqs)
  }

  # --- Process ribbons ---
  ribbon_gap_auto <- is.null(ribbon_params$ribbon_gap)
  ribbonGap  <- process_sequence_param(ribbon_params$ribbon_gap %||% 0.15,
                                       seqs, "link_gap", 0.15)
  ribbon_offset <- ggchord_connection_offsets(ribbon_params, seqs,
    ggchord_path_normal_axis(plot$coordinates))
  if (!collinear) ribbonGap <- ribbonGap + ribbon_offset
  ribbon_color_scheme <- ribbon_params$ribbon_color_scheme %||% "pident"
  ribbon_color_by     <- ribbon_params$ribbon_color_by
  ribbon_color_limits <- ribbon_params$ribbon_color_limits
  ribbon_color_breaks <- ribbon_params$ribbon_color_breaks
  ribbon_color_name   <- ribbon_params$ribbon_color_name
  ribbon_alpha    <- ribbon_params$ribbon_alpha %||% 0.42
  ribbon_alpha_by <- ribbon_params$ribbon_alpha_by
  ribbon_alpha_range <- ribbon_params$ribbon_alpha_range %||% c(0.15, 0.9)
  ribbon_ctrl_pt  <- if (collinear) ribbon_params$ribbon_ctrl_point else
    ribbon_params$ribbon_ctrl_point %||% c(0, 0)
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
    ggchord_stop("link_gap must contain finite numbers")
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
  collinear_transform <- if (collinear) ggchord_collinear_transform(
    seqs, lens, orientation, chord$collinear_evidence %||% ribbon_data %||%
      data_list$ribbon_data, plot$coordinates,
    length_mode = arrangement$length_mode %||% "proportional"
  ) else NULL
  if (collinear) orientation <- collinear_transform$orientation
  ribbon_stat_data <- if (!is.null(ribbon_params$link_stat)) ribbon_data else NULL
  ribbon_stat_report <- if (is.null(ribbon_stat_data)) NULL else data.frame(
    source_rows = I(ggchord_link_source_members(ribbon_stat_data)))
  if (NROW(ribbon_data)) {
    strand <- ggchord_link_orientation(ribbon_data)
    reversed <- strand == "opposite"
    ribbon_data$.ggchord_relative_direction <- ifelse(reversed, "reverse", strand)
    q_minus <- ribbon_data$qstrand %in% "-" & length(ribbon_data$qstrand) > 0L
    if (!length(q_minus)) q_minus <- rep(FALSE, nrow(ribbon_data))
    s_minus <- xor(q_minus, reversed)
    ribbon_data[c("qstart", "qend")] <- list(
      ifelse(q_minus, pmax(ribbon_data$qstart, ribbon_data$qend),
        pmin(ribbon_data$qstart, ribbon_data$qend)),
      ifelse(q_minus, pmin(ribbon_data$qstart, ribbon_data$qend),
        pmax(ribbon_data$qstart, ribbon_data$qend)))
    ribbon_data[c("sstart", "send")] <- list(
      ifelse(s_minus, pmax(ribbon_data$sstart, ribbon_data$send),
        pmin(ribbon_data$sstart, ribbon_data$send)),
      ifelse(s_minus, pmin(ribbon_data$sstart, ribbon_data$send),
        pmax(ribbon_data$sstart, ribbon_data$send)))
  }
  display_rows <- if (collinear && "row_index" %in% names(collinear_transform$table))
    stats::setNames(collinear_transform$table$row_index, seqs) else NULL
  if (identical(sequence_unit, "aa") && !is.null(ribbon_data) && nrow(ribbon_data)) {
    columns <- c("qstart", "qend", "sstart", "send")
    ggchord_require_columns(ribbon_data, columns, "Protein alignments")
    for (column in columns) ggchord_conversion_coordinates(ribbon_data, column, "Protein alignments")

  }
  ggchord_check_interval_unit(ribbon_data, sequence_unit, "Ribbon layer")
  ribbon_selection <- ggchord_select_ribbons(ribbon_data, seqs, ribbon_params, collinear, display_rows,
    cyclic = isTRUE(plot$coordinates$ggchord_radial))
  ribbon_data <- ribbon_selection$data
  ribbon_ctrl_pt <- ggchord_link_control_rows(ribbon_ctrl_pt, ribbon_data$.ggchord_source_row)
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
  gene_w    <- gene_params$gene_width %||% 0.05
  gene_cs   <- gene_params$gene_color_scheme %||% "strand"
  gene_cols <- gene_params$gene_colors
  if (isTRUE(gene_params$is_feature) && !isTRUE(gene_params$user_feature_fill) && is.null(gene_cols) &&
      is.null(gene_params$feature_fill_mapping)) {
    gene_cols <- ggchord_feature_source_palette(
      gene_data_override %||% ggchord_annotation_fallback(data_list))
  }
  if (isTRUE(gene_params$is_feature) && !is.null(names(gene_cols)) && !is.null(gene_data_override$anno))
    gene_cols <- gene_cols[names(gene_cols) %in% gene_data_override$anno]
  ggchord_check_interval_unit(gene_data_override %||% ggchord_annotation_fallback(data_list), sequence_unit, "Feature layer")
  gene_ord  <- gene_params$gene_order
  # Gene label settings come from the dedicated geom_gene_label() layer, with
  # the legacy geom_feature(preset = 'gene') arguments as fallback.
  # The repel layer takes priority over the fixed label layer
  lbl <- if (gene_repel_layer) gene_repel_params else gene_label_params
  if (!is.null(lbl$position_explicit)) {
    lbl$position_supplied <- lbl$position_explicit
  }
  gene_ls   <- gene_label_layer || gene_repel_layer ||
    isTRUE(gene_params$show_label_override) ||
    isTRUE(gene_params$gene_label_show)
  label_theme_element <- if (isTRUE(lbl$is_feature_label)) {
    "ggchord.feature.label"
  } else "ggchord.gene.label"
  gene_lsz  <- lbl$gene_label_size %||%
    gene_params$label_size_override %||%
    gene_params$gene_label_size %||%
    ggchord_theme_text_size(plot, label_theme_element, 2.5)
  gene_lfamily <- lbl$gene_label_family %||% ""
  gene_lfontface <- lbl$gene_label_fontface %||% 1
  gene_llineheight <- lbl$gene_label_lineheight %||% 1.2
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
    if (identical(gene_repel_params$gene_label_layout, "feature")) {
      "feature"
    } else "radial"
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
  feature_label_external <- gene_repel_params$feature_label_external %||% TRUE
  gene_lrepel_fit    <- gene_repel_params$gene_label_fit %||% "wrap"
  gene_label_min_size <- gene_repel_params$gene_label_min_size %||% 1.5
  gene_lrepel_lines  <- gene_repel_params$gene_label_max_lines %||% 2L
  gene_lrepel_side   <- if (gene_repel_layer) {
    if (collinear) gene_repel_params$label_side_requested %||%
      gene_repel_params$gene_label_side %||% "outside" else
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
    ggchord_stop("size must be a finite positive number")
  }
  if (!is.null(gene_lwrap) && (!is.numeric(gene_lwrap) ||
      length(gene_lwrap) != 1 || !is.finite(gene_lwrap) || gene_lwrap < 0)) {
    ggchord_stop("label_wrap must be NULL or a finite non-negative number")
  }
  if (!is.numeric(gene_lrepel_maxov) || length(gene_lrepel_maxov) != 1 ||
      is.na(gene_lrepel_maxov) || gene_lrepel_maxov < 0) {
    ggchord_stop("max_overlaps must be a non-negative number or Inf")
  }

  metric_axis <- ggchord_path_normal_axis(plot$coordinates)
  normal_metric <- ggchord_units_per_mm(metric_axis)
  thickness <- gene_params$feature_mm
  if (!is.null(thickness)) {
    if (!is.null(chord$display_instances))
      thickness <- ggchord_instance_param(thickness, chord$display_instances)
    if (!is.null(names(thickness)) &&
        (anyDuplicated(names(thickness)) || any(!names(thickness) %in% seqs)))
      ggchord_stop("geom_feature(): invalid or unknown thickness keys")
    gene_w <- if (is.list(thickness)) unlist(thickness) else thickness
    gene_w <- gene_w * normal_metric
  }
  geneWidth <- process_gene_param(gene_w, seqs, "thickness",
    2 * normal_metric, FALSE)

  # Feature geometry is resolved before coordinate generation because these
  # values change the actual polygon, not only its appearance. A user-supplied
  # feature-shape scale is cloned and trained here so the layout and legend use
  # exactly the same category-to-shape mapping.
  gene_data_layout <- if (gene_geometry_layer_requested || gene_label_layer ||
      gene_repel_layer) {
    gene_data_override %||% ggchord_annotation_fallback(data_list)
  } else {
    NULL
  }
  if (!is.null(gene_data_layout) && isTRUE(lbl$is_feature_label)) {
    gene_data_layout$.feature_label_repeat <- rep(
      lbl$feature_label_repeat %||% "individual",nrow(gene_data_layout))
  }
  ribbon_obstacles <- plot$ggchord$obstacles %||% ggchord_empty_obstacles()
  if (!is.null(gene_data_layout) && !"anno" %in% names(gene_data_layout)) {
    gene_data_layout$anno <- rep(NA_character_, nrow(gene_data_layout))
  }
  feature_shape_pal <- NULL
  feature_shape_order <- NULL
  if ((isTRUE(gene_params$is_feature) || isTRUE(gene_params$gene_shape_override)) &&
      !is.null(gene_data_layout) &&
      nrow(gene_data_layout) > 0L) {
    raw_shape <- as.character(
      gene_data_layout$.feature_shape_raw %||%
        rep(gene_params$feature_shape %||% "arrow", nrow(gene_data_layout))
    )
    feature_shape_order <- unique(raw_shape)
    if (isTRUE(gene_params$feature_shape_mapped)) {
      shape_scale <- gene_params$shape_scale %||% plot$scales$get_scales("shape")
      if (!is.null(shape_scale)) {
        shape_scale <- shape_scale$clone()
        shape_scale$train(unique(c(chord$facet_feature_shape_order, raw_shape)))
        mapped_shape <- as.character(shape_scale$map(raw_shape))
      } else {
        ggchord_stop("geom_feature(): mapped shape requires an explicit shape scale")
      }
    } else {
      mapped_shape <- raw_shape
    }
    allowed_shape <- c(
      "arrow", "compact_arrow", "promoter_arrow", "primer_arrow",
      "marker", "block", "capped_line", "primer_arc", "chevron", "lollipop"
    )
    if (anyNA(mapped_shape) || any(!mapped_shape %in% allowed_shape)) {
      ggchord_stop(
        "geom_feature(): shape scale must produce feature geometry names"
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
    if (length(gene_params$head_args) && any(!mapped_shape %in% ggchord_feature_heads))
      ggchord_stop("geom_feature(): arrow head parameters require headed shapes")
    biological <- gene_data_layout$.feature_biological_strand
    unknown <- biological == "." & mapped_shape %in% c(ggchord_feature_heads, "chevron")
    if (identical(gene_params$unknown_strand, "error") && any(unknown))
      ggchord_stop("geom_feature(): directional shapes require a known strand")
    mapped_shape[unknown] <- "block"
    gene_data_layout$.feature_shape <- mapped_shape
    gene_data_layout$.feature_width_factor <- rep(1, length(mapped_shape))

  }
  if (!is.null(gene_data_layout) && nrow(gene_data_layout)) {
    if (!is.null(gene_params$head_width_mm))
      gene_data_layout$.feature_head_width <- gene_params$head_width_mm * normal_metric
    gene_data_layout$.feature_width <- vapply(seq_len(nrow(gene_data_layout)),
      function(i) {
        sid <- as.character(gene_data_layout$accver[i])
        strand <- as.character(gene_data_layout$strand[i])
        base_width <- geneWidth[[sid]][strand]
        base_width * as.numeric(
          gene_data_layout$.feature_width_factor[i] %||% 1
        )
      }, numeric(1))
  }
  if (!is.null(gene_data_layout)) {
    ggchord_check_feature_bounds(gene_data_layout, lens)
    active_position <- if (isTRUE(lbl$position_supplied)) {
      lbl$feature_position
    } else {
      gene_params$feature_position %||% ggplot2::position_identity()
    }
    gene_data_layout <- ggchord_apply_feature_position(
      gene_data_layout,
      active_position,
      seqs = seqs, lengths = lens,
      circular = isTRUE(plot$coordinates$ggchord_circular),
      units_per_mm = normal_metric,
      caller = "geom_feature()"
    )
    # Allocate on the complete entity table before selecting/reordering labels;
    # a label subset must not cause the remaining entities to change lanes.
    if (!is.null(chord$matched_feature_rows)) {
      gene_data_layout <- gene_data_layout[chord$matched_feature_rows$rows, , drop = FALSE]
      # Keep the geometry-derived lane, width, shape and label contrast while
      # measuring and drawing the text supplied by this label layer.
      gene_data_layout$.anchor_label_text <- lbl$label_text
      if (!is.null(names(gene_cols)))
        gene_cols <- gene_cols[names(gene_cols) %in% gene_data_layout$anno]
      rownames(gene_data_layout) <- NULL
    }
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
  position_scale <- plot$scales$get_scales("seq_position")
  if (is.null(position_scale) && !is.null(chord$display_instances$window_start))
    position_scale <- scale_seq_position(name = NULL)
  seq_guide <- position_scale$ggchord_axis_guide
  axis_theme <- ggchord_plot_settings(plot)$axis
  guide_theme <- seq_guide$theme
  if (!is.null(guide_theme[["axis.ticks.length"]]))
    axis_theme$ticks.length <- guide_theme[["axis.ticks.length"]]
  if (is.numeric(seq_guide$angle)) axis_theme$text.orientation <- seq_guide$angle
  side_sign <- if (collinear) -1 else 1
  if (identical(seq_guide$side, "negative")) side_sign <- -side_sign
  axis_theme$gap <- side_sign * if (as.numeric(axis_theme$gap) == 0)
    grid::unit(.8, "mm") else axis_theme$gap

  show_axis <- !isTRUE(axis_theme$hidden)
  # Pass physical distances as inches. compute_chord_layout() converts them
  # against the actual curved sequence span once that geometry is available.
  axis_unit_data <- ggchord_unit_inches
  axisGap    <- process_sequence_param(axis_unit_data(axis_theme$gap),
                                       seqs, "axis.gap", 0.04)
  axisMaj    <- process_sequence_param(axis_params$axis_tick_major_number %||% 5,
                                       seqs, "axis_tick_major_number", 5)
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
  axisLabelHide <- position_scale$ggchord_check_overlap %||% isTRUE(axis_theme$text.check.overlap)
  axisLabelOri <- process_axis_orientation(
    axis_theme$text.orientation, seqs
  )
  seq_style <- seq_params$seq_style %||% "auto"
  if (identical(seq_style, "auto")) {
    seq_style <- if (isTRUE(plot$coordinates$ggchord_circular)) "double" else "single"
  }
  seqBackboneOuter <- switch(
    seq_style,
    double = (seq_params$seq_backbone_gap %||% .025) / 2,
    band = (seq_params$seq_backbone_width %||% .035) / 2,
    0
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
  if (!is.null(position_scale)) {
    axis_breaks <- axis_minor_breaks <- axis_labels <- setNames(
      vector("list", length(seqs)), seqs
    )
    for (id in seqs) {
      sc <- position_scale$clone()
      instances <- chord$display_instances
      provenance <- ggchord_focus_source(seq_data)
      at <- match(id, seq_data$accver)
      instance <- if (is.null(instances)) NA_integer_ else match(id, instances$key)
      source_axis <- identical(position_scale$ggchord_position_coordinates, "source")
      origin <- if (source_axis) provenance$origin[at] else
        if (!is.na(instance) && !is.null(instances$window_start)) instances$window_start[instance] else 1
      source_id <- if (!is.na(instance)) instances$accver[instance] else id
      for (field in names(position_scale$ggchord_axis_spec)) {
        values <- position_scale$ggchord_axis_spec[[field]]
        if (is.list(values) && !ggplot2::is_waiver(values)) {
          if (is.null(names(values)) || anyDuplicated(names(values)))
            ggchord_stop("Axis ", field, " must use distinct source or instance keys")
          if (all(names(values) %in% c("source", "instance"))) {
            if (anyDuplicated(names(values$source)) || anyDuplicated(names(values$instance)) ||
                any(!names(values$source) %in% chord$complete_reference_ids) ||
                any(!names(values$instance) %in% (chord$complete_display_instances$instance_id %||% chord$complete_reference_ids)))
              ggchord_stop("Unknown or repeated axis ", field, " key")
            values <- values$instance[[id]] %||% values$source[[source_id]] %||% waiver()
          } else {
            if (any(!names(values) %in% chord$complete_reference_ids))
              ggchord_stop("Unknown source key in axis ", field)
            values <- values[[source_id]] %||% waiver()
          }
        }
        sc[[field]] <- values
      }
      domain <- c(origin, origin + lens[[id]] - 1)
      domain <- suppressWarnings(sc$transform(domain))
      sc$train(domain[is.finite(domain)])
      get_labels <- function(breaks) {
        labels <- sc$labels
        if (!(is.expression(labels) || is.function(labels))) return(sc$get_labels(breaks))
        if (is.function(labels)) labels <- labels(sc$get_transformation()$inverse(breaks))
        if (!is.null(labels) && length(labels) != length(breaks))
          ggchord_stop("scale_seq_position(): labels must match the number of breaks")
        labels
      }
      br <- sc$get_breaks()
      raw_br <- sc$get_transformation()$inverse(br) %||% numeric()
      raw_br <- raw_br - origin + 1
      automatic_breaks <- inherits(sc$breaks, "waiver")
      if (automatic_breaks) {
        domain_limits <- range(sc$get_transformation()$inverse(sc$get_limits()))
        endpoints <- c(1,lens[[id]])
        source_endpoints <- endpoints + origin - 1
        endpoints <- endpoints[source_endpoints>=domain_limits[1] & source_endpoints<=domain_limits[2]]
        raw_br <- sort(unique(c(raw_br[is.finite(raw_br) & raw_br>=1 & raw_br<=lens[[id]]],endpoints)))
        source_breaks <- raw_br + origin - 1
        br <- sc$transform(source_breaks)
      }
      break_labels <- get_labels(br) %||% rep(NA_character_, length(raw_br))
      keep <- is.finite(raw_br) & raw_br >= 0 & raw_br <= lens[[id]]
      br <- br[keep]
      values <- raw_br[keep]
      attr(values,"ggchord_automatic_breaks") <- automatic_breaks
      axis_breaks[id] <- list(values)
      minor <- sc$get_breaks_minor()
      minor <- sc$get_transformation()$inverse(minor) %||% numeric()
      minor <- minor - origin + 1
      axis_minor_breaks[id] <- list(minor[
        is.finite(minor) & minor >= 0 & minor <= lens[[id]]
      ])
      axis_labels[id] <- list(break_labels[keep])
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
      if (collinear) {
        process_sequence_param(1, seqs, "seq_label_hjust", 1)
      } else if (identical(seq_label_orientation, "arc")) {
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
    plot$coordinates$rotation
  } else {
    global$rotation
  }

  layout_strategy <- if (collinear) {
    function(...) compute_collinear_layout(..., gene_label_text = lbl$label_text,
      seq_label_hjust_explicit = !is.null(seq_label_params$seq_label_hjust),
      seq_label_vjust_explicit = !is.null(seq_label_params$seq_label_vjust),
      collinear_transform = collinear_transform,
      collinear_crossrow_styles = ggchord_collinear_row_styles(plot$coordinates),
      collinear_radial = isTRUE(plot$coordinates$ggchord_radial),
      collinear_right_normal = identical(plot$coordinates$collinear_seq_axis, "vertical"),
      collinear_vertical = isTRUE(plot$coordinates$ggchord_sankey) ||
        identical(plot$coordinates$collinear_seq_axis, "vertical"))
  } else function(..., ribbon_offset) compute_chord_layout(..., ring = seq_ring)
  layout <- layout_strategy(
    seqs = seqs, lens = lens, seq_labels = seq_labels,
    length_mode = arrangement$length_mode %||% "proportional",
    seq_colors = seq_colors, seqRadius = seqRadius,
    seq_curvature = seq_curvature, orientation = orientation,
    seq_gap = seq_gap, clockwise = clockwise,
    ribbon_data = ribbon_data, ribbonGap = ribbonGap,
    ribbon_offset = if (isTRUE(plot$coordinates$ggchord_radial) ||
      isTRUE(plot$coordinates$ggchord_sankey)) ribbon_offset * 0 else ribbon_offset,
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
    geneWidth = geneWidth,
    arrow_head_length = (gene_params$head_length_mm %||% 2) *
      ggchord_units_per_mm(if (is.null(metric_axis)) NULL else
        if (metric_axis == "x") "y" else "x"),
    arrow_head_width = gene_params$arrow_head_width %||% 1,
    arrow_head_style = gene_params$arrow_head_style %||% "shouldered",
    short_feature = gene_params$short_feature %||% "auto",
    circular = isTRUE(plot$coordinates$ggchord_circular),
    geneLabelRadialOffset = geneLabelRadialOffset,
    geneLabelCircumOffset = geneLabelCircumOffset,
    geneLabelCircumLimit = geneLabelCircumLimit,
    geneLabelRotation = geneLabelRotation,
    gene_label_show = gene_ls, gene_label_size = gene_lsz,
    gene_label_family = gene_lfamily,
    gene_label_fontface = gene_lfontface,
    gene_label_lineheight = gene_llineheight,
    gene_label_wrap = gene_lwrap,
    gene_label_fit = gene_lrepel_fit,
    gene_label_min_size = gene_label_min_size,
    gene_label_max_lines = gene_lrepel_lines,
    gene_label_orientation = gene_lorientation,
    gene_label_overlap = gene_loverlap,
    gene_label_repel_layer = gene_lrepel_layer,
    gene_label_repel_max_overlaps = gene_lrepel_maxov,
    gene_label_layout = gene_lrepel_layout,
    feature_label_external = feature_label_external,
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
    axis_text_style = {
      style <- ggchord_theme_element(plot, "ggchord.axis.text")
      text_style <- if(inherits(style,"element_text")) list(family=style@family,
        fontface=style@face,lineheight=style@lineheight) else list()
      c(text_style, list(offset_explicit = isTRUE(axis_theme$text.offset.explicit),
        allow_clearance = axisLabelHide))
    },
    axisLabelOrientation = axisLabelOri,
    seqBackboneOuter = seqBackboneOuter,
    axis_breaks = axis_breaks,
    axis_minor_breaks = axis_minor_breaks,
    axis_labels = axis_labels,
    axis_label_hide_overlaps = axisLabelHide,
    show_axis = show_axis,
    rotation = coord_rotation, debug = global$debug,
    geometry_cache = geometry_cache
  )

  layout$seq_ring <- seq_ring
  layout$length_mode <- arrangement$length_mode %||% "proportional"
  layout$seq_ring_radius <- if (is.null(seq_ring)) NULL else seqRadius
  layout$ribbon_stat_data <- ribbon_stat_data
  layout$ribbon_stat_report <- ribbon_stat_report
  layout$ribbon_selection <- ribbon_selection$report
  if (collinear) {
    source_rows <- match(seqs, as.character(seq_data$accver))
    layout$sequence_layout$source_origin <- if (".source_origin" %in% names(seq_data)) {
      seq_data$.source_origin[source_rows]
    } else if (".focus_start" %in% names(seq_data)) seq_data$.focus_start[source_rows] else 1
    layout$sequence_layout$source_accver <- if (".source_accver" %in% names(seq_data)) {
      seq_data$.source_accver[source_rows]
    } else seqs
    layout$sequence_layout$source_end <- layout$sequence_layout$source_origin + unname(lens[seqs]) - 1
    layout$sequence_layout$source_length <- if (".source_length" %in% names(seq_data))
      seq_data$.source_length[source_rows] else unname(lens[seqs])
    for (side in c("left", "right")) {
      field <- paste0(".truncated_", side)
      layout$sequence_layout[[paste0("truncated_", side)]] <- if (field %in% names(seq_data)) {
        seq_data[[field]][source_rows]
      } else FALSE
    }
  }
  layout$backbone_bounds <- if (collinear) data.frame(
    accver = seqs, center = unname(layout$sequence_reference$y[seqs]),
    lower = unname(layout$sequence_reference$y[seqs]) - seqBackboneOuter,
    upper = unname(layout$sequence_reference$y[seqs]) + seqBackboneOuter,
    stringsAsFactors = FALSE
  ) else data.frame(
    accver = names(seqRadius), center = as.numeric(seqRadius),
    inner = as.numeric(seqRadius) - seqBackboneOuter,
    outer = as.numeric(seqRadius) + seqBackboneOuter,
    stringsAsFactors = FALSE
  )
  if (!is.null(seq_ring) && length(layout$seq_arcs)) {
    layout$seq_arcs <- lapply(names(layout$seq_arcs), function(id) {
      arc <- layout$seq_arcs[[id]]
      arc$seq_ring <- unname(seq_ring[[id]])
      arc
    }) |>
      stats::setNames(names(layout$seq_arcs))
  }

  layout$sequence_unit <- sequence_unit
  if (is.data.frame(layout$sequence_layout)) {
    layout$sequence_layout$sequence_unit <- sequence_unit
    layout$sequence_layout$units_per_position <- layout$sequence_layout$units_per_bp
  }
  position_scale <- plot$scales$get_scales("seq_position")
  selection <- ggchord_axis_selection(position_scale$ggchord_axis_sequences,
    chord$complete_reference_ids %||% seqs, chord$complete_display_instances)
  if (!is.null(selection)) layout <- ggchord_filter_sequence_axes(layout, intersect(selection, seqs))
  layout
}
