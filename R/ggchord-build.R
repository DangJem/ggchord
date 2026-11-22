# ====================================================================
# ggplot_build.ggchord: compute the layout, inject data into cloned
# layers, add scales, and set the coordinate system.  Everything is
# driven by the plot object itself, so the plot can be printed, saved
# with ggsave(), or built with ggplot_build() any number of times and
# in any order, without cross-talk between plots.
# ====================================================================

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
    seq_break = layout$seq_breaks %||% fallback,
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
  plot <- ggchord_render_inputs(plot)
  if (isTRUE(plot$coordinates$ggchord_dotplot))
    return(ggchord_cache_layout(plot, ggchord_dotplot_geometry(ggchord_register_layers(plot))))
  plot <- ggchord_prepare_geometry_inputs(plot)
  ggchord_cache_layout(plot, compute_chord_geometry_prepared(plot))
}

ggchord_prepare_geometry_inputs <- function(plot, prepare_instances = TRUE) {
  if (any(vapply(plot$layers, function(x) identical(x$ggchord_type, "alignment"), logical(1))))
    ggchord_stop("geom_alignment_segment() requires coord_dotplot()")
  plot <- ggchord_register_layers(plot)
  for (layer in plot$layers) {
    if (identical(layer$ggchord_params$seq_style, "auto")) {
      circular <- isTRUE(plot$coordinates$ggchord_circular)
      if (!circular && isTRUE(layer$ggchord_params$seq_thickness_supplied))
        ggchord_stop("geom_seq(): line uses linewidth, thickness applies to double or band")
      layer$ggchord_params$seq_style <- if (circular) "double" else "single"
    }
    if (!is.null(layer$ggchord_type) && is.null(layer$ggchord_params$native_geom) &&
        is.null(layer$ggchord_params$ribbon_stat))
      layer$ggchord_params$visual_source <- new.env(parent = emptyenv())
  }
  sequences <- ggchord_sequence_reference_data(plot)
  track_scales <- lapply(plot$scales$scales, function(scale) scale$ggchord_track_id)
  unknown <- setdiff(unlist(track_scales), c("default", names(plot$ggchord$track_bands)))
  if (length(unknown)) ggchord_stop("scale_track_value(): unknown track '", unknown[1L], "'")
  for (scope in ggchord_scale_scopes(plot)) {
    if (!"shape" %in% scope$aesthetics) next
    for (i in scope$selected) if (plot$layers[[i]]$ggchord_type %in% c("gene_poly", "restriction_site"))
      plot$layers[[i]]$ggchord_params$shape_scale <- scope$scale
  }
  for (layer in plot$layers) {
    type <- layer$ggchord_params$type %||% ""
    if (type %in% c("seq_region", "link", "ribbon", "gene", "gene_label", "gene_label_repel", "restriction_site"))
      layer$ggchord_params$reference_ids <- sequences$accver
    if (type %in% c("seq_region", "link", "restriction_site")) {
      input <- ggchord_layer_source(layer, ggchord_layer_fallback(layer, plot$ggchord$data))
      unit <- ggchord_sequence_unit(sequences)
      caller <- ggchord_layer_input_context(layer)
      ggchord_check_interval_unit(input, unit, caller)
    }
  }
  plot <- ggchord_prepare_facet_sources(plot)
  for (layer in plot$layers) {
    if (!is.null(layer$ggchord_params$source_stat)) {
      layer$ggchord_params$source_stat_facet <- plot$ggchord$source_facet_setup
      if (is.null(layer$ggchord_params$source_stat_cache))
        layer$ggchord_params$source_stat_cache <- new.env(parent = emptyenv())
    }
  }
  plot <- ggchord_prepare_track_layers(plot, sequences)
  plot <- ggchord_prepare_link_statistics(plot, sequences)
  plot <- ggchord_bind_domain_labels(plot)
  if (prepare_instances) plot <- ggchord_prepare_instances(plot, sequences)
  plot$scales <- plot$scales$clone()
  plot$scales$scales <- Filter(function(s) is.null(s$ggchord_managed),
                               plot$scales$scales)
  plot$layers <- lapply(plot$layers, ggchord_inherit_ribbon_linetype,
    mapping = plot$mapping)
  plot$layers <- lapply(plot$layers, ggchord_resolve_style_defaults,
    mapping = plot$mapping)
  plot$layers <- lapply(plot$layers, ggchord_bind_visual_mapping, plot = plot)
  ggchord_prepare_site_size_scales(plot)
}

compute_chord_geometry_prepared <- function(plot) {
  if (!inherits(plot$facet, "FacetNull")) return(ggchord_facet_layout(plot))
  ggchord_validate_coord_contract(plot)
  if (isTRUE(plot$ggchord$sequence_layout$auto_row_gap)) {
    # Resolve once using every layer in this panel. Each layer-specific
    # geometry pass must reuse the same reference paths and spacing.
    plot$ggchord$sequence_layout$row_gap <- ggchord_auto_row_gap(plot)
  }
  chord <- plot$ggchord
  if (is.null(chord)) {
    ggchord_stop("Not a valid ggchord object: no data stored on the plot")
  }

  if (isTRUE(plot$coordinates$ggchord_collinear)) {
    evidence <- lapply(plot$layers, function(layer) {
      if (!identical(layer$ggchord_params$type, "ribbon")) return(NULL)
      ggchord_resolve_layer_input(layer, chord$data$ribbon_data)
    })
    evidence <- Filter(function(d) is.data.frame(d) && nrow(d), evidence)
    plot$ggchord$collinear_evidence <- if (length(evidence)) {
      unique(ggchord_rbind_fill(evidence))
    } else chord$data$ribbon_data
  }

  plot <- ggchord_resolve_feature_geometry_context(plot)
  geometry_cache <- new.env(parent = emptyenv())
  avoid_requested <- any(vapply(plot$layers, function(layer) {
    params <- layer$ggchord_params
    !is.null(params) && params$type %in% c("ribbon", "link") &&
      !identical(params$link_avoid %||% "none", "none") &&
      is.null(params$ribbon_gap %||% params$link_gap)
  }, logical(1)))
  has_tracks <- any(vapply(plot$layers, function(layer) identical(layer$ggchord_type, "track"), logical(1)))
  plot$ggchord$obstacles <- if (avoid_requested || has_tracks) {
    ggchord_collect_obstacles(plot, geometry_cache)
  } else {
    ggchord_empty_obstacles()
  }
  primary <- compute_chord_geometry_single(plot, geometry_cache)
  if (isTRUE(plot$coordinates$ggchord_radial)) primary$radial_angles <- setNames(
    (plot$coordinates$radial_start_angle + (if (plot$coordinates$radial_clockwise) -1 else 1) *
      360 * (seq_along(primary$seqs) - 1) / length(primary$seqs)) * pi / 180, primary$seqs)
  primary$display_instances <- plot$ggchord$display_instances
  primary$instance_routing <- plot$ggchord$instance_routing
  primary$instance_omissions <- plot$ggchord$instance_omissions
  if (!is.null(primary$display_instances)) {
    instances <- primary$display_instances
    entities <- unique(instances$accver)
    first <- match(entities, instances$accver)
    primary$legend_seqs <- entities
    primary$legend_seq_colors <- setNames(unname(primary$seq_colors[instances$key[first]]), entities)
    primary$legend_seq_labels <- setNames(entities, entities)
  }
  primary$obstacles <- plot$ggchord$obstacles
  ids <- vapply(plot$layers, function(x) {
    if (is.null(x$ggchord_type)) "" else x$ggchord_layer_id %||% ""
  }, character(1))
  groups <- split(which(nzchar(ids)), ids[nzchar(ids)])

  group_type <- vapply(groups, function(idx) {
    types <- vapply(idx, function(i) plot$layers[[i]]$ggchord_params$type %||% "",
                    character(1))
    main <- types[types %in% c(
      "seq", "ribbon", "link", "gene", "gene_label", "gene_label_repel", "axis",
      "seq_label", "seq_region", "seq_break", "seq_scale", "track", "restriction_site", "domain_label"
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
    if (main_type == "track") next
    needs_own <- nzchar(main_type) && (type_counts[[main_type]] > 1 ||
      main_type %in% c("gene_label", "gene_label_repel"))
    # A label layer owns its Position independently of the entity polygons.
    # Do not reuse a primary layout whose active Position came from a label.
    if (main_type == "gene" &&
        any(group_type %in% c("gene_label", "gene_label_repel"))) needs_own <- TRUE
    sub_layout <- primary
    if (isTRUE(needs_own)) {
      deps <- seq_dep
      matched_feature_rows <- NULL
      if (main_type %in% c("gene_label", "gene_label_repel")) {
        matched_feature_rows <- plot$layers[[tail(idx, 1L)]]$ggchord_params$matched_feature
        if (!is.null(matched_feature_rows))
          deps <- c(deps, groups[[matched_feature_rows$layer_id]])
      }
      if (main_type == "ribbon") deps <- c(deps, gene_geometry_dep)
      sub_plot <- plot
      sub_plot$ggchord$matched_feature_rows <- matched_feature_rows
      sub_plot$layers <- plot$layers[sort(unique(c(deps, idx)))]
      sub_layout <- compute_chord_geometry_single(sub_plot, geometry_cache)
    }
    if (main_type == "link") {
      link_layer <- plot$layers[[idx[1L]]]
      link_input <- ggchord_resolve_layer_input(link_layer, ggchord_layer_fallback(link_layer, chord$data))
      sub_layout$link_lines <- ggchord_attach_input_columns(
        ggchord_link_geometry(link_input, link_layer$ggchord_params, primary), link_input)
    }
    if (main_type == "seq_break") {
      break_layer <- plot$layers[[idx[1L]]]
      sub_layout$seq_breaks <- ggchord_seq_break_geometry(
        ggchord_resolve_layer_input(break_layer, chord$data$seq_data),
        break_layer$ggchord_params, primary)
    }
    if (main_type == "seq_scale") {
      scale_layer <- plot$layers[[idx[1L]]]
      sub_layout$sequence_scales <- ggchord_seq_scale_geometry(
        ggchord_resolve_layer_input(scale_layer, chord$data$seq_data),
        scale_layer$ggchord_params, primary)
    }
    if (main_type == "restriction_site") {
      layer <- plot$layers[[idx[1L]]]
      source <- ggchord_resolve_layer_input(layer, ggchord_layer_fallback(layer, chord$data))
      sub_layout$restriction_sites <- ggchord_site_geometry(layer, source, primary, plot)
    }
    layouts[[id]] <- sub_layout
    registry[[id]] <- list()
    inputs[[id]] <- list()
    for (i in idx) {
      lyr <- plot$layers[[i]]
      component <- lyr$ggchord_type
      registry[[id]][[component]] <- if (identical(component, "seq")) {
        ggchord_seq_style_geometry(
          sub_layout$seq_arcs[intersect(names(sub_layout$seq_arcs), as.character(
            ggchord_resolve_layer_input(lyr, chord$data$seq_data)$accver))],
          modifyList(lyr$ggchord_params, list(
            circular = isTRUE(sub_layout$circular),
            geometry_reference_ids = sub_layout$seqs,
            normal_axis = ggchord_path_normal_axis(plot$coordinates)
          ))
        )
      } else if (identical(component, "seq_scale")) {
        sub_layout$sequence_scales
      } else {
        ggchord_layout_component(
          sub_layout, component, lyr$ggchord_placeholder %||% data.frame()
        )
      }
      fallback <- ggchord_layer_fallback(lyr, chord$data)
      inputs[[id]][[component]] <- if (
          component == "ribbon" && !is.null(sub_layout$ribbon_stat_data)) {
        sub_layout$ribbon_stat_data
      } else {
        ggchord_resolve_layer_input(lyr, fallback)
      }
      registry[[id]][[component]] <- ggchord_attach_input_columns(
        registry[[id]][[component]], inputs[[id]][[component]]
      )
      anchor <- lyr$ggchord_params$matched_feature$layer_id
      if (!is.null(anchor) && is.data.frame(registry[[id]][[component]]))
        registry[[id]][[component]]$anchor_layer <- rep(anchor, nrow(registry[[id]][[component]]))
      if (identical(component, "gene_label_repel") &&
          isTRUE(lyr$ggchord_params$is_primer_label) &&
          identical(lyr$ggchord_params$gene_label_side, "inside") &&
          nrow(registry[[id]][[component]])) {
        registry[[id]][[component]]$.ggchord_primer_inside <- TRUE
      }
      if (identical(component, "restriction_site") &&
          identical(lyr$ggchord_params$label_side, "inside") &&
          nrow(registry[[id]][[component]])) {
        registry[[id]][[component]]$.ggchord_restriction_inside <- TRUE
      }
    }
  }
  primary$layer_inputs <- inputs
  circular_resolved <- FALSE
  external_resolved <- FALSE
  if (has_tracks && isTRUE(plot$coordinates$ggchord_circular)) {
    resolved_annotations <- plot$coordinates$resolve_annotation_registry(registry, primary)
    registry <- resolved_annotations$layer_geometry
    primary$circular_annotation_registry <- resolved_annotations$annotation_registry
    circular_resolved <- TRUE
  }
  if (has_tracks && !isTRUE(primary$collinear) && !circular_resolved) {
    registry <- ggchord_share_external_annotations(registry, primary)
    external_resolved <- TRUE
  }
  if (has_tracks) {
    tracks <- ggchord_layout_tracks(plot, primary, registry, inputs)
    primary <- tracks$layout; registry <- tracks$registry; inputs <- tracks$inputs
    axes <- ggchord_track_axes(primary, registry, plot)
    registry <- axes$registry; primary$track_axes <- axes$axes
    primary$obstacles <- ggchord_rbind_fill(list(primary$obstacles, primary$track_obstacles))
    # Only requested automatic avoidance consumes the newly registered bands.
    for (id in names(groups)) {
      idx <- groups[[id]]; layer <- plot$layers[[idx[1L]]]
      params <- layer$ggchord_params
      if (params$type == "seq_scale") {
        registry[[id]]$seq_scale <- ggchord_attach_input_columns(
          ggchord_seq_scale_geometry(inputs[[id]]$seq_scale, params, primary), inputs[[id]]$seq_scale)
        next
      }
      if (!params$type %in% c("ribbon", "link") ||
          identical(params$link_avoid %||% "none", "none") ||
          !is.null(params$ribbon_gap %||% params$link_gap)) next
      if (params$type == "ribbon") {
        sub <- plot; sub$ggchord$obstacles <- primary$obstacles
        sub$layers <- plot$layers[sort(unique(c(seq_dep, gene_geometry_dep, idx)))]
        updated <- compute_chord_geometry_single(sub, geometry_cache)
        layouts[[id]] <- updated
        registry[[id]]$ribbon <- ggchord_attach_input_columns(updated$ribbon_polys, inputs[[id]]$ribbon)
      } else registry[[id]]$link <- ggchord_attach_input_columns(
        ggchord_link_geometry(inputs[[id]]$link, params, primary), inputs[[id]]$link)
    }
  }
  if (isTRUE(plot$coordinates$ggchord_circular) &&
      is.function(plot$coordinates$resolve_annotation_registry) && !circular_resolved) {
    resolved_annotations <- plot$coordinates$resolve_annotation_registry(
      registry, primary
    )
    registry <- resolved_annotations$layer_geometry
    primary$circular_annotation_registry <-
      resolved_annotations$annotation_registry
  } else if (isTRUE(primary$collinear)) {
    if (!isTRUE(primary$curved_backbone) &&
        !isTRUE(plot$coordinates$ggchord_radial) && !isTRUE(plot$coordinates$ggchord_sankey) &&
        !identical(plot$coordinates$collinear_seq_axis, "vertical")) {
      registry <- ggchord_collinear_share_labels(registry, primary)
      primary$collinear_label_diagnostics <- attr(registry, "collinear_label_diagnostics")
      attr(registry, "collinear_label_diagnostics") <- NULL
    }
  } else if (!circular_resolved && !external_resolved) {
    registry <- ggchord_share_external_annotations(registry, primary)
  }
  for (id in names(registry)) for (component in names(registry[[id]])) {
    d <- registry[[id]][[component]]
    if (identical(component,"gene_label_repel") && isTRUE(primary$circular) &&
        is.data.frame(d) && all(c("feature_label_mode","anchor_x","anchor_y") %in% names(d))) {
      text <- d[d$.component %in% "text",,drop=FALSE]
      backbone <- unname(primary$seqRadius[as.character(text$accver)])
      internal_exit <- text$feature_label_mode %in% "external" &
        sqrt(text$anchor_x^2+text$anchor_y^2) < backbone -
          .1*(primary$text_units_per_inch %||% .35)/25.4
      # Perimeter solvers rebuild segment vertices without text anchor fields.
      # Resolve style from the final owning label, never from a segment's NAs.
      d$.feature_internal_exit <- internal_exit[match(d$source_row,text$source_row)]
    }
    if (is.data.frame(d) && "label_size_requested" %in% names(d))
      d <- ggchord_label_fit_metadata(d)
    registry[[id]][[component]] <- d
  }
  primary$layer_geometry <- registry
  primary$layer_inputs <- inputs
  primary$layer_layouts <- layouts
  connection_reports <- lapply(names(layouts), function(id) {
    if (!identical(group_type[[id]], "ribbon")) return(NULL)
    d <- layouts[[id]]$ribbon_selection
    if (is.null(d) || !nrow(d)) return(NULL)
    d$layer_id <- id
    d
  })
  primary$connection_report <- ggchord_rbind_fill(connection_reports) %||% data.frame()

  # Plot limits must see all independent layers, not only the compatibility
  # fields in the primary layout.
  collect <- function(component, columns = NULL) {
    values <- lapply(registry, `[[`, component)
    values <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, values)
    if (!is.null(columns)) values <- lapply(values, `[`, columns)
    if (length(values)) ggchord_rbind_fill(values) else data.frame()
  }
  for (pair in list(
    c("ribbon_polys", "ribbon"), c("link_lines", "link"), c("gene_polys", "gene_poly"),
    c("gene_labels", "gene_text"), c("gene_label_segments", "gene_label_segment"),
    c("seq_labels_df", "seq_label"), c("region_polys", "seq_region"),
    c("seq_breaks", "seq_break"),
    c("sequence_scales", "seq_scale"),
    c("axis_lines", "axis_line"), c("axis_ticks", "axis_seg"),
    c("restriction_sites", "restriction_site"),
    c("seq_center_label", "seq_center_label")
  )) {
    combined <- collect(pair[2])
    if (nrow(combined) > 0) primary[[pair[1]]] <- combined
  }
  repel_labels <- collect("gene_text_repel")
  fixed_labels <- collect("gene_text")
  callout_labels <- collect("gene_label_repel")
  if (isTRUE(primary$collinear)) {
    primary$gene_labels <- data.frame()
    segments <- list(collect("gene_label_segment"))
    if (nrow(callout_labels) && ".component" %in% names(callout_labels)) {
      segments[[2L]] <- callout_labels[
        callout_labels$.component %in% "segment", , drop = FALSE
      ]
    }
    segments <- Filter(function(d) nrow(d) > 0L, segments)
    primary$gene_label_segments <- if (length(segments)) {
      ggchord_rbind_fill(segments)
    } else data.frame()
  }
  if (nrow(callout_labels) && ".component" %in% names(callout_labels)) {
    callout_labels <- callout_labels[
      callout_labels$.component %in% "text", , drop = FALSE
    ]
  }
  if (nrow(callout_labels)) {
    callout_labels$text <- if ("feature_label" %in% names(callout_labels)) {
      ifelse(!is.na(callout_labels$feature_label) &
        nzchar(callout_labels$feature_label), callout_labels$feature_label,
        callout_labels$label)
    } else callout_labels$label
    callout_labels$text_angle <- callout_labels$angle %||% 0
  }
  all_gene_labels <- Filter(function(x) nrow(x) > 0,
                            list(fixed_labels, repel_labels, callout_labels))
  if (length(all_gene_labels)) {
    primary$gene_labels <- ggchord_rbind_fill(all_gene_labels)
  }

  palettes <- lapply(layouts, function(x) x$gene_pal)
  palettes <- Filter(function(x) !is.null(x) && length(x) > 0, palettes)
  if (length(palettes)) {
    pal <- do.call(c, unname(palettes))
    primary$gene_pal <- pal[!duplicated(names(pal))]
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
  if (isTRUE(primary$curved_backbone)) {
    primary <- ggchord_bend_local_layout(primary)
    registry <- primary$layer_geometry
  }
  primary$extremes <- get_plot_extremes(
    allRibbon = ggchord_rbind_fill(Filter(Negate(is.null), list(
      primary$ribbon_polys, primary$link_lines, primary$restriction_sites,
      if (isTRUE(primary$collinear) && nrow(primary$region_polys)) primary$region_polys,
      if (isTRUE(primary$collinear) && nrow(collect("seq"))) collect("seq"),
      if (isTRUE(primary$collinear) && nrow(collect("seq_break"))) collect("seq_break")
      , collect("track", c("x", "y"))
      , collect("domain_label", c("x", "y"))
      , if (nrow(collect("seq_scale"))) collect("seq_scale")
    ))),
    seqArcs = primary$seq_arcs,
    axisLines = primary$axis_lines,
    axisTicks = primary$axis_ticks,
    gene_polys = primary$gene_polys,
    gene_arrows = primary$gene_labels,
    seq_labels = primary$seq_labels_df,
    show_axis = primary$show_axis
  )

  if (identical(primary$sequence_unit, "aa") && is.data.frame(primary$sequence_layout)) {
    primary$sequence_layout$units_per_bp <- NA_real_
    if ("gap_bp" %in% names(primary$sequence_layout)) {
      primary$sequence_layout$gap_positions <- primary$sequence_layout$gap_bp
      primary$sequence_layout$gap_bp <- NA_real_
    }
  }
  if (!is.null(plot$ggchord$ref)) {
    if (isTRUE(plot$coordinates$ggchord_radial)) primary <- ggchord_project_radial_layout(primary, plot) else
    if (isTRUE(plot$coordinates$ggchord_sankey)) primary <- ggchord_project_sankey_layout(primary,plot) else
      if (isTRUE(primary$collinear) && identical(plot$coordinates$collinear_seq_axis, "vertical"))
      primary <- ggchord_project_vertical_layout(primary)
  }
  projected_labels <- any(vapply(primary$layer_geometry, function(owner) any(vapply(
    owner[intersect(names(owner), c("gene_text", "gene_text_repel", "gene_label_repel"))],
    function(d) is.data.frame(d) && nrow(d) > 0L, logical(1L))), logical(1L)))
  if (projected_labels && (isTRUE(primary$curved_backbone) ||
      isTRUE(primary$radial_projected) || isTRUE(primary$sankey_projected) ||
      isTRUE(primary$vertical_projected))) {
    primary$layer_geometry <- ggchord_collinear_share_labels(primary$layer_geometry, primary)
    primary$collinear_label_diagnostics <- attr(primary$layer_geometry, "collinear_label_diagnostics")
    attr(primary$layer_geometry, "collinear_label_diagnostics") <- NULL
    labels <- lapply(primary$layer_geometry, function(owner) {
      d <- owner$gene_label_repel
      if (is.data.frame(d) && nrow(d)) d <- d[d$.component %in% "text", , drop = FALSE]
      ggchord_rbind_fill(Filter(function(x) is.data.frame(x) && nrow(x),
        list(owner$gene_text, owner$gene_text_repel, d)))
    })
    primary$gene_labels <- ggchord_rbind_fill(labels) %||% data.frame()
    segments <- lapply(primary$layer_geometry, function(owner) {
      d <- owner$gene_label_repel
      if (is.data.frame(d) && nrow(d)) d <- d[d$.component %in% "segment", , drop = FALSE]
      ggchord_rbind_fill(Filter(function(x) is.data.frame(x) && nrow(x),
        list(owner$gene_label_segment, d)))
    })
    primary$gene_label_segments <- ggchord_rbind_fill(segments) %||% data.frame()
    # Discard pre-placement text extents; they must not enlarge the final panel.
    geometry <- c(unlist(primary$layer_geometry, recursive = FALSE),
      list(primary$axis_lines, primary$axis_ticks))
    xx <- yy <- numeric()
    for (d in geometry) if (is.data.frame(d)) for (pair in list(c("x", "y"),
        c("x0", "y0"), c("x1", "y1"), c("xend", "yend"), c("text_x", "text_y"),
        c("label_x", "label_y"))) {
      if (all(pair %in% names(d))) { xx <- c(xx, d[[pair[1]]]); yy <- c(yy, d[[pair[2]]]) }
    }
    primary$extremes <- list(x_min = min(xx, na.rm = TRUE), x_max = max(xx, na.rm = TRUE),
      y_min = min(yy, na.rm = TRUE), y_max = max(yy, na.rm = TRUE))
  }
  primary$axis_ticks <- ggchord_space_axis_endpoints(primary)
  label_ids <- names(groups)[group_type == "domain_label"]
  label_specs <- list()
  shared_labels <- NULL
  if (length(label_ids)) {
    shared_labels <- new.env(parent = emptyenv())
    shared_labels$metric <- ggchord_annotation_metric(primary)
    shared_labels$scale <- if (isTRUE(primary$circular)) c(1, 1) else shared_labels$metric * 25.4
    shared_labels$boxes <- ggchord_annotation_text_obstacles(primary$layer_geometry, primary,
      include_axis = FALSE, scale = shared_labels$scale)
    shared_labels$endpoints <- ggchord_annotation_text_obstacles(list(), primary,
      scale = shared_labels$scale, endpoints_only = TRUE)
    shared_labels$polygons <- ggchord_annotation_polygons(primary$layer_geometry, shared_labels$scale)
    shared_labels$polygons_bounds <- ggchord_polygon_bounds(shared_labels$polygons)
    shared_labels$leaders <- data.frame()
    if (!isTRUE(primary$circular)) {
      shared_labels$backbones <- ggchord_rbind_fill(
        ggchord_annotation_backbones(primary, shared_labels$scale)) %||% data.frame()
    }
    fixed <- vapply(label_ids, function(id) !inherits(
      plot$layers[[groups[[id]][1L]]]$ggchord_params$label_layout, "LabelLayoutRepel"),
      logical(1L))
    label_ids <- label_ids[order(!fixed)]
  }
  for (id in label_ids) {
    layer <- plot$layers[[groups[[id]][1L]]]
    p <- layer$ggchord_params
    owner <- p$anchor_layer
    component <- if (p$label_kind == "site") "restriction_site" else "gene_poly"
    target <- if (!is.null(owner)) inputs[[owner]][[component]] else NULL
    target_geometry <- if (!is.null(owner)) primary$layer_geometry[[owner]][[component]] else NULL
    input <- inputs[[id]]$domain_label
    if (!is.null(owner) && NROW(input) && !NROW(target))
      ggchord_stop("Label anchor '", owner, "' is absent from this panel")
    p$bold_callout <- isTRUE(primary$circular) && ggchord_domain_bold_callout(layer, plot)
    geometry <- ggchord_domain_label_geometry(layer, input, primary,
      target, target_geometry, chord$display_instances,
      radial = isTRUE(plot$coordinates$ggchord_radial),
      vertical = identical(plot$coordinates$collinear_seq_axis, "vertical"),
      shared = shared_labels, bold_callout = p$bold_callout)
    primary$layer_geometry[[id]]$domain_label <- ggchord_attach_input_columns(geometry, input)
    label_specs[[id]] <- p
  }
  if (!is.null(shared_labels)) {
    if (isTRUE(primary$circular))
      primary$layer_geometry <- ggchord_place_circular_domain_labels(primary, label_specs)
    primary$layer_geometry <- ggchord_mark_leader_crossings(primary$layer_geometry, primary)
  }
  primary$axis_heading <- ggchord_axis_heading(plot, primary)
  primary$axis_geometry <- ggchord_axis_geometry(primary)
  if (isTRUE(chord$sequence_layout$auto_row_gap) ||
      (isTRUE(plot$coordinates$ggchord_sankey) && is.null(plot$coordinates$sankey_column_gap_mm)))
    primary$reference_spacing <- ggchord_reference_spacing(primary)
  primary <- ggchord_check_displayed_hsps(primary, plot)
  primary
}


# ====================================================================
# Shared helpers used by ggplot_build.ggchord() and layout preparation.
# ====================================================================
