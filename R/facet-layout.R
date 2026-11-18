# Resolve sample membership before lanes and label collisions. ggplot2 owns
# the grid, strips and common visual scales.
ggchord_facet_vars <- function(plot) {
  facet <- plot$facet
  if (inherits(facet, "FacetNull")) return(character())
  if (!inherits(facet, c("FacetWrap", "FacetGrid")))
    ggchord_stop("Sequence layouts support facet_wrap() and facet_grid()")
  p <- facet$params
  if (any(unlist(p$free)) || !identical(p$margins %||% FALSE, FALSE))
    ggchord_stop("Sequence panels require fixed scales and no facet margins")
  vars <- c(p$facets, p$rows, p$cols)
  plain <- vapply(vars, function(x) rlang::is_symbol(rlang::get_expr(x)), logical(1))
  if (!all(plain) || !identical(unname(vapply(vars, rlang::as_name, character(1))), names(vars)))
    ggchord_stop("Create facet columns in the source tables, then use vars(column)")
  reserved <- intersect(names(vars), c("PANEL", "ROW", "COL", "group", "source_row",
    "x", "y", "xend", "yend", "layer_id", "component"))
  if (length(reserved)) ggchord_stop("Facet column conflicts with a geometry field: ",
    paste(reserved, collapse = ", "), ". Copy it to a sample column before faceting.")
  names(vars)
}

# Panel numbers come from the whole plot, including observations that a
# layer's layout later excludes. A selected track still indexes its full input.
ggchord_facet_frames <- function(plot) {
  lapply(plot$layers, function(layer) {
    type <- layer$ggchord_params$type
    if (is.null(type)) return(if (is.function(layer$data)) layer$data(plot$data) else layer$data)
    if (identical(type, "track") && !is.null(layer$ggchord_params$track_sources))
      return(layer$ggchord_params$track_sources)
    fallback <- ggchord_layer_fallback(layer, plot$ggchord$data)
    if (!is.null(layer$ggchord_params$native_geom)) {
      return(ggchord_layer_source(layer, fallback))
    }
    ggchord_resolve_layer_input(layer, fallback)
  })
}

ggchord_facet_panels <- function(plot, sources) {
  frames <- Filter(is.data.frame, c(list(plot$data), sources))
  params <- plot$facet$setup_params(frames, plot$facet$params)
  list(panels = plot$facet$compute_layout(frames, params), params = params)
}

ggchord_check_panel_layout <- function(layer, reference = FALSE) {
  layout <- layer$layout
  if (is.null(layout) || (is.numeric(layout) && all(is.finite(layout) & layout == floor(layout))))
    return(invisible(NULL))
  if (reference && is.character(layout) && length(layout) == 1L &&
      layout %in% c("fixed", "fixed_rows", "fixed_cols")) return(invisible(NULL))
  if (reference) ggchord_stop("Reference layer layout needs integer panel numbers or a fixed facet layout")
  ggchord_stop("Observation layers accept integer panel numbers in layout. ",
    "Repeating observations across samples would change their ownership; supply sample-keyed rows.")
}

# Train the default feature roles on observations admitted to any panel.
# Shape affects the polygons themselves, so use the common scale before laying
# out the panels and retain the same keys for the final guides.
ggchord_facet_feature_styles <- function(plot, sources, members) {
  fills <- explicit_colours <- shapes <- list()
  strand_fill <- FALSE
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    if (!identical(layer$ggchord_type, "gene_poly") || !NROW(sources[[i]])) next
    map <- members[[i]]
    if (is.numeric(layer$layout)) map <- map[as.integer(map$PANEL) %in% layer$layout, , drop = FALSE]
    rows <- sort(unique(map$.ggchord_input_index))
    data <- sources[[i]][rows, , drop = FALSE]
    strand_fill <- strand_fill || (nrow(data) > 0L &&
      !is.null(layer$mapping$gene_fill) && identical(layer$ggchord_params$gene_color_scheme, "strand"))
    if (isTRUE(layer$ggchord_params$feature_shape_mapped)) shapes[[i]] <- data$.feature_shape_raw
    fill <- layer$mapping$feature_fill
    if (is.null(fill) || "feature_fill" %in% names(layer$aes_params) ||
        "feature_fill" %in% names(layer$ggchord_input_mapping)) next
    values <- as.character(data[[rlang::as_name(fill)]])
    fills[[i]] <- values
    colour <- data$feature_color
    if (length(colour)) {
      valid <- !is.na(values) & !is.na(colour) & nzchar(colour)
      explicit_colours[[i]] <- stats::setNames(as.character(colour[valid]), values[valid])
    }
  }
  result <- list(scales = list())
  if (strand_fill) {
    result$gene_color_scheme <- "strand"
    result$gene_pal <- process_strand_colors(NULL)
  }
  fill_order <- unique(unlist(fills, use.names = FALSE))
  fill_order <- fill_order[!is.na(fill_order)]
  if (length(fill_order)) {
    colours <- do.call(c, explicit_colours)
    colours <- colours[!duplicated(names(colours))]
    result$feature_fill_order <- fill_order
    result$feature_fill_pal <- process_manual_colors(colours, fill_order, NULL)
    if (!plot$scales$has_scale("feature_fill")) {
      result$scales$fill <- if (ggchord_uses_plasmid_feature_scale(plot$coordinates$ggchord_circular, fill_order)) {
        scale_feature_fill_plasmid(guide = "none")
      } else scale_feature_fill_manual(values = result$feature_fill_pal, guide = "none")
    }
  }
  shape_order <- unique(unlist(shapes, use.names = FALSE))
  if (length(shape_order)) {
    shape_scale <- plot$scales$get_scales("feature_shape")
    if (is.null(shape_scale)) {
      shape_scale <- if (ggchord_uses_plasmid_feature_scale(plot$coordinates$ggchord_circular, shape_order)) {
        scale_feature_shape_plasmid()
      } else scale_feature_shape_manual(values = ggchord_default_feature_shapes(shape_order))
      result$scales$shape <- shape_scale
    }
    trained <- shape_scale$clone()
    trained$train(shape_order)
    mapped_shapes <- as.character(trained$map(shape_order))
    mapped_shapes[is.na(mapped_shapes)] <- "arrow"
    result$feature_shape_order <- shape_order
    result$feature_shape_pal <- stats::setNames(mapped_shapes, shape_order)
  }
  result
}

ggchord_facet_layout <- function(plot) {
  shared <- c("seq", "seq_label", "seq_break", "seq_scale", "axis", "seq_center_label")
  types <- vapply(plot$layers, function(layer) layer$ggchord_params$type %||% "", character(1))
  has_layout <- any(vapply(plot$layers, function(layer)
    !is.null(layer$ggchord_type) && !is.null(layer$layout), logical(1)))
  references_only <- all(!nzchar(types) | types %in% shared) && !has_layout
  if (references_only) {
    # Keep native facets (including expressions) on a shared frame. When a
    # backbone itself supplies a fixed facet key it defines a frame per panel.
    keys <- tryCatch(ggchord_facet_vars(plot), error = function(e) character())
    keyed_backbone <- any(vapply(plot$layers[types == "seq"], function(layer)
      any(keys %in% names(layer$ggchord_params$backbone_input %||%
        ggchord_resolve_layer_input(layer, plot$ggchord$data$seq_data))), logical(1)))
    if (!keyed_backbone) {
      plot$facet <- ggplot2::facet_null()
      return(compute_chord_geometry_prepared(plot))
    }
  }
  keys <- ggchord_facet_vars(plot)
  sources <- lapply(plot$layers, function(layer) {
    type <- layer$ggchord_params$type
    if (is.null(type)) {
      input <- if (is.function(layer$data)) layer$data(plot$data) else layer$data
      return(if (is.data.frame(input)) input else data.frame())
    }
    ggchord_check_panel_layout(layer, reference = type %in% shared)
    if (!type %in% c(shared, "track", "gene", "gene_label", "gene_label_repel", "link", "ribbon", "seq_region"))
      ggchord_stop("Faceting this layer requires a panel ownership rule: ", type)
    if (type == "ribbon" && !is.null(layer$ggchord_params$ribbon_stat))
      ggchord_stop("Facet source ribbons directly; for aggregated ribbons compute each panel's table and retain its source members")
    fallback <- ggchord_layer_fallback(layer, plot$ggchord$data)
    input <- ggchord_resolve_layer_input(layer, fallback)
    absent <- setdiff(keys, names(input))
    if (!type %in% shared && length(absent) && NROW(input) > 0L)
      ggchord_stop("Layer '", layer$ggchord_layer_id, "' needs facet column(s): ",
        paste(absent, collapse = ", "), ". Supply observations with their sample keys.")
    input
  })
  full <- lapply(plot$layers, function(layer) layer$ggchord_params$track_sources)
  facet_sources <- Map(function(layer, data) {
    layer$ggchord_params$instance_facet_source %||% data
  }, plot$layers, sources)
  setup <- ggchord_facet_panels(plot, c(facet_sources, full))
  params <- setup$params
  panels <- setup$panels
  members <- Map(function(data, layer) {
    if (is.null(data)) return(NULL)
    at <- layer$ggchord_params$instance_facet_rows
    rows <- if (is.null(at)) data[intersect(keys, names(data))] else
      layer$ggchord_params$instance_facet_source[at, , drop = FALSE]
    attr(rows, "layout") <- layer$layout
    rows$.ggchord_input_index <- seq_len(nrow(data))
    plot$facet$map_data(rows, panels, params)
  }, sources, plot$layers)
  plot$ggchord$shared_track_scales <- ggchord_facet_track_scales(plot, sources, members, panels)
  feature_styles <- ggchord_facet_feature_styles(plot, sources, members)
  if (length(feature_styles$scales)) {
    plot$scales <- plot$scales$clone()
    for (scale in feature_styles$scales) plot$scales$add(scale)
  }
  plot$ggchord$facet_feature_shape_order <- feature_styles$feature_shape_order
  with_panel <- function(data, panel) {
    if (is.data.frame(data)) {
      data <- ggchord_restore_instance_data(data, plot$ggchord$display_instances, fields = keys)
      for (key in keys) {
        value <- rep(panels[[key]][panel], nrow(data))
        present <- !is.na(data[[key]])
        if (any(present) && !identical(as.character(data[[key]][present]), as.character(value[present])))
          ggchord_stop("Facet column conflicts with layout output: ", key, ". Copy the facet key to a sample column.")
        data[[key]] <- value
      }
      data$PANEL <- rep.int(panel, nrow(data))
    }
    data
  }
  layouts <- lapply(seq_len(nrow(panels)), function(panel) {
    if (references_only && !any(vapply(members[types == "seq"], function(x)
        panel %in% as.integer(x$PANEL), logical(1)))) return(NULL)
    sub <- plot
    sub$facet <- ggplot2::facet_null()
    selected <- which(vapply(plot$layers, function(layer) is.null(layer$ggchord_type) ||
      !is.numeric(layer$layout) || panel %in% layer$layout, logical(1)))
    sub$layers <- lapply(selected, function(i) {
      layer <- plot$layers[[i]]
      if (is.null(layer$ggchord_type)) return(layer)
      copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
      map <- members[[i]]
      rows <- map$.ggchord_input_index[as.integer(map$PANEL) == panel]
      panel_input <- sources[[i]][rows, , drop = FALSE]
      source <- attr(panel_input, "ggchord_source")
      if (!is.null(source)) {
        source$rows <- source$rows[rows]
        attr(panel_input, "ggchord_source") <- source
      }
      copy$ggchord_params$panel_input <- panel_input
      copy
    })
    backbone <- selected[types[selected] == "seq"]
    arrangement <- sub$ggchord$sequence_layout
    if (length(backbone) && !is.null(arrangement)) {
      i <- tail(backbone, 1L)
      all_ids <- as.character(sources[[i]]$accver)
      original <- as.character(plot$ggchord$data$seq_data$accver)
      ordered <- arrangement$order %||% c(intersect(original, all_ids), setdiff(all_ids, original))
      if (!all(ordered %in% all_ids)) ggchord_stop(
        "layout_sequences(): order must contain known sequence IDs from the displayed reference")
      reference <- sub$layers[[match(i, selected)]]$ggchord_params$panel_input
      local <- intersect(ordered, reference$accver)
      if (!length(local)) return(NULL)
      if (length(local) < length(ordered)) {
        arrangement$order <- local
        for (field in c("orientation", "gap")) if (!is.null(arrangement[[field]])) {
          default <- if (field == "orientation") 1 else
            if (isTRUE(plot$coordinates$ggchord_collinear)) 0.3 else 0.03
          arrangement[[field]] <- process_sequence_param(arrangement[[field]], ordered,
            paste0("layout_sequences(", field, ")"), default)[local]
        }
        sub$ggchord$sequence_layout <- arrangement
      }
    }
    metrics <- ggchord_render_context()$facet_contexts
    layout <- if (is.null(metrics)) compute_chord_geometry_prepared(sub) else
      ggchord_with_render_context(metrics[[panel]], compute_chord_geometry_prepared(sub))
    for (i in seq_along(sub$layers)) {
      layer <- sub$layers[[i]]
      id <- layer$ggchord_layer_id
      if (is.null(id)) next
      map <- members[[selected[i]]]
      rows <- map$.ggchord_input_index[as.integer(map$PANEL) == panel]
      component <- layer$ggchord_type
      geom <- layout$layer_geometry[[id]][[component]]
      if (is.data.frame(geom)) {
        if ("source_row" %in% names(geom)) geom$source_row <- rows[geom$source_row]
        layout$layer_geometry[[id]][[component]] <- with_panel(geom, panel)
      }
      for (field in c("connection_report", "track_diagnostics", "collinear_label_diagnostics", "circular_annotation_registry")) {
        table <- layout[[field]]
        if (!is.data.frame(table) || !"source_row" %in% names(table)) next
        owner <- table$layer_id %||% table$owner
        idx <- which(owner == id)
        if (length(idx)) {
          table$source_row[idx] <- rows[table$source_row[idx]]
          layout[[field]] <- table
        }
      }
    }
    # Circular label clearance and projected axes need this panel's own reference.
    layout$axis_geometry <- ggchord_axis_geometry(layout)
    layout$sequence_display <- ggchord_sequence_display(layout)
    layout
  })
  nonempty <- Filter(Negate(is.null), layouts)
  result <- nonempty[[1L]]
  style_fields <- setdiff(names(feature_styles), "scales")
  result[style_fields] <- feature_styles[style_fields]
  for (layout in layouts) {
    extra <- setdiff(names(layout$track_scales), names(result$track_scales))
    if (length(extra)) result$track_scales[extra] <- layout$track_scales[extra]
  }
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    id <- layer$ggchord_layer_id
    if (is.null(id)) next
    type <- layer$ggchord_type
    result$layer_inputs[[id]][[type]] <- sources[[i]]
    result$layer_geometry[[id]][[type]] <- ggchord_rbind_fill(lapply(layouts, function(x)
      x$layer_geometry[[id]][[type]])) %||% data.frame()
    if (identical(type, "track")) {
      result$track_statistics[[id]] <- ggchord_rbind_fill(lapply(layouts, function(x)
        x$track_statistics[[id]])) %||% sources[[i]][FALSE, , drop = FALSE]
      result$track_source_data[[id]] <- layer$ggchord_params$track_sources %||% sources[[i]]
    }
  }
  components <- c(ribbon_polys = "ribbon", link_lines = "link", region_polys = "seq_region")
  for (field in names(components)) {
    tables <- lapply(result$layer_geometry, `[[`, components[[field]])
    result[[field]] <- ggchord_rbind_fill(Filter(function(x) NROW(x) > 0L, tables))
  }
  if (NROW(result$ribbon_polys) && !NROW(layouts[[1L]]$ribbon_polys)) {
    drawn <- Find(function(x) NROW(x$ribbon_polys) > 0L, layouts)
    fields <- intersect(c("ribbon_color_scheme", "ribbon_colors", "ribbon_color_name",
      "ribbon_color_limits", "ribbon_color_breaks"), names(drawn))
    result[fields] <- drawn[fields]
  }
  if (isTRUE(result$collinear)) {
    fields <- c("seqs", "sequence_layout", "sequence_reference", "sequence_unit",
      "seq_axis", "display_instances", "sankey_projected", "sankey_frames", "radial_frames")
    result$panel_references <- lapply(layouts, function(x) x[intersect(fields, names(x))])
  }
  for (field in c("axis_geometry", "axis_lines", "axis_ticks", "axis_selection",
      "sequence_display", "sequence_layout", "orientation_diagnostics", "radial_frames",
      "radial_endpoints", "gene_polys", "gene_labels", "gene_label_segments", "track_registry",
      "connection_report", "track_axes", "track_diagnostics", "circular_annotation_registry", "collinear_label_diagnostics")) {
    tables <- lapply(seq_along(layouts), function(i) with_panel(layouts[[i]][[field]], i))
    if (any(vapply(tables, is.data.frame, logical(1)))) result[[field]] <- ggchord_rbind_fill(tables)
  }
  for (axis in c("x", "y")) for (bound in c("min", "max")) {
    field <- paste(axis, bound, sep = "_")
    values <- vapply(nonempty, function(x) x$extremes[[field]], numeric(1))
    result$extremes[[field]] <- if (bound == "min") min(values) else max(values)
  }
  result$facet_layout <- panels
  if (has_layout) result$facet_keys <- keys
  ggchord_check_displayed_hsps(result, plot)
}
