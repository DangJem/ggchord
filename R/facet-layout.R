# Resolve sample membership before lanes and label collisions. ggplot2 owns
# the grid, strips and common visual scales.
ggchord_facet_vars <- function(plot) {
  facet <- plot$facet
  if (inherits(facet, "FacetNull")) return(character())
  if (!inherits(facet, c("FacetWrap", "FacetGrid")))
    ggchord_stop("Sequence layouts support facet_wrap() and facet_grid()")
  p <- facet$params
  ggchord_validate_coord_contract(plot)
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
    if (is.null(type)) return(layer$layer_data(plot$data))
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
  frames <- Filter(is.data.frame, c(list(plot$ggchord$facet_reference %||%
    plot$ggchord$data$seq_data), sources))
  params <- plot$facet$setup_params(frames, plot$facet$params)
  list(panels = plot$facet$compute_layout(frames, params), params = params)
}

# Determine panel membership from source data before statistics create rows.
ggchord_prepare_facet_sources <- function(plot) {
  if (!is.null(plot$ggchord$source_facet_setup) || inherits(plot$facet, "FacetNull"))
    return(plot)
  frames <- lapply(plot$layers, function(layer) {
    if (is.null(layer$ggchord_type)) return(layer$layer_data(plot$data))
    if (!is.null(layer$ggchord_params$track_sources)) return(layer$ggchord_params$track_sources)
    if (!is.null(layer$ggchord_params$anchor_layer)) return(layer$ggchord_input_data)
    ggchord_layer_source(layer, ggchord_layer_fallback(layer, plot$ggchord$data))
  })
  plot$ggchord$source_facet_setup <- ggchord_facet_panels(plot, frames)
  plot$ggchord$source_facet_setup$facet <- plot$facet
  plot$ggchord$source_facet_setup$keys <- ggchord_facet_vars(plot)
  plot
}

ggchord_facet_stat_input <- function(data, plot, layer) {
  ggchord_source_stat_members(data, plot$ggchord$source_facet_setup, layer)
}

ggchord_source_stat_members <- function(data, setup, layer) {
  if (is.null(setup)) return(NULL)
  rows <- data[intersect(setup$keys, names(data))]
  rows$.ggchord_input_index <- seq_len(nrow(data))
  attr(rows, "layout") <- layer$layout
  setup$facet$map_data(rows, setup$panels, setup$params)
}

# Draw against the same panel table used to select and project source rows.
# The original Facet still maps observations, assigns scales and draws strips.
ggchord_precomputed_facet <- function(facet, panels) {
  force(panels)
  ggplot2::ggproto(NULL, facet,
    compute_layout = function(data, params) panels,
    map_data = function(data, layout, params) {
      if (!".ggchord_panel" %in% names(data)) return(facet$map_data(data, layout, params))
      data$PANEL <- layout$PANEL[match(data$.ggchord_panel, as.integer(layout$PANEL))]
      data$.ggchord_panel <- NULL
      data
    })
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

ggchord_facet_backbones <- function(plot, panels, params) {
  data <- plot$ggchord$facet_reference %||% plot$ggchord$data$seq_data
  ids <- plot$ggchord$data$seq_data$accver
  keys <- ggchord_facet_vars(plot)
  rows <- data[intersect(keys, names(data))]
  rows$.ggchord_input_index <- seq_len(nrow(data))
  membership <- plot$facet$map_data(rows, panels, params)
  order <- plot$ggchord$sequence_layout$order %||% ids
  lapply(seq_len(nrow(panels)), function(panel) {
    local <- membership$.ggchord_input_index[as.integer(membership$PANEL) == panel]
    list(order = order, displayed = intersect(order, ids[local]))
  })
}

# Train feature styles only on rows whose reference is displayed in that panel.
ggchord_facet_feature_styles <- function(plot, sources, members, backbones) {
  fills <- explicit_colours <- shapes <- list()
  strand_fill <- FALSE
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    if (!identical(layer$ggchord_type, "gene_poly") || !NROW(sources[[i]])) next
    map <- members[[i]]
    if (is.numeric(layer$layout)) map <- map[as.integer(map$PANEL) %in% layer$layout, , drop = FALSE]
    rows <- split(map$.ggchord_input_index, as.integer(map$PANEL))
    rows <- Map(function(idx, panel) idx[sources[[i]]$accver[idx] %in% backbones[[panel]]$displayed],
      rows, as.integer(names(rows)))
    rows <- sort(unique(as.integer(unlist(rows, use.names = FALSE))))
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
  shared <- c("seq", "seq_label", "seq_break", "seq_scale", "axis", "seq_center_label", "domain_label")
  types <- vapply(plot$layers, function(layer) layer$ggchord_params$type %||% "", character(1))
  keys <- ggchord_facet_vars(plot)
  sources <- lapply(plot$layers, function(layer) {
    type <- layer$ggchord_params$type
    if (is.null(type)) {
      input <- layer$layer_data(plot$data)
      return(if (is.data.frame(input)) input else data.frame())
    }
    ggchord_check_panel_layout(layer, reference = type %in% shared)
    if (!type %in% c(shared, "track", "gene", "gene_label", "gene_label_repel", "link", "ribbon", "seq_region", "restriction_site"))
      ggchord_stop("Faceting this layer requires a panel ownership rule: ", type)
    fallback <- ggchord_layer_fallback(layer, plot$ggchord$data)
    ggchord_resolve_layer_input(layer, fallback)
  })
  full <- lapply(plot$layers, function(layer) layer$ggchord_params$track_sources)
  facet_sources <- Map(function(layer, data) {
    layer$ggchord_params$instance_facet_source %||% data
  }, plot$layers, sources)
  setup <- plot$ggchord$source_facet_setup %||% ggchord_facet_panels(plot, c(facet_sources, full))
  params <- setup$params
  panels <- setup$panels
  members <- Map(function(data, layer) {
    if (is.null(data)) return(NULL)
    if (".ggchord_stat_panel" %in% names(data)) {
      return(data.frame(.ggchord_input_index = seq_len(nrow(data)),
        PANEL = panels$PANEL[match(data$.ggchord_stat_panel, as.integer(panels$PANEL))]))
    }
    at <- layer$ggchord_params$instance_facet_rows
    rows <- if (is.null(at)) data[intersect(keys, names(data))] else
      layer$ggchord_params$instance_facet_source[at, , drop = FALSE]
    attr(rows, "layout") <- layer$layout
    rows$.ggchord_input_index <- seq_len(nrow(data))
    plot$facet$map_data(rows, panels, params)
  }, sources, plot$layers)
  plot$ggchord$shared_track_scales <- ggchord_facet_track_scales(plot, sources, members, panels)
  backbones <- ggchord_facet_backbones(plot, panels, params)
  if (isTRUE(plot$coordinates$ggchord_circular)) {
    reference <- plot$ggchord$data$seq_data
    source <- plot$ggchord$display_instances
    originals <- if (is.null(source)) stats::setNames(reference$accver, reference$accver) else
      stats::setNames(source$accver, source$key)
    bad <- which(vapply(backbones, function(x) length(unique(originals[x$displayed])) != 1L, logical(1)))
    if (length(bad)) ggchord_stop("coord_circular() requires one reference per panel (panel ", bad[1L], ")")
  }
  feature_styles <- ggchord_facet_feature_styles(plot, sources, members, backbones)
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
        data[[key]] <- value
      }
      data$PANEL <- data$.ggchord_panel <- rep.int(panel, nrow(data))
    }
    data
  }
  layouts <- lapply(seq_len(nrow(panels)), function(panel) {
    local <- backbones[[panel]]$displayed
    if (!length(local)) return(NULL)
    sub <- plot
    sub$facet <- ggplot2::facet_null()
    sub$ggchord$data$seq_data <- sub$ggchord$data$seq_data[
      match(local, sub$ggchord$data$seq_data$accver), , drop = FALSE]
    if (!is.null(sub$ggchord$display_instances))
      sub$ggchord$display_instances <- sub$ggchord$display_instances[
        sub$ggchord$display_instances$key %in% local, , drop = FALSE]
    selected <- which(vapply(plot$layers, function(layer) is.null(layer$ggchord_type) ||
      !is.numeric(layer$layout) || panel %in% layer$layout, logical(1)))
    panel_rows <- vector("list", length(plot$layers))
    sub$layers <- lapply(selected, function(i) {
      layer <- plot$layers[[i]]
      if (is.null(layer$ggchord_type)) return(layer)
      copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
      map <- members[[i]]
      rows <- map$.ggchord_input_index[as.integer(map$PANEL) == panel]
      source_data <- sources[[i]]
      if ("accver" %in% names(source_data)) rows <- rows[source_data$accver[rows] %in% local]
      if (all(c("qaccver", "saccver") %in% names(source_data)))
        rows <- rows[source_data$qaccver[rows] %in% local & source_data$saccver[rows] %in% local]
      panel_rows[[i]] <<- rows
      copy$ggchord_params$visual_source <- ggchord_visual_rows(layer$ggchord_params$visual_source, rows)
      panel_input <- sources[[i]][rows, , drop = FALSE]
      source <- attr(panel_input, "ggchord_source")
      if (!is.null(source)) {
        source$rows <- source$rows[rows]
        attr(panel_input, "ggchord_source") <- source
      }
      copy$ggchord_params$panel_input <- panel_input
      if (types[i] %in% c("link", "ribbon")) {
        field <- paste0(types[i], "_ctrl_point")
        copy$ggchord_params[[field]] <- ggchord_link_control_rows(layer$ggchord_params[[field]], rows)
      }
      copy
    })
    arrangement <- sub$ggchord$sequence_layout
    ordered <- backbones[[panel]]$order
    if (!is.null(arrangement) && length(local) < length(ordered)) {
      arrangement$order <- local
      for (field in c("orientation", "gap", "row_gap", "radius", "curvature")) {
        if (!is.null(arrangement[[field]])) arrangement[[field]] <- arrangement[[field]][local]
      }
      sub$ggchord$sequence_layout <- arrangement
    }
    metrics <- ggchord_render_context()$facet_contexts
    layout <- if (is.null(metrics)) compute_chord_geometry_prepared(sub) else
      ggchord_with_render_context(metrics[[panel]], compute_chord_geometry_prepared(sub))
    for (i in seq_along(sub$layers)) {
      layer <- sub$layers[[i]]
      id <- layer$ggchord_layer_id
      if (is.null(id)) next
      rows <- panel_rows[[selected[i]]]
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
  result <- if (length(nonempty)) nonempty[[1L]] else {
    # Layer and sequence selections can exclude every backbone. Keep the
    # registered inputs and native panels without inventing a reference frame.
    coord <- plot$coordinates
    reference <- sources[[tail(which(types == "seq"), 1L)]]
    list(collinear = isTRUE(coord$ggchord_collinear), circular = isTRUE(coord$ggchord_circular),
      seq_axis = coord$collinear_seq_axis, seqs = character(), sequence_layout = data.frame(),
      sequence_display = data.frame(), sequence_unit = ggchord_sequence_unit(reference),
      length_mode = plot$ggchord$sequence_layout$length_mode, radial_frames = data.frame(key = character()),
      display_instances = plot$ggchord$display_instances, instance_routing = plot$ggchord$instance_routing,
      instance_omissions = plot$ggchord$instance_omissions)
  }
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
  fields <- c("seqs", "sequence_layout", "sequence_reference", "sequence_unit",
    "collinear", "circular", "rotation", "seq_axis", "display_instances",
    "sankey_projected", "sankey_frames", "radial_projected", "radial_frames")
  result$panel_references <- lapply(layouts, function(x) {
    if (is.null(x)) return(list(seqs = character()))
    x[intersect(fields, names(x))]
  })
  for (field in c("axis_geometry", "axis_lines", "axis_ticks", "axis_selection",
      "sequence_display", "sequence_layout", "orientation_diagnostics", "radial_frames",
      "radial_endpoints", "gene_polys", "gene_labels", "gene_label_segments", "track_registry",
      "connection_report", "track_axes", "track_diagnostics", "circular_annotation_registry", "collinear_label_diagnostics")) {
    tables <- lapply(seq_along(layouts), function(i) with_panel(layouts[[i]][[field]], i))
    if (any(vapply(tables, is.data.frame, logical(1)))) result[[field]] <- ggchord_rbind_fill(tables)
  }
  if (length(nonempty)) for (axis in c("x", "y")) for (bound in c("min", "max")) {
    field <- paste(axis, bound, sep = "_")
    values <- vapply(nonempty, function(x) x$extremes[[field]], numeric(1))
    result$extremes[[field]] <- if (bound == "min") min(values) else max(values)
  }
  result$facet_layout <- panels
  result$facet_keys <- keys
  ggchord_check_displayed_hsps(result, plot)
}
