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

ggchord_facet_layout <- function(plot) {
  shared <- c("seq", "seq_label", "seq_break", "seq_scale", "axis", "seq_center_label")
  types <- vapply(plot$layers, function(layer) layer$ggchord_params$type %||% "", character(1))
  if (all(!nzchar(types) | types %in% shared)) {
    # Reference geometry is already keyed by its input. Let ggplot2 distribute
    # it with native layers, including data functions and facet expressions.
    plot$facet <- ggplot2::facet_null()
    return(compute_chord_geometry_prepared(plot))
  }
  keys <- ggchord_facet_vars(plot)
  sources <- lapply(plot$layers, function(layer) {
    type <- layer$ggchord_params$type
    if (is.null(type)) {
      input <- if (is.function(layer$data)) layer$data(plot$data) else layer$data
      return(if (is.data.frame(input)) input else data.frame())
    }
    fallback <- if (type %in% shared || type == "track") plot$ggchord$data$seq_data else
      if (type %in% c("gene", "gene_label", "gene_label_repel"))
        ggchord_annotation_fallback(plot$ggchord$data) else
          ggchord_stop("Faceting this layer requires a panel ownership rule: ", type)
    input <- ggchord_resolve_layer_input(layer, fallback)
    absent <- setdiff(keys, names(input))
    if (!type %in% shared && length(absent) && NROW(input) > 0L)
      ggchord_stop("Layer '", layer$ggchord_layer_id, "' needs facet column(s): ",
        paste(absent, collapse = ", "), ". Supply observations with their sample keys.")
    input
  })
  frames <- Filter(is.data.frame, c(list(plot$data), sources))
  params <- plot$facet$setup_params(frames, plot$facet$params)
  panels <- plot$facet$compute_layout(frames, params)
  members <- lapply(sources, function(data) {
    if (is.null(data)) return(NULL)
    rows <- data[intersect(keys, names(data))]
    rows$.ggchord_input_index <- seq_len(nrow(data))
    plot$facet$map_data(rows, panels, params)
  })
  with_panel <- function(data, panel) {
    if (is.data.frame(data)) for (key in keys) data[[key]] <- rep(panels[[key]][panel], nrow(data))
    data
  }
  layouts <- lapply(seq_len(nrow(panels)), function(panel) {
    sub <- plot
    sub$facet <- ggplot2::facet_null()
    sub$layers <- lapply(seq_along(plot$layers), function(i) {
      layer <- plot$layers[[i]]
      if (is.null(layer$ggchord_type)) return(layer)
      copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
      map <- members[[i]]
      rows <- map$.ggchord_input_index[as.integer(map$PANEL) == panel]
      copy$ggchord_params$panel_input <- sources[[i]][rows, , drop = FALSE]
      if (identical(layer$ggchord_type, "track")) copy$ggchord_params$track_domain_data <- sources[[i]]
      copy
    })
    layout <- compute_chord_geometry_prepared(sub)
    for (i in seq_along(sub$layers)) {
      layer <- sub$layers[[i]]
      id <- layer$ggchord_layer_id
      if (is.null(id)) next
      map <- members[[i]]
      rows <- map$.ggchord_input_index[as.integer(map$PANEL) == panel]
      component <- layer$ggchord_type
      geom <- layout$layer_geometry[[id]][[component]]
      if (!is.data.frame(geom)) next
      if ("source_row" %in% names(geom)) geom$source_row <- rows[geom$source_row]
      layout$layer_geometry[[id]][[component]] <- with_panel(geom, panel)
      for (field in c("track_diagnostics", "collinear_label_diagnostics", "circular_annotation_registry")) {
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
    layout
  })
  result <- layouts[[1L]]
  for (i in seq_along(plot$layers)) {
    layer <- plot$layers[[i]]
    id <- layer$ggchord_layer_id
    if (is.null(id)) next
    type <- layer$ggchord_type
    result$layer_inputs[[id]][[type]] <- sources[[i]]
    result$layer_geometry[[id]][[type]] <- ggchord_rbind_fill(lapply(layouts, function(x)
      x$layer_geometry[[id]][[type]]))
    if (identical(type, "track")) {
      result$track_statistics[[id]] <- ggchord_rbind_fill(lapply(layouts, function(x) x$track_statistics[[id]]))
      result$track_source_data[[id]] <- layer$ggchord_params$track_sources %||% sources[[i]]
    }
  }
  for (field in c("gene_polys", "gene_labels", "gene_label_segments", "track_registry",
      "track_axes", "track_diagnostics", "circular_annotation_registry", "collinear_label_diagnostics")) {
    tables <- lapply(seq_along(layouts), function(i) with_panel(layouts[[i]][[field]], i))
    if (any(vapply(tables, is.data.frame, logical(1)))) result[[field]] <- ggchord_rbind_fill(tables)
  }
  for (axis in c("x", "y")) for (bound in c("min", "max")) {
    field <- paste(axis, bound, sep = "_")
    values <- vapply(layouts, function(x) x$extremes[[field]], numeric(1))
    result$extremes[[field]] <- if (bound == "min") min(values) else max(values)
  }
  result$facet_layout <- panels
  result
}
