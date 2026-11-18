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
    layout = lyr$layout, show.legend = lyr$show.legend,
    check.aes = FALSE
  )
  new$name <- lyr$name
  # Preserve the ggchord custom fields on the reconstructed layer
  for (fld in c(
    "ggchord_type", "ggchord_params", "ggchord_placeholder",
    "ggchord_layer_id", "ggchord_input_data", "ggchord_input_mapping", "ggchord_annotation_origin",
    "ggchord_role_aes", "ggchord_resolved_input", "ggchord_theme_element",
    "ggchord_theme_components", "ggchord_input_transform", "ggchord_obstacle_provider"
  )) {
    if (!is.null(lyr[[fld]])) new[[fld]] <- lyr[[fld]]
  }
  new
}

# Register a plot-local copy, including repeated uses of the same layer object.
# IDs describe component order in this plot, never mutate the user's ggproto.
ggchord_register_layers <- function(plot) {
  next_id <- 0L
  declared <- lapply(plot$layers, function(layer) layer$ggchord_params$id)
  declared <- Filter(Negate(is.null), declared)
  if (any(!vapply(declared, function(id) is.character(id) && length(id) == 1L &&
      !is.na(id) && nzchar(trimws(id)), logical(1))))
    ggchord_stop("Layer id must be a non-empty string")
  reserved <- unlist(declared, use.names = FALSE)
  if (anyDuplicated(reserved)) ggchord_stop("Duplicate layer id '", reserved[anyDuplicated(reserved)], "'")
  plot$layers <- lapply(plot$layers, function(layer) {
    if (is.null(layer$ggchord_type)) return(layer)
    copy <- if (ggchord_uses_inherited_line_xy(layer, plot$mapping)) {
      ggchord_native_track_line(layer)
    } else reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    next_id <<- next_id + 1L
    id <- copy$ggchord_params$id
    if (is.null(id)) {
      id <- sprintf("layer-%04d", next_id)
      while (id %in% reserved) {
        next_id <<- next_id + 1L
        id <- sprintf("layer-%04d", next_id)
      }
    }
    copy$ggchord_layer_id <- id
    type <- copy$ggchord_params$type
    kind <- if (type %in% c("ribbon", "alignment")) "links" else
      if (type %in% c("gene", "gene_label", "gene_label_repel") ||
          identical(copy$ggchord_params$native_stat, "coverage")) "features"
    if (is.null(copy$ggchord_params$input_description)) {
      input <- copy$ggchord_input_data
      copy$ggchord_params$input_description <- if (is.function(input)) "layer data callback" else
        if (!is.null(input)) paste0("layer data (", class(input)[1L], ")") else
        if (length(kind)) paste(plot$ggchord$candidates[[kind]], collapse = ", ") else
        if (type == "seq_region") "region data" else "reference data"
    }
    if (length(kind) && (is.null(copy$ggchord_input_data) || is.function(copy$ggchord_input_data))) {
      candidates <- plot$ggchord$candidates[[kind]]
      if (length(candidates) > 1L)
        ggchord_stop("Layer '", id, "' has multiple eligible inputs: ", paste(candidates, collapse = ", "),
          ". Select one with data = prepared$", candidates[length(candidates)], ".")
      input <- if (kind == "links") plot$ggchord$data$ribbon_data else
        plot$ggchord$data$feature_data %||% plot$ggchord$data$gene_data
      if (is.null(input) && is.null(copy$ggchord_input_data))
        ggchord_stop("Layer '", id, "' requires ", kind, ". Supply a layer data table or prepare_ggchord_data().")
    }
    if (isTRUE(copy$inherit.aes) && length(plot$mapping) && !identical(type, "track")) {
      copy <- ggchord_inherit_ribbon_linetype(copy, plot$mapping)
      local <- copy$ggchord_input_mapping %||% ggplot2::aes()
      fixed <- names(copy$aes_params)
      defaults <- copy$ggchord_params$style_defaults
      fixed <- setdiff(fixed, names(defaults)[defaults])
      if ("seq_colour" %in% fixed) fixed <- c(fixed, "colour")
      if (any(c("ribbon_fill", "feature_fill", "gene_fill") %in% fixed)) fixed <- c(fixed, "fill")
      for (nm in setdiff(names(plot$mapping), c(names(local), fixed))) local[[nm]] <- plot$mapping[[nm]]
      copy$ggchord_input_mapping <- local
      copy$inherit.aes <- FALSE
    }
    copy
  })
  plot
}
