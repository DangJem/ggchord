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
              seq_region = integer(0), restriction_site = integer(0),
              seq_center_label = integer(0))
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    type <- lyr$ggchord_type %||% ""
    if (type %in% names(idx)) idx[[type]] <- c(idx[[type]], i)
  }
  idx
}

# Register a plot-local copy, including repeated uses of the same layer object.
# IDs describe component order in this plot, never mutate the user's ggproto.
ggchord_register_layers <- function(plot) {
  next_id <- 0L
  plot$layers <- lapply(plot$layers, function(layer) {
    copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    if (!is.null(copy$ggchord_type)) {
      next_id <<- next_id + 1L
      copy$ggchord_layer_id <- sprintf("layer-%04d", next_id)
    }
    copy
  })
  plot
}
