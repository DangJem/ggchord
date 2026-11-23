#' Reconstruct a layer with the given data (and optional remapped mapping).
#'
#' LayerInstance objects cannot be cloned with \code{ggproto(NULL, .)}, so the
#' layer is rebuilt through \code{layer()} with the same geom/stat/mapping/params.
#' @noRd
reconstruct_layer <- function(lyr, data, mapping = NULL) {
  params <- c(lyr$geom_params, lyr$stat_params, lyr$aes_params)
  params <- params[!duplicated(names(params))]
  new <- layer(
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

#' Assign a plot-local layer identifier automatically
#'
#' IDs identify layers for [scale_for()], label binding and export. Ordinary
#' visual scales are shared across layers by aesthetic.
#' @return A declaration resolved when the layer is added.
#' @export
#' @examples
#' ggchord(seq_data_example) + geom_seq(id = auto_id())
auto_id <- function() {
  ggchord_local_console_errors()
  structure(list(), class = "ggchord_auto_id")
}

ggchord_shared_id <- function(id) {
  is.null(id) || (is.atomic(id) && !is.object(id) && is.null(dim(id)) &&
    length(id) == 1L && is.na(id) && !isTRUE(is.nan(id)))
}

ggchord_declared_id <- function(id) {
  if (inherits(id, "ggchord_auto_id") || ggchord_shared_id(id)) return(NULL)
  ggchord_check_layer_id(id)
  id
}

ggchord_check_layer_id <- function(id) {
  if (!is.character(id) || length(id) != 1L || is.na(id) || !nzchar(trimws(id)))
    ggchord_stop("Layer id must be a non-empty string")
  invisible(id)
}

# Freeze identity on a plot-local copy. Native layers retain their methods.
ggchord_assign_layer_id <- function(layer, occupied, index) {
  id <- layer$ggchord_layer_id %||% ggchord_declared_id(layer$ggchord_params$id)
  if (!is.null(id) && id %in% occupied) ggchord_stop("Duplicate layer id '", id, "'")
  if (is.null(id)) {
    id <- sprintf("layer-%04d", index)
    while (id %in% occupied) {
      index <- index + 1L
      id <- sprintf("layer-%04d", index)
    }
  }
  copy <- if (is.null(layer$ggchord_type)) ggproto(NULL, layer) else
    reconstruct_layer(layer, layer$data, mapping = layer$mapping)
  copy$ggchord_layer_id <- id
  copy
}

# Prepare disposable drawing copies. Older plot objects can lack frozen IDs.
ggchord_register_layers <- function(plot) {
  ids <- lapply(plot$layers, function(layer)
    layer$ggchord_layer_id %||% ggchord_declared_id(layer$ggchord_params$id))
  occupied <- unlist(ids, use.names = FALSE)
  if (anyDuplicated(occupied))
    ggchord_stop("Duplicate layer id '", occupied[anyDuplicated(occupied)], "'")
  plot$layers <- lapply(seq_along(plot$layers), function(i) {
    layer <- plot$layers[[i]]
    if (is.null(layer$ggchord_layer_id)) {
      layer <- ggchord_assign_layer_id(layer, setdiff(occupied, ids[[i]]), i)
      occupied <<- c(occupied, layer$ggchord_layer_id)
    }
    if (is.null(layer$ggchord_type)) return(layer)
    copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    type <- copy$ggchord_params$type
    copy$ggchord_params$plot_data <- list(if (!is_waiver(plot$data)) plot$data)
    copy$ggchord_params$reference <- plot$ggchord$data$seq_data
    input <- copy$ggchord_input_data
    copy$ggchord_params$input_description <- if (is.function(input)) "layer data callback" else
      if (is.null(input)) "default layer data" else "explicit layer data"

    if (!identical(type, "track")) {
      local <- copy$ggchord_input_mapping %||% aes()
      if (isTRUE(copy$inherit.aes)) {
        for (name in setdiff(names(plot$mapping), names(local))) local[[name]] <- plot$mapping[[name]]
      }
      defaults <- copy$ggchord_params$style_defaults
      fixed <- setdiff(names(copy$aes_params), names(defaults)[defaults])
      local[intersect(names(local), fixed)] <- NULL
      copy$ggchord_input_mapping <- local
      copy$inherit.aes <- FALSE
    }
    ggchord_inherit_feature_shape(copy)
  })
  ggchord_resolve_auto_anchors(plot)
}
