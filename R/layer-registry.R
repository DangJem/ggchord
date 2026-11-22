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

#' Choose a layer identifier and default scale scope
#'
#' In a ggchord plot, layers that map the same visual aesthetic train their
#' own scales by default. A character ID gives a layer a stable name for
#' [scale_for()], label binding and export. An ordinary explicit scale applies
#' to all layers mapping its aesthetic except those selected by `scale_for()`.
#' Use `id = NULL` or a scalar `NA` on each layer that should share the
#' ordinary ggplot2 scale. Fixed styles and position scales are independent
#' of this choice. Native ggplot2 plots use their ordinary scales.
#' @return A declaration resolved when the layers of a plot are registered.
#' @export
#' @examples
#' library(ggplot2)
#' seqs <- data.frame(accver = c("A", "B"), length = c(1000, 800))
#' genes <- data.frame(accver = c("A", "B"), start = c(100, 150),
#'   end = c(240, 380), strand = c("+", "-"), label = c("a", "b"))
#' ggchord(seqs) +
#'   geom_seq(aes(colour = accver)) +
#'   geom_gene(aes(colour = strand), data = genes)
#' ggchord(seqs) +
#'   geom_seq(aes(colour = accver), id = NULL) +
#'   geom_gene(aes(colour = strand), data = genes, id = NA)
auto_id <- function() structure(list(), class = "ggchord_auto_id")

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

# Register a plot-local copy without changing the caller's layer object.
ggchord_register_layers <- function(plot) {
  next_id <- 0L
  declared <- lapply(plot$layers, function(layer) {
    if (is.null(layer$ggchord_type)) layer$ggchord_layer_id else ggchord_declared_id(layer$ggchord_params$id)
  })
  declared <- Filter(Negate(is.null), declared)
  for (id in declared) ggchord_check_layer_id(id)
  reserved <- unlist(declared, use.names = FALSE)
  if (anyDuplicated(reserved)) ggchord_stop("Duplicate layer id '", reserved[anyDuplicated(reserved)], "'")
  plot$layers <- lapply(plot$layers, function(layer) {
    if (is.null(layer$ggchord_type)) return(layer)
    copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    next_id <<- next_id + 1L
    id <- ggchord_declared_id(copy$ggchord_params$id)
    if (is.null(id)) {
      id <- sprintf("layer-%04d", next_id)
      while (id %in% reserved) {
        next_id <<- next_id + 1L
        id <- sprintf("layer-%04d", next_id)
      }
    }
    copy$ggchord_layer_id <- id
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
