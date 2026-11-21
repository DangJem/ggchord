#' Highlight reference intervals
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param shape Default shape. An explicitly supplied fixed shape takes priority over its mapping.
#' @param thickness Physical thickness in mm or absolute grid units. Named values select references.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param alpha Opacity from zero to one.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Provide one-based closed intervals with `accver`, `start` and `end`. A strand or disjoint segments can be supplied where appropriate. The gene shortcut defaults to an arrow and strand placement. The primer shortcut defaults to a directional narrow path, and the region shortcut defaults to a translucent block.
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' ggchord(seqs) +
#'   geom_seq_region(data = genes, fill = "#E5B66A", thickness = 4) +
#'   geom_seq() +
#'   coord_collinear()
geom_seq_region <- function(mapping = NULL, data = NULL, shape = "block",
    thickness = 2, stat = "identity", position = position_offset(0),
    alpha = 0.25, ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  layer <- ggchord_feature_layer(mapping, data, shape, !missing(shape), thickness,
    stat, position, c(list(alpha = alpha), list(...)), na.rm, show.legend, inherit.aes, id)
  layer$ggchord_params$style_defaults["alpha"] <- missing(alpha)
  layer
}
