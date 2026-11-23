#' Draw primer binding intervals
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param shape Default shape. An explicitly supplied fixed shape takes priority over its mapping.
#' @param thickness Physical thickness in mm or absolute grid units. Named values select references.
#' @param stat A compatible ggplot2 Stat or its name. See [geom_feature()] for
#'   source grouping and reuse during rendering.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer for addition to a plot.
#' @details Provide one-based closed intervals with `accver`, `start` and `end`. A strand or disjoint segments can be supplied where appropriate. The gene shortcut defaults to an arrow and strand placement. The primer shortcut defaults to a directional narrow path, and the region shortcut defaults to a translucent block.
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' primers <- data.frame(
#'   accver = "A", feature_id = c("p1", "p2"),
#'   start = c(80, 700), end = c(100, 720),
#'   strand = c("+", "-"), label = c("Forward", "Reverse")
#' )
#' ggchord(one_seq) +
#'   geom_seq(geometry = "double", thickness = 3) +
#'   geom_primer(data = primers, thickness = 0.6, id = "primers") +
#'   coord_circular()
geom_primer <- function(mapping = NULL, data = NULL, shape = "primer_arc",
    thickness = 2, stat = "identity", position = position_offset(0), ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_annotation_layer(mapping, data, "interval", shape, !missing(shape), thickness,
    stat, position, list(...), na.rm, show.legend, inherit.aes, id,
    feature_role = "primer")
}
