#' Label the centre of a circular reference
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional reference table or function selecting from the registered reference.
#' @param position Position adjustment. For bound labels, leave this unset to follow the target layer.
#' @param label_layout A fixed or repelled label layout declaration.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID. Use a name for explicit references.
#' @return A ggplot2 layer for addition to a plot.
#' @details The source is the complete reference table with `accver` and `length`. A backbone can select an authoritative subset. An explicit complete table in geom_seq() can register a reference when the plot has no global data.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' ggchord(one_seq) +
#'   geom_seq() +
#'   geom_seq_center_label(aes(label = accver), size = 5) +
#'   coord_circular()
geom_seq_center_label <- function(mapping = NULL, data = NULL, position = "identity",
    label_layout = label_layout_fixed(orientation = "horizontal"), ...,
    na.rm = FALSE, show.legend = FALSE, inherit.aes = TRUE, id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  geom_seq_label(mapping, data, anchor = "center", position = position,
    label_layout = label_layout, ..., na.rm = na.rm, show.legend = show.legend,
    inherit.aes = inherit.aes, id = id)
}
