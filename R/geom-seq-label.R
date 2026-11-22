#' Label reference sequences
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional reference table or function selecting from the registered reference.
#' @param anchor Position on the entity or reference path that anchors the text.
#' @param position Position adjustment. For bound labels, leave this unset to follow the target layer.
#' @param label_layout A fixed or repelled label layout declaration.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Labels read the registered reference. Use
#'   `label_layout_fixed(orientation = "tangent")` to align each label with
#'   its own backbone while keeping the text upright after a rotation or
#'   direction reversal. Set `angle` for a fixed text angle.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' ggchord(seqs) +
#'   geom_seq() +
#'   geom_seq_label(aes(label = accver), anchor = "start") +
#'   coord_collinear()
#'
#' ggchord(seqs) + geom_seq_arrow() +
#'   geom_seq_label(label_layout = label_layout_fixed(orientation = "tangent", offset = 4)) +
#'   coord_chord(start_angle = 40)
geom_seq_label <- function(mapping = NULL, data = NULL, anchor = "middle",
    position = "identity", label_layout = label_layout_fixed(orientation = "horizontal"),
    ..., na.rm = FALSE, show.legend = FALSE, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  anchor <- ggchord_match_arg(anchor, c("start", "middle", "end", "center"))
  ggchord_domain_label(mapping, data, "identity", NULL, anchor, position,
    label_layout, "sequence", list(...), na.rm, show.legend, inherit.aes, id)
}
