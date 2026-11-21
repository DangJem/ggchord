#' Draw HSP alignments as Dotplot segments
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Observation table or a data function. For ggchord layers, provide the required data-role columns.
#' @param stat A compatible ggplot2 Stat or its name. The default uses the supplied records.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details Source references are `qaccver` and `saccver`. Point links consume `qpos` and `spos`. Interval links and alignment segments consume closed interval endpoints on both references. Supply reliable `qstrand` and `sstrand` or explicit relative direction when the geometry requires it.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' ggchord(seqs) +
#'   geom_alignment_segment(aes(colour = pident), data = links) +
#'   layout_dotplot(x_order = "A", y_order = "B") +
#'   coord_dotplot(ratio = 1) +
#'   scale_colour_gradient(limits = c(0, 100))
geom_alignment_segment <- function(mapping = NULL, data = NULL, stat = "identity",
    position = "identity", ...,
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  caller <- "geom_alignment_segment()"
  if (!ggchord_identity_stat(stat)) ggchord_stop(caller, ": use an identity Stat")
  if (!ggchord_is_data_identity(position)) ggchord_stop(caller, ": use an identity Position")
  dots <- ggchord_colour_dots(list(...), caller)
  ggchord_check_geom_params(dots, GeomAlignmentSegment, caller)
  layer <- ggplot2::layer(
    data = data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()),
    mapping = aes(x = x, y = y, xend = xend, yend = yend),
    geom = GeomAlignmentSegment, stat = "identity", position = "identity",
    show.legend = show.legend, inherit.aes = inherit.aes, check.aes = FALSE,
    params = c(list(na.rm = na.rm), dots)
  )
  layer$ggchord_type <- "alignment"
  layer$ggchord_params <- list(id = id, type = "alignment")
  ggchord_capture_layer_input(layer, data, mapping,
    c("qaccver", "qstart", "qend", "saccver", "sstart", "send", "qstrand", "sstrand", "direction"))
}

GeomAlignmentSegment <- ggproto("GeomAlignmentSegment", ggplot2::GeomSegment,
  default_aes = aes(colour = "#59636D", alpha = .6, linewidth = .4, linetype = 1))
