# Primer binding-site geometry.

#' Draw primer binding sites
#'
#' A biological primer-binding layer built on [geom_feature()] and its compact
#' directional `primer_arrow`. Binding discovery remains separate: pass the
#' output of [find_primer_bindings()] here, and add
#' [geom_feature_label_repel()] with the same data when labels are wanted.
#'
#' @inheritParams geom_feature
#' @param feature_width Primer-band width.
#' @param feature_fill Default primer fill colour.
#' @return A ggplot2 layer.
#' @export
#' @examples
#' primer <- data.frame(
#'   accver = "circle", start = 100, end = 120,
#'   strand = "+", name = "sequencing primer"
#' )
#' ggchord(data.frame(accver = "circle", length = 1000)) +
#'   geom_seq() + geom_primer(data = primer) + coord_circular()
geom_primer <- function(mapping = NULL, data = NULL,
                        feature_width = .058,
                        feature_fill = "#A020F0",
                        arrow_head_length = .042,
                        arrow_head_width = 1.18,
                        position = NULL,
                        show.legend = FALSE,
                        inherit.aes = FALSE, ...) {
  if (is.null(position)) {
    position <- position_feature_stack(
      spacing = .10, base_position = position_plasmid()
    )
  }
  layer <- geom_feature(
    mapping = mapping, data = data, feature_shape = "primer_arrow",
    feature_width = feature_width, arrow_head_length = arrow_head_length,
    arrow_head_width = arrow_head_width, arrow_head_style = "shouldered",
    short_feature = "auto", position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    feature_fill = feature_fill, ...
  )
  inherited_transform <- layer$ggchord_input_transform
  layer$ggchord_input_transform <- function(x) {
    if (!"feature_type" %in% names(x)) x$feature_type <- "primer"
    if (!"feature_label" %in% names(x) && "name" %in% names(x)) {
      x$feature_label <- x$name
    }
    inherited_transform(x)
  }
  layer$ggchord_params$feature_role <- "primer"
  layer$ggchord_params$is_primer <- TRUE
  layer
}
