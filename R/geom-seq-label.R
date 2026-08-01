# geom-seq-label.R - sequence label layer
# Places sequence labels on (or outside) the sequence arcs.  The positions are
# computed at build time from the chord layout and injected into this layer.

#' Add a sequence label layer
#'
#' Places sequence labels at the midpoint of each sequence arc, radially
#' offset from the arc. Labels can be styled via their radial offset, rotation
#' and font size.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param seq_label_radius Optional numeric/vector. Radial position of the
#'   labels as a multiplier of the sequence arc radius (e.g. 1.15 = 15% outside
#'   the arc), default NULL (1.15)
#' @param seq_label_rotation Optional numeric/vector. Additional label rotation
#'   (degrees), default NULL (0)
#' @param seq_label_size Optional numeric/vector. Label font size, default NULL (3)
#' @param seq_labels Optional character vector. Override the label texts
#'   (defaults to the sequence labels from \code{geom_seq()} or the sequence IDs)
#' @param show_legend Whether to show the legend, default FALSE
#' @param ... Additional arguments passed to \code{geom_text()}
#'
#' @return A list of ggplot2 layers
#' @export
geom_seq_label <- function(mapping = NULL, data = NULL,
                           seq_label_radius = NULL,
                           seq_label_rotation = NULL,
                           seq_label_size = NULL,
                           seq_labels = NULL,
                           show_legend = FALSE,
                           ...) {
  lyr <- ggplot2::geom_text(
    data = data.frame(text_x = numeric(0), text_y = numeric(0),
                      label = character(0), text_angle = numeric(0),
                      size = numeric(0), hjust = numeric(0),
                      vjust = numeric(0), seq_label = character(0)),
    mapping = aes(x = text_x, y = text_y, label = label,
                  angle = text_angle, hjust = hjust, vjust = vjust,
                  size = size),
    inherit.aes = FALSE,
    show.legend = show_legend,
    ...
  )
  lyr$ggchord_params <- list(
    type              = "seq_label",
    seq_label_radius  = seq_label_radius,
    seq_label_rotation = seq_label_rotation,
    seq_label_size    = seq_label_size,
    seq_labels        = seq_labels
  )
  list(lyr)
}
