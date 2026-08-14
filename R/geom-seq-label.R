# geom-seq-label.R - sequence label layer
# Places sequence labels on (or outside) the sequence arcs.  The positions are
# computed at build time from the chord layout and injected into this layer.

#' Add a sequence label layer
#'
#' Places sequence labels at the midpoint of each sequence arc, radially
#' offset from the arc. Labels can be styled via their radial offset, rotation,
#' font size, text orientation and justification.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param seq_label_radius Optional numeric/vector. Radial position of the
#'   labels as a multiplier of the sequence arc radius: \code{1} is on the arc,
#'   \code{> 1} places the label outside (away from the chord center) and
#'   \code{< 1} places it inside, default 1
#' @param seq_label_rotation Optional numeric/vector. Additional label rotation
#'   (degrees) on top of the arc-aligned orientation, default NULL (0). Ignored
#'   when \code{seq_label_orientation = "horizontal"}.
#' @param seq_label_size Optional numeric/vector. Label font size, default NULL (3)
#' @param seq_labels Optional character vector. Override the label texts
#'   (defaults to the sequence labels from \code{geom_seq()} or the sequence IDs)
#' @param seq_label_orientation Character, default "arc". Label text
#'   orientation: \code{"arc"} rotates the text along the sequence arc (and
#'   keeps it readable), \code{"horizontal"} draws every label horizontally,
#'   extending away from the chord center.
#' @param seq_label_hjust Optional numeric/vector. Horizontal justification of
#'   the labels, default NULL. The default arc orientation uses -0.2 so the
#'   text sits just inside the sequence; with
#'   \code{seq_label_orientation = "horizontal"} the justification is chosen
#'   automatically so the text extends away from the chord center.
#' @param seq_label_vjust Optional numeric/vector. Vertical justification of
#'   the labels, default NULL (0.5).
#' @param check_overlap Logical, default FALSE. When TRUE, labels that would
#'   overlap a previously drawn label are skipped (ggplot2's
#'   \code{geom_text()} option).
#' @param show_legend Whether to show the legend, default FALSE
#' @param ... Additional arguments passed to \code{geom_text()}
#'
#' @return A list of ggplot2 layers
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + geom_seq() + geom_seq_label()
#' p
geom_seq_label <- function(mapping = NULL, data = NULL,
                           seq_label_radius = 1,
                           seq_label_rotation = NULL,
                           seq_label_size = NULL,
                           seq_labels = NULL,
                           seq_label_orientation = c("arc", "horizontal"),
                           seq_label_hjust = NULL,
                           seq_label_vjust = NULL,
                           check_overlap = FALSE,
                           show_legend = FALSE,
                           ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  seq_label_orientation <- match.arg(seq_label_orientation)
  lyr <- ggplot2::geom_text(
    data = data.frame(text_x = numeric(0), text_y = numeric(0),
                      label = character(0), text_angle = numeric(0),
                      size = numeric(0), hjust = numeric(0),
                      vjust = numeric(0)),
    mapping = aes(x = text_x, y = text_y, label = label,
                  angle = text_angle, hjust = hjust, vjust = vjust,
                  size = size),
    inherit.aes = FALSE,
    show.legend = show_legend,
    check_overlap = check_overlap,
    ...
  )
  lyr$ggchord_type <- "seq_label"
  lyr$ggchord_params <- list(
    type              = "seq_label",
    seq_label_radius  = seq_label_radius,
    seq_label_rotation = seq_label_rotation,
    seq_label_size    = seq_label_size,
    seq_labels        = seq_labels,
    seq_label_orientation = seq_label_orientation,
    seq_label_hjust   = seq_label_hjust,
    seq_label_vjust   = seq_label_vjust
  )
  list(lyr)
}
