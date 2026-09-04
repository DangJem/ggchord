# geom-seq.R - sequence arc layer
# Fetches pre-computed sequence arc data from the package environment and renders it with geom_path
# Sequence layout parameters are specified in this layer and stored for use at print time

seq_geom <- rename_geom_aes(ggplot2::GeomPath, renames = c(colour = "seq_colour"))

#' Add a sequence arc layer
#'
#' Draws arcs (or straight lines, depending on the curvature setting) representing sequences in the chord diagram.
#' Sequence layout parameters (order, orientation, radius, curvature, colors, etc.) are specified here.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param seq_order Optional character vector. Specifies the drawing order of sequences
#' @param seq_orientation Optional numeric (1 or -1). Sequence orientation, default 1
#' @param seq_gap Optional numeric. Gap proportion between sequences, default 0.03
#' @param seq_radius Optional numeric (> 0). Sequence arc radius, default 1.0
#' @param seq_curvature Optional numeric. Arc curvature (0=straight, 1=standard arc, >1=more curved), default 1.0
#' @param linewidth Arc line width, default 0.9
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_path()}
#'
#' @return A ggplot2 layer.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + geom_seq()
#' p
geom_seq <- function(mapping = NULL, data = NULL,
                     seq_order = NULL,
                     seq_orientation = NULL,
                     seq_gap = NULL,
                     seq_radius = NULL,
                     seq_curvature = NULL,
                     linewidth = 0.9,
                     position = "identity",
                     show.legend = TRUE,
                     inherit.aes = FALSE,
                     ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  dots <- ggchord_colour_dots(dots, "geom_seq()")
  ggchord_reject_retired(dots, "geom_seq()", c(
    seq_labels = "geom_seq_label(labels = ...) or a label scale",
    seq_colors = "scale_seq_colour_manual(values = ...)",
    show_legend = "show.legend",
    legend_position = "guides(seq_colour = guide_ggchord_legend(position = ...))"
  ))
  removed_group_args <- intersect(
    names(dots),
    c("seq_group", "seq_group_gap", "seq_group_labels",
      "seq_group_label_radius", "seq_group_colors")
  )
  if (length(removed_group_args) > 0L) {
    ggchord_stop(
      "Sequence grouping was removed in v0.10.0; remove argument(s): ",
      paste(removed_group_args, collapse = ", ")
    )
  }
  removed_group_aes <- intersect(
    names(mapping), c("seq_group", "group_colour")
  )
  if (length(removed_group_aes) > 0L) {
    ggchord_stop(
      "Sequence grouping was removed in v0.10.0; remove aesthetic(s): ",
      paste(removed_group_aes, collapse = ", ")
    )
  }

  # The layout is computed at build time (ggplot_build.ggchord). The
  # parameters are attached to the layer itself so that the plot object is
  # fully self-contained.
  lyr <- ggplot2::layer(
    data        = data.frame(x = numeric(0), y = numeric(0),
                             seq_id = character(0)),
    mapping     = ggplot2::aes(x = x, y = y, group = seq_id, seq_colour = seq_id),
    stat        = "identity",
    geom        = seq_geom,
    position    = position,
    show.legend = if (identical(show.legend, TRUE)) {
                    c(seq_colour = TRUE, fill = FALSE)
                  } else show.legend,
    inherit.aes = inherit.aes,
    check.param = FALSE,
    key_glyph   = key_glyph_seq,
    params      = c(list(
      linewidth = linewidth,
      arrow = grid::arrow(type = "closed", length = grid::unit(2.4, "mm"))
    ), dots)
  )
  lyr$ggchord_type <- "seq"
  lyr$ggchord_params <- list(
    type                  = "seq",
    seq_order             = seq_order,
    seq_labels            = NULL,
    seq_orientation       = seq_orientation,
    seq_gap               = seq_gap,
    seq_radius            = seq_radius,
    seq_curvature         = seq_curvature,
    seq_colors            = NULL,
    legend_position       = NULL
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping, c("seq_id", "length", "seq_ring")
  )
  lyr
}
