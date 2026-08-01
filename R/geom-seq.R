# geom-seq.R - sequence arc layer
# Fetches pre-computed sequence arc data from the package environment and renders it with geom_path
# Sequence layout parameters are specified in this layer and stored for use at print time

#' Add a sequence arc layer
#'
#' Draws arcs (or straight lines, depending on the curvature setting) representing sequences in the chord diagram.
#' Sequence layout parameters (order, orientation, radius, curvature, colors, etc.) are specified here.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param seq_order Optional character vector. Specifies the drawing order of sequences
#' @param seq_labels Optional character vector or named vector. Sequence labels
#' @param seq_orientation Optional numeric (1 or -1). Sequence orientation, default 1
#' @param seq_gap Optional numeric. Gap proportion between sequences, default 0.03
#' @param seq_radius Optional numeric (> 0). Sequence arc radius, default 1.0
#' @param seq_curvature Optional numeric. Arc curvature (0=straight, 1=standard arc, >1=more curved), default 1.0
#' @param seq_colors Optional color vector or named vector. Sequence colors
#' @param linewidth Arc line width, default 1.2
#' @param show_legend Whether to show the legend for this layer, default TRUE
#' @param ... Additional arguments passed to \code{geom_path()}
#'
#' @return A list of ggplot2 layers
#' @export
geom_seq <- function(mapping = NULL, data = NULL,
                     seq_order = NULL,
                     seq_labels = NULL,
                     seq_orientation = NULL,
                     seq_gap = NULL,
                     seq_radius = NULL,
                     seq_curvature = NULL,
                     seq_colors = NULL,
                     linewidth = 1.2,
                     show_legend = TRUE,
                     ...) {
  set_seq_params(list(
    seq_order      = seq_order,
    seq_labels     = seq_labels,
    seq_orientation = seq_orientation,
    seq_gap        = seq_gap,
    seq_radius     = seq_radius,
    seq_curvature  = seq_curvature,
    seq_colors     = seq_colors
  ))

  # Only return the layer (scales are set uniformly by print.ggchord)
  list(
    ggplot2::layer(
      data        = data.frame(x = numeric(0), y = numeric(0),
                               seq_id = character(0)),
      mapping     = aes(x = x, y = y, group = seq_id, color = seq_id),
      stat        = "identity",
      geom        = GeomPath,
      position    = "identity",
      show.legend = if (identical(show_legend, TRUE)) {
                      c(colour = TRUE, fill = FALSE)
                    } else show_legend,
      inherit.aes = FALSE,
      key_glyph   = key_glyph_seq,
      params      = list(
        linewidth = linewidth,
        arrow = grid::arrow(type = "closed", length = grid::unit(3, "mm")),
        ...
      )
    )
  )
}
