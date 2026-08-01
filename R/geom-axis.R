# geom-axis.R - axis layer
# Fetches pre-computed axis lines, tick marks, and label data from the package environment
# Axis parameters are specified in this layer and stored for use at print time

#' Add an axis layer
#'
#' Draws axes for each sequence in the chord diagram (including axis lines, major/minor ticks, and labels).
#' Axis parameters (spacing, tick count/length, label size/orientation, etc.) are specified here.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param show_axis Logical. Whether to show the axis, default TRUE
#' @param axis_gap Optional numeric/vector. Spacing between sequence and axis, default 0.04
#' @param axis_tick_major_number Optional integer/vector. Number of major ticks, default 5
#' @param axis_tick_major_length Optional numeric/vector. Major tick length ratio, default 0.02
#' @param axis_tick_minor_number Optional integer/vector. Number of minor ticks, default 4
#' @param axis_tick_minor_length Optional numeric/vector. Minor tick length ratio, default 0.01
#' @param axis_label_size Optional numeric/vector. Tick label font size, default 3
#' @param axis_label_offset Optional numeric/vector. Label offset ratio, default 1.5
#' @param axis_label_orientation Optional character/numeric/vector. Label orientation, default "horizontal"
#' @param show_legend Whether to show the legend, default FALSE (axes do not participate in legends)
#' @param ... Additional arguments passed to geom_path/geom_segment/geom_text
#'
#' @return A list of ggplot2 layers
#' @export
geom_axis <- function(mapping = NULL, data = NULL,
                      show_axis = NULL,
                      axis_gap = NULL,
                      axis_tick_major_number = NULL,
                      axis_tick_major_length = NULL,
                      axis_tick_minor_number = NULL,
                      axis_tick_minor_length = NULL,
                      axis_label_size = NULL,
                      axis_label_offset = NULL,
                      axis_label_orientation = NULL,
                      show_legend = FALSE,
                      ...) {
  empty_id <- data.frame(x = numeric(0), y = numeric(0),
                         seq_id = character(0))
  empty_seg <- data.frame(x0 = numeric(0), y0 = numeric(0),
                          x1 = numeric(0), y1 = numeric(0),
                          label = character(0), label_x = numeric(0),
                          label_y = numeric(0), size = numeric(0),
                          seq_id = character(0))

  path_layer <- geom_path(data = empty_id,
                          mapping = aes(x = x, y = y, group = seq_id),
                          color = "black", linewidth = 0.3,
                          inherit.aes = FALSE, show.legend = show_legend, ...)
  path_layer$ggchord_params <- list(
    type                    = "axis",
    show_axis               = show_axis,
    axis_gap                = axis_gap,
    axis_tick_major_number  = axis_tick_major_number,
    axis_tick_major_length  = axis_tick_major_length,
    axis_tick_minor_number  = axis_tick_minor_number,
    axis_tick_minor_length  = axis_tick_minor_length,
    axis_label_size         = axis_label_size,
    axis_label_offset       = axis_label_offset,
    axis_label_orientation  = axis_label_orientation
  )

  list(
    path_layer,
    geom_segment(data = empty_seg,
                 mapping = aes(x = x0, y = y0,
                               xend = x1, yend = y1),
                 color = "black", linewidth = 0.3,
                 inherit.aes = FALSE, show.legend = show_legend, ...),
    geom_text(data = empty_seg[integer(0), ],
              mapping = aes(x = label_x, y = label_y,
                            label = label, size = size),
              inherit.aes = FALSE, color = "black",
              angle = 0, show.legend = show_legend, ...)
  )
}
