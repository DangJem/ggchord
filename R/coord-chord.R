# coord-chord.R - lightweight chord diagram coordinate system
# Uses ggplot2's built-in coord_fixed
# print.ggchord replaces the coord dynamically; this function only provides an initial placeholder

#' Chord diagram coordinate system
#'
#' A lightweight Coord. It creates placeholder coordinates in \code{ggchord()},
#' which are replaced at print time with the actual extents computed from the layout.
#'
#' @param layout Chord layout object (passed internally by ggchord(), may be NULL)
#'
#' @return A Coord object for ggplot2 \code{+} composition
#' @export
coord_chord <- function(layout = NULL) {
  # Initial placeholder range; replaced by print.ggchord at print time
  coord_fixed(
    ratio = 1,
    xlim  = c(-5, 5),
    ylim  = c(-5, 5),
    clip  = "off"
  )
}
