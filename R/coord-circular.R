# Dedicated coordinate contract for one circular sequence.

CoordCircular <- ggplot2::ggproto("CoordCircular", CoordGgchord)

#' Single-sequence circular coordinate system
#'
#' @param rotation Angle in degrees at which the genomic origin is placed.
#' @param direction Direction of increasing genomic coordinates.
#' @param gap Opening in degrees. Zero draws a closed circle.
#' @inheritParams coord_chord
#' @return A Coord object for ggplot2 composition.
#' @export
coord_circular <- function(
    rotation = 0,
    direction = c("clockwise", "counterclockwise"),
    gap = 0,
    ratio = 1,
    xlim = NULL, ylim = NULL,
    expand = FALSE, clip = "off",
    fit = c("labels", "geometry", "manual")) {
  direction <- match.arg(direction)
  fit <- match.arg(fit)
  if (!is.numeric(rotation) || length(rotation) != 1L || !is.finite(rotation)) {
    ggchord_stop("coord_circular(): rotation must be one finite number")
  }
  if (!is.numeric(gap) || length(gap) != 1L || !is.finite(gap) ||
      gap < 0 || gap >= 360) {
    ggchord_stop("coord_circular(): gap must be one finite number in [0, 360)")
  }
  validated <- coord_chord(
    rotation = rotation, ratio = ratio, xlim = xlim, ylim = ylim,
    expand = expand, clip = clip, fit = fit
  )
  coord <- new_ggchord_coord(
    CoordCircular, rotation, ratio, xlim, ylim, expand, clip, fit,
    validated$user_xlim, validated$user_ylim
  )
  coord$ggchord_circular <- TRUE
  coord$circular_gap <- as.numeric(gap)
  coord$circular_direction <- direction
  coord
}
