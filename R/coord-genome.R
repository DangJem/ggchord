# Dedicated coordinate contract for one circular genome.

#' Single-genome circular coordinate system
#'
#' Provides an explicit layout for one complete sequence. The sequence can be
#' drawn as a closed circle or with a controlled opening, while existing gene,
#' feature, region, axis and label layers remain reusable.
#'
#' @param gap Opening in degrees, from 0 (closed) up to but not including 180.
#' @param rotation Angle in degrees at which the sequence origin is placed.
#'   The default, 90, places it at the top of the circle.
#' @param direction Genomic direction around the circle.
#' @inheritParams coord_chord
#'
#' @return A Coord object for ggplot2 composition.
#' @export
#'
#' @examples
#' genome <- data.frame(accver = "genome", length = 5000)
#' ggchord(genome) + geom_seq() + coord_genome(gap = 12)
coord_genome <- function(gap = 0, rotation = 90,
                         direction = c("clockwise", "counterclockwise"),
                         ratio = 1, xlim = NULL, ylim = NULL,
                         expand = FALSE, clip = "off",
                         fit = c("labels", "geometry", "manual")) {
  direction <- match.arg(direction)
  if (!is.numeric(gap) || length(gap) != 1L || !is.finite(gap) ||
      gap < 0 || gap >= 180) {
    ggchord_stop("coord_genome(): gap must be one finite number in [0, 180)")
  }
  coord <- coord_chord(
    rotation = rotation, ratio = ratio, xlim = xlim, ylim = ylim,
    expand = expand, clip = clip, fit = fit
  )
  class(coord) <- unique(c("CoordGenome", class(coord)))
  coord$ggchord_genome <- TRUE
  coord$genome_gap <- gap
  coord$genome_direction <- direction
  coord
}
