# coord-chord.R - chord diagram coordinate system

#' Chord diagram coordinate system
#'
#' Controls global rotation, aspect ratio, coordinate limits, clipping and the
#' strategy used to fit chord geometry and labels.
#'
#' @param rotation Global clockwise layout rotation in degrees, default 45.
#' @param ratio Fixed y/x aspect ratio, default 1.
#' @param xlim,ylim Optional user limits. Explicit limits take priority over
#'   automatically fitted limits.
#' @param expand Logical. Expand coordinate limits, default TRUE.
#' @param clip Whether drawing is clipped to the panel, default \code{"off"}.
#' @param fit Fitting strategy: \code{"labels"} includes measured label boxes,
#'   \code{"geometry"} fits geometric elements only, and \code{"manual"}
#'   requires explicit \code{xlim} and \code{ylim}.
#'
#' @return A Coord object for ggplot2 \code{+} composition
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + coord_chord() + geom_seq()
#' p
coord_chord <- function(rotation = 45, ratio = 1,
                        xlim = NULL, ylim = NULL,
                        expand = TRUE, clip = "off",
                        fit = c("labels", "geometry", "manual")) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  fit <- match.arg(fit)
  if (!is.numeric(rotation) || length(rotation) != 1 || !is.finite(rotation)) {
    ggchord_stop("coord_chord(): rotation must be a finite number")
  }
  if (!is.numeric(ratio) || length(ratio) != 1 || !is.finite(ratio) || ratio <= 0) {
    ggchord_stop("coord_chord(): ratio must be a finite positive number")
  }
  validate_limit <- function(x, name) {
    if (is.null(x)) return(invisible(NULL))
    if (!is.numeric(x) || length(x) != 2 || any(!is.finite(x)) || x[1] == x[2]) {
      ggchord_stop("coord_chord(): ", name,
                   " must be NULL or two distinct finite numbers")
    }
  }
  validate_limit(xlim, "xlim")
  validate_limit(ylim, "ylim")
  if (!is.logical(expand) || length(expand) != 1 || is.na(expand)) {
    ggchord_stop("coord_chord(): expand must be TRUE or FALSE")
  }
  if (!clip %in% c("on", "off")) {
    ggchord_stop("coord_chord(): clip must be 'on' or 'off'")
  }
  if (fit == "manual" && (is.null(xlim) || is.null(ylim))) {
    ggchord_stop("coord_chord(): fit = 'manual' requires xlim and ylim")
  }

  coord <- coord_fixed(
    ratio = ratio, xlim = xlim, ylim = ylim,
    expand = expand, clip = clip
  )
  coord$ggchord_coord <- TRUE
  coord$rotation <- rotation
  coord$fit <- fit
  coord$user_xlim <- xlim
  coord$user_ylim <- ylim
  coord
}
