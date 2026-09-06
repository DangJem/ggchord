# geom-region.R - sequence-region highlight layer (v0.9.0)

region_geom <- rename_geom_aes(
  ggplot2::GeomPolygon, renames = c(fill = "region_fill")
)

#' Highlight regions along sequence arcs
#'
#' Draws rectangular bands on sequence arcs for one or more coordinate
#' intervals. This is useful for marking loci, repeats, CRISPR arrays, or other
#' user-defined regions without turning them into gene arrows.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data data.frame with at least \code{accver}, \code{start},
#'   \code{end}; optional \code{label}, \code{category} and \code{color}.
#' @param fill Character. Default fill colour for regions, default
#'   \code{"#F59E0B"}.
#' @param colour Outline colour, default \code{"#B45309"}.
#' @param alpha Numeric (0-1). Region alpha, default 0.25.
#' @param region_width Numeric. Band width in chord radius units, default 0.08.
#' @param region_offset Numeric. Radial offset from the sequence arc, default 0.
#' @param region_side Character. \code{"inside"}, \code{"outside"} or
#'   \code{"auto"}; \code{"auto"} places regions inside the chord when possible.
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A ggplot2 layer.
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' regions <- data.frame(accver = "MT108731.1",
#'                       start = 1000, end = 4000,
#'                       color = "orange")
#' p <- ggchord(seq_data_example) + geom_seq() +
#'   geom_seq_region(data = regions)
#' p
geom_seq_region <- function(mapping = NULL, data = NULL,
                            fill = "#F59E0B",
                            colour = "#B45309",
                            alpha = 0.25,
                            region_width = 0.08,
                            region_offset = 0,
                            region_side = c("inside", "outside", "auto"),
                            position = "identity",
                            show.legend = FALSE,
                            inherit.aes = FALSE,
                            ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  dots <- list(...)
  alias <- ggchord_colour_alias(colour, dots, "geom_seq_region()", sys.call())
  colour <- alias$colour
  dots <- alias$dots
  ggchord_reject_retired(dots, "geom_seq_region()", c(
    regions = "data",
    region_fill = "fill",
    region_color = "colour",
    region_alpha = "alpha",
    show_legend = "show.legend"
  ))
  region_side <- match.arg(region_side)
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      !is.finite(alpha) || alpha < 0 || alpha > 1) {
    ggchord_stop("alpha must be in [0, 1]")
  }
  if (!is.numeric(region_width) || length(region_width) != 1 ||
      !is.finite(region_width) || region_width <= 0) {
    ggchord_stop("region_width must be a finite positive number")
  }
  if (!is.numeric(region_offset) || length(region_offset) != 1 ||
      !is.finite(region_offset)) {
    ggchord_stop("region_offset must be a finite number")
  }

  empty_polys <- data.frame(
    x = numeric(0), y = numeric(0),
    group = integer(0),
    zregionfill = character(0),
    colour = character(0),
    alpha = numeric(0),
    stringsAsFactors = FALSE
  )

  lyr <- ggplot2::layer(
    data        = empty_polys,
    mapping     = ggplot2::aes(x = x, y = y, group = group,
                      region_fill = zregionfill, colour = colour,
                      alpha = alpha),
    stat        = "identity",
    geom        = region_geom,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    check.aes   = FALSE,
    check.param = FALSE,
    params      = dots
  )
  lyr$ggchord_type <- "seq_region"
  lyr$ggchord_params <- list(
    type           = "seq_region",
    regions        = data,
    region_fill    = fill,
    region_color   = colour,
    region_alpha   = alpha,
    region_width   = region_width,
    region_offset  = region_offset,
    region_side    = region_side
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping, c("accver", "start", "end")
  )
  lyr
}
