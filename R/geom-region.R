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
#' @param data Optional data.frame; an alias for \code{regions} when
#'   \code{regions} is NULL.
#' @param regions data.frame with at least \code{seq_id}, \code{start},
#'   \code{end}; optional \code{label}, \code{category} and \code{color}.
#' @param region_fill Character. Default fill colour for regions, default
#'   \code{"#F59E0B"}.
#' @param region_color Character. Outline colour, default \code{"#B45309"}.
#' @param region_alpha Numeric (0-1). Region alpha, default 0.25.
#' @param region_width Numeric. Band width in chord radius units, default 0.08.
#' @param region_offset Numeric. Radial offset from the sequence arc, default 0.
#' @param region_side Character. \code{"inside"}, \code{"outside"} or
#'   \code{"auto"}; \code{"auto"} places regions inside the chord when possible.
#' @param show_legend Logical. Whether to show a legend for category colours,
#'   default FALSE.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A list of ggplot2 layers
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' regions <- data.frame(seq_id = "MT108731.1",
#'                       start = 1000, end = 4000,
#'                       color = "orange")
#' p <- ggchord(seq_data_example) + geom_seq() + geom_seq_region(regions)
#' p
geom_seq_region <- function(mapping = NULL, data = NULL,
                            regions = NULL,
                            region_fill = "#F59E0B",
                            region_color = "#B45309",
                            region_alpha = 0.25,
                            region_width = 0.08,
                            region_offset = 0,
                            region_side = c("inside", "outside", "auto"),
                            show_legend = FALSE,
                            ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  # Preserve the documented pre-v0.9 positional `geom_seq_region(regions)`
  # call while keeping the standard mapping/data argument order.
  if (is.data.frame(mapping) && is.null(data) && is.null(regions)) {
    regions <- mapping
    mapping <- NULL
  }
  regions <- regions %||% data
  region_side <- match.arg(region_side)
  if (!is.numeric(region_alpha) || length(region_alpha) != 1 ||
      !is.finite(region_alpha) || region_alpha < 0 || region_alpha > 1) {
    ggchord_stop("region_alpha must be in [0, 1]")
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
    position    = "identity",
    show.legend = show_legend,
    inherit.aes = FALSE,
    check.aes   = FALSE,
    check.param = FALSE,
    params      = list(...)
  )
  lyr$ggchord_type <- "seq_region"
  lyr$ggchord_params <- list(
    type           = "seq_region",
    regions        = regions,
    region_fill    = region_fill,
    region_color   = region_color,
    region_alpha   = region_alpha,
    region_width   = region_width,
    region_offset  = region_offset,
    region_side    = region_side
  )
  lyr <- ggchord_capture_layer_input(
    lyr, regions, mapping, c("seq_id", "start", "end")
  )
  list(lyr)
}
