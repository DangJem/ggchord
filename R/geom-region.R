# geom-region.R - sequence-region highlight layer (v0.9.0)

region_geom <- ggplot2::ggproto(
  "GeomChordRegion", GeomPolygon
)
aes_names <- names(region_geom$default_aes)
aes_names[aes_names == "fill"] <- "zregionfill"
names(region_geom$default_aes) <- aes_names
region_geom$handle_na <- function(self, data, params) {
  colnames(data)[colnames(data) == "zregionfill"] <- "fill"
  GeomPolygon$handle_na(data, params)
}
region_geom$draw_key <- function(data, params, size) {
  colnames(data)[colnames(data) == "zregionfill"] <- "fill"
  GeomPolygon$draw_key(data, params, size)
}

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
    alpha = numeric(0),
    stringsAsFactors = FALSE
  )

  lyr <- ggplot2::layer(
    data        = empty_polys,
    mapping     = aes(x = x, y = y, group = group,
                      zregionfill = zregionfill, alpha = alpha),
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
  list(lyr)
}
