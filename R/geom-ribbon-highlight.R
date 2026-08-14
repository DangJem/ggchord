# geom-ribbon-highlight.R - ribbon highlight layer (v0.9.0)

#' Highlight selected alignment ribbons
#'
#' Draws a second polygon layer on top of selected ribbons so they can be
#' emphasized without changing the underlying Identity(%) legend. Selection is
#' done with safe, explicit filters (row numbers, query/subject IDs, pident and
#' length ranges) or a predicate function.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param ribbon_ids Optional integer vector of original ribbon row numbers to
#'   highlight.
#' @param qaccver Optional character vector; only ribbons whose query ID is in
#'   this set are highlighted.
#' @param saccver Optional character vector; only ribbons whose subject ID is in
#'   this set are highlighted.
#' @param min_pident Optional numeric. Minimum percent identity.
#' @param max_pident Optional numeric. Maximum percent identity.
#' @param min_length Optional numeric. Minimum alignment length.
#' @param max_length Optional numeric. Maximum alignment length.
#' @param predicate Optional function taking the ribbon data.frame and returning
#'   a logical vector with one element per row. Evaluated safely (no string
#'   parsing).
#' @param highlight_color Character. Highlight fill colour, default \code{"#E11D48"}.
#' @param highlight_alpha Numeric (0-1). Highlight alpha, default 0.8.
#' @param highlight_outline_color Optional outline colour, default \code{NULL}.
#' @param highlight_outline_width Numeric. Outline width, default 0.3.
#' @param highlight_expand Numeric. Fraction by which to expand the highlight
#'   polygon (currently reserved; 0 keeps the exact ribbon geometry).
#' @param show_legend Logical. Whether to show a legend, default FALSE.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A list of ggplot2 layers
#' @export
geom_ribbon_highlight <- function(mapping = NULL, data = NULL,
                                  ribbon_ids = NULL,
                                  qaccver = NULL,
                                  saccver = NULL,
                                  min_pident = NULL,
                                  max_pident = NULL,
                                  min_length = NULL,
                                  max_length = NULL,
                                  predicate = NULL,
                                  highlight_color = "#E11D48",
                                  highlight_alpha = 0.8,
                                  highlight_outline_color = NULL,
                                  highlight_outline_width = 0.3,
                                  highlight_expand = 0,
                                  show_legend = FALSE,
                                  ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!is.null(predicate) && !is.function(predicate)) {
    ggchord_stop("predicate must be a function taking ribbon_data and returning a logical vector")
  }
  if (!is.numeric(highlight_alpha) || length(highlight_alpha) != 1 ||
      !is.finite(highlight_alpha) || highlight_alpha < 0 || highlight_alpha > 1) {
    ggchord_stop("highlight_alpha must be in [0, 1]")
  }

  empty_polys <- data.frame(
    x = numeric(0), y = numeric(0), group = integer(0),
    stringsAsFactors = FALSE
  )
  lyr <- ggplot2::layer(
    data        = empty_polys,
    mapping     = aes(x = x, y = y, group = group),
    stat        = "identity",
    geom        = GeomPolygon,
    position    = "identity",
    show.legend = show_legend,
    inherit.aes = FALSE,
    check.aes   = FALSE,
    check.param = FALSE,
    params      = c(list(...),
                    list(colour = highlight_outline_color %||% NA,
                         linewidth = highlight_outline_width))
  )
  lyr$ggchord_type <- "ribbon_highlight"
  lyr$ggchord_params <- list(
    type                     = "ribbon_highlight",
    ribbon_ids               = ribbon_ids,
    qaccver                  = qaccver,
    saccver                  = saccver,
    min_pident               = min_pident,
    max_pident               = max_pident,
    min_length               = min_length,
    max_length               = max_length,
    predicate                = predicate,
    highlight_color          = highlight_color,
    highlight_alpha          = highlight_alpha,
    highlight_outline_color  = highlight_outline_color,
    highlight_outline_width  = highlight_outline_width,
    highlight_expand         = highlight_expand
  )
  list(lyr)
}
