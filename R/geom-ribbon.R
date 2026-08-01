# geom-ribbon.R - alignment ribbon layer
# Fetches pre-computed ribbon polygon data from the package environment and renders it with geom_polygon
# Ribbon parameters are specified in this layer and stored for use at print time

#' Add an alignment ribbon layer
#'
#' Draws colored ribbons corresponding to BLAST alignment results. Color scheme and spacing parameters are specified here.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param ribbon_color_scheme Character. Color scheme "pident", "query" or "single", default "pident"
#' @param ribbon_colors Optional color vector. Ribbon color parameters
#' @param ribbon_alpha Numeric (0-1). Ribbon transparency, default 0.35
#' @param ribbon_ctrl_point Optional vector/list. Bezier control points, default c(0,0)
#' @param ribbon_gap Optional numeric/vector. Spacing between sequences and ribbons, default 0.15
#' @param alpha Ribbon transparency (overrides ribbon_alpha), defaults to the value used in the layout
#' @param show_legend Whether to show the legend, default TRUE
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A list of ggplot2 layers
#' @export
geom_ribbon <- function(mapping = NULL, data = NULL,
                        ribbon_color_scheme = NULL,
                        ribbon_colors = NULL,
                        ribbon_alpha = NULL,
                        ribbon_ctrl_point = NULL,
                        ribbon_gap = NULL,
                        alpha = NULL,
                        show_legend = TRUE,
                        ...) {
  set_ribbon_params(list(
    ribbon_color_scheme = ribbon_color_scheme,
    ribbon_colors       = ribbon_colors,
    ribbon_alpha        = alpha %||% ribbon_alpha,
    ribbon_ctrl_point   = ribbon_ctrl_point,
    ribbon_gap          = ribbon_gap
  ))

  # Placeholder data
  empty_polys <- data.frame(
    x = numeric(0), y = numeric(0),
    group = integer(0),
    # ggnewscale builds the plot as soon as a new fill scale is added.
    # Keep both possible mapped columns on the placeholder so that this
    # preliminary build succeeds before print.ggchord() injects real data.
    pident = numeric(0),
    fill = character(0),
    alpha = numeric(0),
    stringsAsFactors = FALSE
  )

  # Detect ribbon_color_scheme to decide the fill mapping variable
  # The real data is injected at print time and the fill column is replaced
  scheme <- ribbon_color_scheme %||% "pident"

  if (scheme == "pident") {
    fill_aes <- aes(x = x, y = y, group = group, fill = pident, alpha = alpha)
  } else {
    fill_aes <- aes(x = x, y = y, group = group, fill = fill, alpha = alpha)
  }

  list(
    ggplot2::layer(
      data        = empty_polys,
      mapping     = fill_aes,
      stat        = "identity",
      geom        = GeomPolygon,
      position    = "identity",
      show.legend = if (identical(show_legend, TRUE)) {
                      c(fill = TRUE, colour = FALSE)
                    } else show_legend,
      inherit.aes = FALSE,
      key_glyph   = key_glyph_ribbon,
      params      = list(...)
    )
  )
}
