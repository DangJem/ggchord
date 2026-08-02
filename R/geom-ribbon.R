# geom-ribbon.R - alignment ribbon layer
# Fetches pre-computed ribbon polygon data from the package environment and renders it with geom_polygon
# Ribbon parameters are specified in this layer and stored for use at print time

# ---------------------------------------------------------------------------
# Internal helper: clone a polygon geom whose fill aesthetic is exposed under a
# different name.  The ribbon layer uses this so that its fill scale can be
# independent from the gene layer's fill scale (the gene layer keeps the plain
# "fill" aesthetic).  This replicates what the ggnewscale package used to do.
# ---------------------------------------------------------------------------
rename_fill_geom <- function(geom = GeomPolygon, old_aes = "fill",
                             new_aes = "zfill") {
  new_geom <- ggplot2::ggproto(
    paste0("GeomChord", sub("^Geom", "", class(geom)[1])), geom
  )

  # Rename the fill entry in the geom's default aesthetics so that no plain
  # "fill" default is injected into the ribbon data.
  aes_names <- names(new_geom$default_aes)
  aes_names[aes_names == old_aes] <- new_aes
  names(new_geom$default_aes) <- aes_names

  # The draw code (and legend keys) look for the original aesthetic name;
  # rename the data column back just before drawing.
  old_handle_na <- geom$handle_na
  new_geom$handle_na <- function(self, data, params) {
    colnames(data)[colnames(data) == new_aes] <- old_aes
    old_handle_na(data, params)
  }
  old_draw_key <- geom$draw_key
  new_geom$draw_key <- function(data, params, size) {
    colnames(data)[colnames(data) == new_aes] <- old_aes
    old_draw_key(data, params, size)
  }
  new_geom
}

# The geom used by the ribbon layer; its fill aesthetic is "fill_ribbon".
ribbon_geom <- rename_fill_geom()

#' Add an alignment ribbon layer
#'
#' Draws colored ribbons corresponding to alignment results. Color scheme and spacing parameters are specified here.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param ribbon_color_scheme Character. Color scheme "pident", "query" or "single", default "pident"
#' @param ribbon_colors Optional color vector. Ribbon color parameters
#' @param ribbon_alpha Numeric (0-1). Ribbon transparency, default 0.35
#' @param ribbon_ctrl_point Optional vector/list. Bezier control points, default c(0,0)
#' @param ribbon_gap Optional numeric/vector. Spacing between sequences and ribbons, default 0.15
#' @param alpha Ribbon transparency (overrides ribbon_alpha), defaults to the value used in the layout
#' @param ribbon_outline_color Character. Color of the ribbon outline (border), default "black"
#' @param ribbon_outline_width Numeric. Line width of the ribbon outline, default 0.05
#' @param ribbon_outline_linetype Numeric or character. Line type of the ribbon outline, default 1 (solid); see \code{linetype} in ggplot2 for options
#' @param show_legend Whether to show the legend, default TRUE
#' @param legend_position Position of this layer's legend (the Identity(%)
#'   colourbar): one of "left", "right", "top", "bottom" or "inside", default
#'   "left". Pass NULL to let the legend follow
#'   \code{theme(legend.position = ...)} together with the other legends. Can also be set with
#'   \code{theme(legend.position.ribbon = ...)}.
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
                        ribbon_outline_color = "black",
                        ribbon_outline_width = 0.05,
                        ribbon_outline_linetype = 1,
                        show_legend = TRUE,
                        legend_position = "left",
                        ...) {
  ribbon_alpha <- alpha %||% ribbon_alpha

  # Placeholder data
  empty_polys <- data.frame(
    x = numeric(0), y = numeric(0),
    group = integer(0),
    # Keep both possible mapped columns on the placeholder so that the layer
    # is valid until print.ggchord() injects the real data.
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

  # Ribbon outline (border) styling: passed as fixed colour/linewidth
  # aesthetics of the polygon layer when the user provides them.
  outline_params <- list()
  if (!is.null(ribbon_outline_color)) outline_params$colour <- ribbon_outline_color
  if (!is.null(ribbon_outline_width)) outline_params$linewidth <- ribbon_outline_width
  if (!is.null(ribbon_outline_linetype)) outline_params$linetype <- ribbon_outline_linetype

  lyr <- ggplot2::layer(
    data        = empty_polys,
    mapping     = fill_aes,
    stat        = "identity",
    geom        = ribbon_geom,
    position    = "identity",
    show.legend = if (identical(show_legend, TRUE)) {
                    c(fill = TRUE, colour = FALSE)
                  } else show_legend,
    inherit.aes = FALSE,
    check.aes   = FALSE,
    check.param = FALSE,
    key_glyph   = key_glyph_ribbon,
    params      = c(list(...), outline_params)
  )
  lyr$ggchord_type <- "ribbon"
  lyr$ggchord_params <- list(
    type                    = "ribbon",
    ribbon_color_scheme     = ribbon_color_scheme,
    ribbon_colors           = ribbon_colors,
    ribbon_alpha            = ribbon_alpha,
    ribbon_ctrl_point       = ribbon_ctrl_point,
    ribbon_gap              = ribbon_gap,
    ribbon_outline_color    = ribbon_outline_color,
    ribbon_outline_width    = ribbon_outline_width,
    ribbon_outline_linetype = ribbon_outline_linetype,
    legend_position         = legend_position
  )
  list(lyr)
}
