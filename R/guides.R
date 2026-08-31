# guides.R - thin ggplot2 guide wrappers

#' Guides for ggchord role aesthetics
#'
#' Thin wrappers around [ggplot2::guide_legend()] and
#' [ggplot2::guide_colourbar()] with compact defaults suitable for
#' chord diagrams. All ordinary ggplot2 guide arguments remain available.
#'
#' @param title Guide title.
#' @param theme Optional guide-specific theme.
#' @param position,direction Guide position and direction.
#' @param override.aes A list of legend-key aesthetic overrides.
#' @param nrow,ncol Legend key layout.
#' @param reverse Whether to reverse the key order.
#' @param order Guide order.
#' @param ... Additional arguments passed to the corresponding ggplot2 guide.
#'
#' @return A ggplot2 guide object.
#' @export
guide_ggchord_legend <- function(
    title = waiver(), theme = NULL, position = NULL, direction = NULL,
    override.aes = list(), nrow = NULL, ncol = NULL, reverse = FALSE,
    order = 0, ...) {
  compact <- ggplot2::theme(
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(5.5, "mm"),
    legend.text = element_text(margin = margin(l = 1.2))
  )
  if (!is.null(theme)) compact <- compact + theme
  guide_legend(
    title = title, theme = compact, position = position,
    direction = direction, override.aes = override.aes,
    nrow = nrow, ncol = ncol, reverse = reverse, order = order, ...
  )
}

#' @rdname guide_ggchord_legend
#' @param nbin Number of colourbar bins.
#' @param display Colourbar display method.
#' @param alpha Colourbar alpha override.
#' @param draw.ulim,draw.llim Whether to draw upper/lower limit ticks.
#' @param angle Tick-label angle.
#' @param available_aes Aesthetics accepted by the guide.
#' @export
guide_ggchord_colourbar <- function(
    title = waiver(), theme = NULL, nbin = NULL, display = "raster",
    alpha = NA, draw.ulim = TRUE, draw.llim = TRUE, angle = NULL,
    position = NULL, direction = NULL, reverse = FALSE, order = 0,
    available_aes = c(
      "colour", "color", "fill", "ribbon_fill", "gene_fill",
      "feature_fill", "region_fill"
    ), ...) {
  horizontal <- identical(direction, "horizontal") ||
    (!is.null(position) && position %in% c("top", "bottom"))
  compact <- ggplot2::theme(
    legend.title.position = "top",
    legend.key.width = unit(if (horizontal) 42 else 3, "mm"),
    legend.key.height = unit(if (horizontal) 3 else 42, "mm"),
    legend.ticks.length = unit(1, "mm"),
    legend.text = element_text(margin = margin(l = 1.2)),
    legend.title = element_text(margin = margin(b = 1.5))
  )
  if (!is.null(theme)) compact <- compact + theme
  guide_colourbar(
    title = title, theme = compact, nbin = nbin, display = display,
    alpha = alpha, draw.ulim = draw.ulim, draw.llim = draw.llim,
    angle = angle, position = position, direction = direction,
    reverse = reverse, order = order, available_aes = available_aes, ...
  )
}
