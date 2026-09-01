# guides.R - thin ggplot2 guide wrappers

#' Scale legend furniture relative to the active export device
#' @noRd
ggchord_device_scale <- function(device_inches = NULL,
                                 reference = c(width = 8, height = 6),
                                 limits = c(0.5, 1.6)) {
  if (is.null(device_inches)) {
    device_inches <- tryCatch(
      grDevices::dev.size("in"),
      error = function(e) c(NA_real_, NA_real_)
    )
  }
  if (!is.numeric(device_inches) || length(device_inches) != 2L ||
      any(!is.finite(device_inches)) || any(device_inches <= 0)) {
    return(1)
  }
  value <- min(device_inches / unname(reference))
  min(limits[2], max(limits[1], value))
}

ggchord_check_guide_scale <- function(size_scale, caller) {
  if (!is.numeric(size_scale) || length(size_scale) != 1L ||
      !is.finite(size_scale) || size_scale <= 0) {
    ggchord_stop(caller, ": size_scale must be one positive finite number")
  }
  size_scale
}

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
#' @param size_scale Positive multiplier for keys, spacing and typography.
#'   ggchord's automatically generated guides derive this from the active
#'   export device; direct calls default to 1.
#' @param ... Additional arguments passed to the corresponding ggplot2 guide.
#'
#' @return A ggplot2 guide object.
#' @export
guide_ggchord_legend <- function(
    title = ggplot2::waiver(), theme = NULL, position = NULL, direction = NULL,
    override.aes = list(), nrow = NULL, ncol = NULL, reverse = FALSE,
    order = 0, size_scale = 1, ...) {
  size_scale <- ggchord_check_guide_scale(
    size_scale, "guide_ggchord_legend()"
  )
  compact <- ggplot2::theme(
    legend.key = ggplot2::element_rect(fill = NA, colour = NA),
    legend.key.height = grid::unit(4 * size_scale, "mm"),
    legend.key.width = grid::unit(5.5 * size_scale, "mm"),
    legend.text = ggplot2::element_text(
      size = 8 * size_scale,
      margin = ggplot2::margin(l = 1.2 * size_scale)
    ),
    legend.title = ggplot2::element_text(size = 9 * size_scale),
    legend.margin = ggplot2::margin(
      t = size_scale, b = size_scale, unit = "mm"
    )
  )
  if (!is.null(theme)) compact <- compact + theme
  ggplot2::guide_legend(
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
    title = ggplot2::waiver(), theme = NULL, nbin = NULL, display = "raster",
    alpha = NA, draw.ulim = TRUE, draw.llim = TRUE, angle = NULL,
    position = NULL, direction = NULL, reverse = FALSE, order = 0,
    available_aes = c(
      "colour", "color", "fill", "ribbon_fill", "gene_fill",
      "feature_fill", "region_fill"
    ), size_scale = 1, ...) {
  size_scale <- ggchord_check_guide_scale(
    size_scale, "guide_ggchord_colourbar()"
  )
  horizontal <- identical(direction, "horizontal") ||
    (!is.null(position) && position %in% c("top", "bottom"))
  compact <- ggplot2::theme(
    legend.title.position = "top",
    legend.key.width = grid::unit(
      (if (horizontal) 50 else 3.6) * size_scale, "mm"
    ),
    legend.key.height = grid::unit(
      (if (horizontal) 3.6 else 50) * size_scale, "mm"
    ),
    legend.ticks.length = grid::unit(size_scale, "mm"),
    legend.text = ggplot2::element_text(
      size = 8 * size_scale,
      margin = ggplot2::margin(l = 1.2 * size_scale)
    ),
    legend.title = ggplot2::element_text(
      size = 9 * size_scale,
      margin = ggplot2::margin(b = 2 * size_scale, unit = "mm")
    ),
    legend.margin = ggplot2::margin(
      t = size_scale, b = size_scale, unit = "mm"
    )
  )
  if (!is.null(theme)) compact <- compact + theme
  ggplot2::guide_colourbar(
    title = title, theme = compact, nbin = nbin, display = display,
    alpha = alpha, draw.ulim = draw.ulim, draw.llim = draw.llim,
    angle = angle, position = position, direction = direction,
    reverse = reverse, order = order, available_aes = available_aes, ...
  )
}
