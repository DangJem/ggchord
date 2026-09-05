# guides.R - thin ggplot2 guide wrappers

#' Scale legend furniture relative to the active export device
#' @noRd
ggchord_device_scale <- function(device_inches = NULL,
                                 reference = c(width = 8, height = 6),
                                 limits = c(0.5, 1.6)) {
  if (is.null(device_inches)) {
    device_inches <- if (grDevices::dev.cur() == 1L) {
      c(NA_real_, NA_real_)
    } else tryCatch(
      grDevices::dev.size("in"), error = function(e) c(NA_real_, NA_real_)
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

#' @rdname guide_ggchord_legend
#' @export
guide_ggchord_colorbar <- guide_ggchord_colourbar

#' Resolve a guide theme and placement for one chord role
#' @noRd
ggchord_role_guide_spec <- function(plot, role, colourbar = FALSE) {
  settings <- ggchord_plot_settings(plot)
  common <- settings$legend
  item <- settings$legends[[role]]
  # Keep legends at the outer sides of the plotting area by default. The
  # continuous ribbon/Identity guide uses the left edge; compact categorical
  # guides use the right. An explicitly supplied common or role position still
  # takes precedence, including an ordinary theme(legend.position = ...).
  default_position <- if (identical(role, "ribbon")) "left" else "right"
  resolved_position <- item$position %||% common$position %||%
    default_position
  global_position <- common$position %||% plot$theme$legend.position
  hidden <- identical(global_position, "none") || isTRUE(item$hidden) ||
    identical(item$position, "none")
  if (hidden) return(list(hidden = TRUE))

  size_scale <- ggchord_device_scale()
  field_map <- c(
    background = "legend.background", key = "legend.key",
    key.size = "legend.key.size", key.width = "legend.key.width",
    key.height = "legend.key.height", text = "legend.text",
    text.position = "legend.text.position", title = "legend.title",
    title.position = "legend.title.position", margin = "legend.margin",
    spacing = "legend.spacing"
  )
  vals <- list()
  for (field in names(field_map)) {
    value <- item[[field]] %||% common[[field]]
    if (!is.null(value)) vals[[field_map[[field]]]] <- value
  }
  # key.width/key.height take precedence over key.size. Responsive defaults
  # are used only when neither the role nor the common theme supplied units.
  key_size <- item$key.size %||% common$key.size
  explicit_width <- item$key.width %||% common$key.width %||%
    plot$theme$legend.key.width
  explicit_height <- item$key.height %||% common$key.height %||%
    plot$theme$legend.key.height
  if (!is.null(explicit_width)) vals$legend.key.width <- explicit_width
  if (!is.null(explicit_height)) vals$legend.key.height <- explicit_height
  if (is.null(explicit_width) && !is.null(key_size))
    vals$legend.key.width <- key_size
  if (is.null(explicit_height) && !is.null(key_size))
    vals$legend.key.height <- key_size
  resolved_direction <- item$direction %||% common$direction
  horizontal <- colourbar && (
    identical(resolved_direction, "horizontal") ||
    resolved_position %in% c("top", "bottom") ||
      identical(common$box, "horizontal")
  )
  if (is.null(vals$legend.key.width))
    vals$legend.key.width <- grid::unit(
      if (colourbar && horizontal) 54 else if (colourbar) 3.8 else 5.5, "mm"
    ) * size_scale
  if (is.null(vals$legend.key.height))
    vals$legend.key.height <- grid::unit(
      if (colourbar && !horizontal) 54 else if (colourbar) 3.8 else 4, "mm"
    ) * size_scale
  raw_text <- item$text %||% common$text
  raw_title <- item$title %||% common$title
  plot_text <- plot$theme$legend.text
  plot_title <- plot$theme$legend.title
  if (is.null(raw_text) && inherits(plot_text, "element_text") &&
      !is.null(plot_text@size) && !inherits(plot_text@size, "rel")) raw_text <- plot_text
  if (is.null(raw_title) && inherits(plot_title, "element_text") &&
      !is.null(plot_title@size) && !inherits(plot_title@size, "rel")) raw_title <- plot_title
  responsive_text <- is.null(raw_text) ||
    (inherits(raw_text, "element_text") &&
       !is.null(raw_text@size) && inherits(raw_text@size, "rel"))
  responsive_title <- is.null(raw_title) ||
    (inherits(raw_title, "element_text") &&
       !is.null(raw_title@size) && inherits(raw_title@size, "rel"))
  if (responsive_text) {
    vals$legend.text <- ggplot2::element_text(
      size = ggchord_theme_point_size(plot, "legend.text", 8) * size_scale
    )
  }
  if (responsive_title) {
    vals$legend.title <- ggplot2::element_text(
      size = ggchord_theme_point_size(plot, "legend.title", 9) * size_scale
    )
  }
  if (colourbar) {
    if (!is.null(item[["ticks"]])) vals$legend.ticks <- item[["ticks"]]
    if (!is.null(item[["ticks.length"]]))
      vals$legend.ticks.length <- item[["ticks.length"]]
    if (!is.null(item[["axis.line"]]))
      vals$legend.axis.line <- item[["axis.line"]]
  }
  list(
    hidden = FALSE,
    position = resolved_position,
    direction = resolved_direction %||% NULL,
    theme = if (length(vals)) do.call(ggplot2::theme, vals) else NULL,
    size_scale = size_scale
  )
}

#' Construct the default guide for one chord role
#' @noRd
ggchord_role_guide <- function(plot, role, colourbar = FALSE, order = 0,
                               override.aes = list()) {
  spec <- ggchord_role_guide_spec(plot, role, colourbar)
  if (isTRUE(spec$hidden)) return("none")
  if (colourbar) {
    guide_ggchord_colourbar(
      position = spec$position, direction = spec$direction,
      theme = spec$theme, size_scale = spec$size_scale, order = order,
      available_aes = "ribbon_fill"
    )
  } else {
    guide_ggchord_legend(
      position = spec$position, direction = spec$direction,
      theme = spec$theme, size_scale = spec$size_scale, order = order,
      override.aes = override.aes
    )
  }
}
