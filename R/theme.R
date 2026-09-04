# Complete ggchord themes and chord-specific appearance settings.

#' ggchord themes
#'
#' Complete themes for chord diagrams. Besides ordinary ggplot2 theme
#' elements, these functions control automatic sequence axes, chord labels and
#' the appearance and placement of guides belonging to each ggchord role.
#' Palettes, breaks, labels and guide ordering remain scale/guide concerns.
#'
#' @param text,plot.title,plot.subtitle,plot.caption,plot.tag Text elements.
#' @param plot.title.position,plot.caption.position,plot.tag.position Standard
#'   ggplot2 title, caption and tag positions.
#' @param plot.background,panel.background Background elements.
#' @param plot.margin Plot margin.
#' @param axis `NULL` to show automatic sequence axes, or
#'   [ggplot2::element_blank()] to hide the complete axis decoration.
#' @param axis.line,axis.ticks,axis.minor.ticks Axis line elements.
#' @param axis.text Axis text element.
#' @param axis.gap,axis.ticks.length,axis.minor.ticks.length,axis.text.offset
#'   Scalar [grid::unit()] objects controlling physical axis spacing.
#' @param axis.text.orientation One of `"horizontal"`, `"parallel"`,
#'   `"perpendicular"`, or a finite numeric angle.
#' @param axis.text.check.overlap Whether overlapping axis labels are hidden.
#' @param seq.label,gene.label Text elements for existing label layers.
#' @param gene.label.segment Line element for existing gene leader lines.
#' @param legend.position,legend.justification,legend.direction,legend.box,legend.box.just,legend.box.margin,legend.box.background,legend.box.spacing,legend.background,legend.margin,legend.spacing,legend.spacing.x,legend.spacing.y,legend.key,legend.key.size,legend.key.width,legend.key.height,legend.text,legend.text.position,legend.title,legend.title.position Standard legend theme settings inherited by roles.
#' @param legend.seq,legend.ribbon,legend.gene,legend.feature,legend.region
#'   `NULL` to show and inherit common settings, or [ggplot2::element_blank()]
#'   to hide every guide for that role.
#' @param legend.seq.position,legend.ribbon.position,legend.gene.position,legend.feature.position,legend.region.position Per-role guide positions.
#' @param legend.seq.direction,legend.ribbon.direction,legend.gene.direction,legend.feature.direction,legend.region.direction Per-role directions.
#' @param legend.seq.background,legend.ribbon.background,legend.gene.background,legend.feature.background,legend.region.background Per-role backgrounds.
#' @param legend.seq.key,legend.ribbon.key,legend.gene.key,legend.feature.key,legend.region.key Per-role key elements.
#' @param legend.seq.key.size,legend.seq.key.width,legend.seq.key.height,legend.ribbon.key.size,legend.ribbon.key.width,legend.ribbon.key.height,legend.gene.key.size,legend.gene.key.width,legend.gene.key.height,legend.feature.key.size,legend.feature.key.width,legend.feature.key.height,legend.region.key.size,legend.region.key.width,legend.region.key.height Per-role key dimensions.
#' @param legend.seq.text,legend.ribbon.text,legend.gene.text,legend.feature.text,legend.region.text Per-role text elements.
#' @param legend.seq.text.position,legend.ribbon.text.position,legend.gene.text.position,legend.feature.text.position,legend.region.text.position Per-role text positions.
#' @param legend.seq.title,legend.ribbon.title,legend.gene.title,legend.feature.title,legend.region.title Per-role title elements.
#' @param legend.seq.title.position,legend.ribbon.title.position,legend.gene.title.position,legend.feature.title.position,legend.region.title.position Per-role title positions.
#' @param legend.seq.margin,legend.ribbon.margin,legend.gene.margin,legend.feature.margin,legend.region.margin Per-role margins.
#' @param legend.seq.spacing,legend.ribbon.spacing,legend.gene.spacing,legend.feature.spacing,legend.region.spacing Per-role spacing.
#' @param legend.ribbon.ticks,legend.ribbon.axis.line Ribbon colourbar lines.
#' @param legend.ribbon.ticks.length Ribbon colourbar tick length.
#' @param ... Additional arguments passed only to [ggplot2::theme()]. Removed
#'   underscore-style arguments produce a migration error.
#' @return A complete ggplot2 theme carrying ggchord appearance settings.
#' @export
theme_ggchord <- function(
    text = NULL,
    plot.title = NULL, plot.subtitle = NULL, plot.caption = NULL,
    plot.tag = NULL, plot.title.position = NULL,
    plot.caption.position = NULL, plot.tag.position = NULL,
    plot.background = NULL, panel.background = NULL, plot.margin = NULL,
    axis = NULL, axis.line = NULL, axis.ticks = NULL,
    axis.minor.ticks = NULL, axis.text = NULL,
    axis.gap = NULL, axis.ticks.length = NULL,
    axis.minor.ticks.length = NULL, axis.text.offset = NULL,
    axis.text.orientation = NULL, axis.text.check.overlap = NULL,
    seq.label = NULL, gene.label = NULL, gene.label.segment = NULL,
    legend.position = NULL, legend.justification = NULL,
    legend.direction = NULL, legend.box = NULL, legend.box.just = NULL,
    legend.box.margin = NULL, legend.box.background = NULL,
    legend.box.spacing = NULL, legend.background = NULL,
    legend.margin = NULL, legend.spacing = NULL, legend.spacing.x = NULL,
    legend.spacing.y = NULL, legend.key = NULL, legend.key.size = NULL,
    legend.key.width = NULL, legend.key.height = NULL, legend.text = NULL,
    legend.text.position = NULL, legend.title = NULL,
    legend.title.position = NULL,
    legend.seq = NULL, legend.seq.position = NULL,
    legend.seq.direction = NULL, legend.seq.background = NULL,
    legend.seq.key = NULL, legend.seq.key.size = NULL,
    legend.seq.key.width = NULL, legend.seq.key.height = NULL,
    legend.seq.text = NULL, legend.seq.text.position = NULL,
    legend.seq.title = NULL, legend.seq.title.position = NULL,
    legend.seq.margin = NULL, legend.seq.spacing = NULL,
    legend.ribbon = NULL, legend.ribbon.position = NULL,
    legend.ribbon.direction = NULL, legend.ribbon.background = NULL,
    legend.ribbon.key = NULL, legend.ribbon.key.size = NULL,
    legend.ribbon.key.width = NULL, legend.ribbon.key.height = NULL,
    legend.ribbon.text = NULL, legend.ribbon.text.position = NULL,
    legend.ribbon.title = NULL, legend.ribbon.title.position = NULL,
    legend.ribbon.margin = NULL, legend.ribbon.spacing = NULL,
    legend.ribbon.ticks = NULL, legend.ribbon.ticks.length = NULL,
    legend.ribbon.axis.line = NULL,
    legend.gene = NULL, legend.gene.position = NULL,
    legend.gene.direction = NULL, legend.gene.background = NULL,
    legend.gene.key = NULL, legend.gene.key.size = NULL,
    legend.gene.key.width = NULL, legend.gene.key.height = NULL,
    legend.gene.text = NULL, legend.gene.text.position = NULL,
    legend.gene.title = NULL, legend.gene.title.position = NULL,
    legend.gene.margin = NULL, legend.gene.spacing = NULL,
    legend.feature = NULL, legend.feature.position = NULL,
    legend.feature.direction = NULL, legend.feature.background = NULL,
    legend.feature.key = NULL, legend.feature.key.size = NULL,
    legend.feature.key.width = NULL, legend.feature.key.height = NULL,
    legend.feature.text = NULL, legend.feature.text.position = NULL,
    legend.feature.title = NULL, legend.feature.title.position = NULL,
    legend.feature.margin = NULL, legend.feature.spacing = NULL,
    legend.region = NULL, legend.region.position = NULL,
    legend.region.direction = NULL, legend.region.background = NULL,
    legend.region.key = NULL, legend.region.key.size = NULL,
    legend.region.key.width = NULL, legend.region.key.height = NULL,
    legend.region.text = NULL, legend.region.text.position = NULL,
    legend.region.title = NULL, legend.region.title.position = NULL,
    legend.region.margin = NULL, legend.region.spacing = NULL,
    ...) {
  ggchord_theme_from_environment("default", environment(), list(...))
}

#' @noRd
ggchord_theme_variant <- function(preset) {
  out <- theme_ggchord
  body(out) <- substitute(
    ggchord_theme_from_environment(PRESET, environment(), list(...)),
    list(PRESET = preset)
  )
  environment(out) <- environment(theme_ggchord)
  out
}

#' @rdname theme_ggchord
#' @export
theme_ggchord_minimal <- ggchord_theme_variant("minimal")
#' @rdname theme_ggchord
#' @export
theme_ggchord_dark <- ggchord_theme_variant("dark")
#' @rdname theme_ggchord
#' @export
theme_ggchord_publication <- ggchord_theme_variant("publication")

#' @noRd
ggchord_theme_from_environment <- function(preset, env, dots) {
  old <- c(
    base_size = "text = element_text(size = ...)",
    base_family = "text = element_text(family = ...)",
    axis_line = "axis.line", axis_ticks = "axis.ticks",
    axis_text = "axis.text", seq_label = "seq.label",
    gene_label = "gene.label", gene_label_segment = "gene.label.segment"
  )
  caller <- paste0("theme_ggchord", if (preset == "default") "" else paste0("_", preset), "()")
  ggchord_reject_retired(dots, caller, old)
  nms <- setdiff(names(formals(theme_ggchord)), "...")
  values <- stats::setNames(lapply(nms, get, envir = env), nms)
  base_size <- if (preset == "publication") 9 else 11
  dark <- preset == "dark"
  fg <- if (dark) "#F1F3F5" else "#30353A"
  bg <- if (dark) "#17191C" else "white"
  out <- ggplot2::theme_void(base_size = base_size, base_family = "") +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size, family = "", colour = fg),
      plot.title = ggplot2::element_text(hjust = .5, size = ggplot2::rel(1.27),
        face = "bold", colour = if (dark) fg else "#202428",
        margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(.9)),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(.72)),
      plot.background = ggplot2::element_rect(fill = bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = bg, colour = NA),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5),
      legend.background = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA),
      legend.box.spacing = grid::unit(if (preset == "minimal") 3 else 4, "mm"),
      legend.spacing = grid::unit(2, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.spacing.y = grid::unit(4, "mm"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(.73), colour = fg),
      legend.title = ggplot2::element_text(size = ggplot2::rel(.82), face = "bold", colour = fg),
      ggchord.axis.line = ggplot2::element_line(
        colour = if (dark) "#ADB5BD" else if (preset == "minimal") "#90969B" else "#737A80",
        linewidth = if (preset == "minimal") .25 else .3),
      ggchord.axis.ticks = ggplot2::element_line(
        colour = if (dark) "#ADB5BD" else if (preset == "minimal") "#A0A5AA" else "#8A9096",
        linewidth = if (preset == "minimal") .2 else .25),
      ggchord.axis.minor.ticks = ggplot2::element_line(
        colour = if (dark) "#8F979F" else "#A5AAAF", linewidth = .2),
      ggchord.axis.text = ggplot2::element_text(
        colour = if (dark) fg else "#3D4348", size = ggplot2::rel(.72)),
      ggchord.seq.label = ggplot2::element_text(
        colour = if (dark) fg else "#252A2E", size = ggplot2::rel(.8)),
      ggchord.gene.label = ggplot2::element_text(
        colour = if (dark) fg else "#2B3035", size = ggplot2::rel(.64)),
      ggchord.gene.label.segment = ggplot2::element_line(
        colour = if (dark) "#CED4DA" else if (preset == "minimal") "#A1A6AA" else "#858B91",
        linewidth = if (preset == "minimal") .2 else .25)
    )
  custom <- list(
    ggchord.axis.line = values$axis.line,
    ggchord.axis.ticks = values$axis.ticks,
    ggchord.axis.minor.ticks = values$axis.minor.ticks,
    ggchord.axis.text = values$axis.text,
    ggchord.seq.label = values$seq.label,
    ggchord.gene.label = values$gene.label,
    ggchord.gene.label.segment = values$gene.label.segment
  )
  standard <- c("text", "plot.title", "plot.subtitle", "plot.caption", "plot.tag",
    "plot.title.position", "plot.caption.position", "plot.tag.position",
    "plot.background", "panel.background", "plot.margin", "legend.position",
    "legend.justification", "legend.direction", "legend.box", "legend.box.just",
    "legend.box.margin", "legend.box.background", "legend.box.spacing",
    "legend.background", "legend.margin", "legend.spacing", "legend.spacing.x",
    "legend.spacing.y", "legend.key", "legend.key.size", "legend.key.width",
    "legend.key.height", "legend.text", "legend.text.position", "legend.title",
    "legend.title.position")
  overrides <- c(values[standard], custom, dots)
  overrides <- overrides[!vapply(overrides, is.null, logical(1))]
  if (length(overrides)) out <- out + do.call(ggplot2::theme, overrides)
  attr(out, "ggchord.settings") <- ggchord_theme_settings(values)
  out
}

#' @noRd
ggchord_theme_settings <- function(values) {
  if (!is.null(values$axis) && !inherits(values$axis, "element_blank"))
    ggchord_stop("theme_ggchord(): axis must be NULL or element_blank()")
  unit_nms <- c("axis.gap", "axis.ticks.length", "axis.minor.ticks.length",
                "axis.text.offset", "legend.ribbon.ticks.length")
  for (nm in unit_nms) {
    x <- values[[nm]]
    if (!is.null(x) && (!grid::is.unit(x) || length(x) != 1L))
      ggchord_stop("theme_ggchord(): ", nm, " must be one grid::unit()")
    if (!is.null(x) && !grid::unitType(x) %in%
        c("inches", "cm", "mm", "points", "bigpts", "picas"))
      ggchord_stop("theme_ggchord(): ", nm, " must use an absolute physical unit")
  }
  orientation <- values$axis.text.orientation %||% "parallel"
  if (is.character(orientation)) {
    if (length(orientation) != 1L || is.na(orientation) ||
        !orientation %in% c("horizontal", "parallel", "perpendicular"))
      ggchord_stop("theme_ggchord(): invalid axis.text.orientation")
  } else if (!is.numeric(orientation) || length(orientation) != 1L || !is.finite(orientation)) {
    ggchord_stop("theme_ggchord(): axis.text.orientation must be a supported name or finite angle")
  }
  overlap <- values$axis.text.check.overlap %||% FALSE
  if (!is.logical(overlap) || length(overlap) != 1L || is.na(overlap))
    ggchord_stop("theme_ggchord(): axis.text.check.overlap must be TRUE or FALSE")
  fields <- c("position", "direction", "background", "key", "key.size",
    "key.width", "key.height", "text", "text.position", "title",
    "title.position", "margin", "spacing")
  roles <- c("seq", "ribbon", "gene", "feature", "region")
  role_settings <- stats::setNames(vector("list", length(roles)), roles)
  allowed <- c("left", "right", "top", "bottom", "inside", "none")
  for (role in roles) {
    parent <- values[[paste0("legend.", role)]]
    if (!is.null(parent) && !inherits(parent, "element_blank"))
      ggchord_stop("theme_ggchord(): legend.", role, " must be NULL or element_blank()")
    item <- list(hidden = inherits(parent, "element_blank"))
    for (field in fields) item[[field]] <- values[[paste0("legend.", role, ".", field)]]
    if (!is.null(item$position) && (!is.character(item$position) ||
        length(item$position) != 1L || is.na(item$position) || !item$position %in% allowed))
      ggchord_stop("theme_ggchord(): invalid legend.", role, ".position")
    role_settings[[role]] <- item
  }
  role_settings$ribbon$ticks <- values$legend.ribbon.ticks
  role_settings$ribbon$ticks.length <- values$legend.ribbon.ticks.length
  role_settings$ribbon$axis.line <- values$legend.ribbon.axis.line
  common_names <- c("position", "justification", "direction", "box",
    "box.just", "box.margin", "box.background", "box.spacing",
    "background", "margin", "spacing", "spacing.x", "spacing.y", "key",
    "key.size", "key.width", "key.height", "text", "text.position",
    "title", "title.position")
  common <- stats::setNames(lapply(common_names, function(x) {
    values[[paste0("legend.", x)]]
  }), common_names)
  list(axis = list(
    hidden = inherits(values$axis, "element_blank"),
    gap = values$axis.gap %||% grid::unit(.8, "mm"),
    ticks.length = values$axis.ticks.length %||% grid::unit(.5, "mm"),
    minor.ticks.length = values$axis.minor.ticks.length %||% grid::unit(.25, "mm"),
    text.offset = values$axis.text.offset %||% grid::unit(.8, "mm"),
    text.orientation = orientation, text.check.overlap = overlap),
    legend = common, legends = role_settings)
}

#' Convert an absolute unit to inches without opening a graphics device
#' @noRd
ggchord_unit_inches <- function(x) {
  value <- as.numeric(x)
  switch(grid::unitType(x), inches = value, cm = value / 2.54,
    mm = value / 25.4, points = value / 72.27,
    bigpts = value / 72, picas = value / 6,
    ggchord_stop("ggchord axis units must be absolute physical units"))
}

#' @noRd
ggchord_plot_settings <- function(plot) {
  found <- attr(plot$theme, "ggchord.settings")
  if (!is.null(found)) return(found)
  nms <- setdiff(names(formals(theme_ggchord)), "...")
  ggchord_theme_settings(stats::setNames(rep(list(NULL), length(nms)), nms))
}

#' @noRd
ggchord_theme_element <- function(plot, name) {
  tryCatch(ggplot2::calc_element(name, plot$theme), error = function(e) NULL)
}
#' @noRd
ggchord_theme_text_size <- function(plot, name, fallback) {
  el <- ggchord_theme_element(plot, name)
  if (is.null(el) || inherits(el, "element_blank") || is.null(el@size) ||
      length(el@size) != 1L || !is.finite(el@size)) return(fallback)
  el@size / ggplot2::.pt
}
#' @noRd
ggchord_theme_point_size <- function(plot, name, fallback) {
  el <- ggchord_theme_element(plot, name)
  if (is.null(el) || inherits(el, "element_blank") || is.null(el@size) ||
      length(el@size) != 1L || !is.finite(el@size)) return(fallback)
  el@size
}

#' @noRd
ggchord_apply_theme_styles <- function(plot) {
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    components <- lyr$ggchord_theme_components
    if (!is.null(components)) {
      for (param_name in names(components)) {
        el <- ggchord_theme_element(plot, unname(components[[param_name]]))
        if (is.null(el)) next
        current <- lyr$geom_params[[param_name]] %||% list()
        if (inherits(el, "element_blank")) current$alpha <- 0
        else if (inherits(el, "element_line")) {
          vals <- list(colour = el@colour, linewidth = el@linewidth,
            linetype = el@linetype, lineend = el@lineend, linejoin = el@linejoin)
          for (nm in names(vals)) if (is.null(current[[nm]]) && !is.null(vals[[nm]])) current[[nm]] <- vals[[nm]]
        } else if (inherits(el, "element_text")) {
          vals <- list(colour = el@colour, family = el@family,
            fontface = el@face, lineheight = el@lineheight)
          for (nm in names(vals)) if (is.null(current[[nm]]) && !is.null(vals[[nm]])) current[[nm]] <- vals[[nm]]
        }
        lyr$geom_params[[param_name]] <- current
      }
      plot$layers[[i]] <- lyr
      next
    }
    name <- lyr$ggchord_theme_element
    if (is.null(name)) next
    el <- ggchord_theme_element(plot, name)
    if (is.null(el)) next
    if (inherits(el, "element_blank")) {
      if (is.null(lyr$aes_params$alpha)) lyr$aes_params$alpha <- 0
      plot$layers[[i]] <- lyr
      next
    }
    if (inherits(el, "element_line")) {
      vals <- list(colour = el@colour, linewidth = el@linewidth, linetype = el@linetype)
      geom_vals <- list(lineend = el@lineend, linejoin = el@linejoin)
    } else if (inherits(el, "element_text")) {
      vals <- list(colour = el@colour, family = el@family,
        fontface = el@face, lineheight = el@lineheight)
      geom_vals <- list()
    } else next
    for (nm in names(vals)) if (is.null(lyr$aes_params[[nm]]) && !is.null(vals[[nm]])) lyr$aes_params[[nm]] <- vals[[nm]]
    for (nm in names(geom_vals)) if (is.null(lyr$geom_params[[nm]]) && !is.null(geom_vals[[nm]])) lyr$geom_params[[nm]] <- geom_vals[[nm]]
    plot$layers[[i]] <- lyr
  }
  plot
}
