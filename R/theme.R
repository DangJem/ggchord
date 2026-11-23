#' Common appearance for sequence visualizations
#'
#' A warm ivory (#FFFCF5), borderless panel and centred bold plot title keep the focus on
#' the sequence geometry. Standard ggplot2 theme elements can override them.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @param ... Must be empty.
#' @details Add `theme()` to configure registered `ggchord.*` elements or
#'   standard ggplot2 theme elements.
#' @return A ggplot2 theme.
#' @export
#' @examples
#' library(ggplot2)
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' p_basic + theme_ggchord(base_size = 12) +
#'   theme(legend.position = "bottom",
#'   ggchord.axis.text = element_text(colour = "#40545A"))
theme_ggchord <- function(base_size = 11, base_family = "", ...) {
  ggchord_local_console_errors()
  ggchord_check_theme_arguments(...)
  ggchord_make_theme(base_size, base_family, "default")
}

#' Minimal sequence-plot appearance
#' @inheritParams theme_ggchord
#' @return An ordinary ggplot2 theme.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' p_basic + theme_ggchord_minimal(base_size = 11)
theme_ggchord_minimal <- function(base_size = 11, base_family = "", ...) {
  ggchord_local_console_errors()
  ggchord_check_theme_arguments(...)
  ggchord_make_theme(base_size, base_family, "minimal")
}

#' Dark sequence-plot appearance
#' @inheritParams theme_ggchord
#' @return An ordinary ggplot2 theme.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' p_basic + theme_ggchord_dark(base_size = 11)
theme_ggchord_dark <- function(base_size = 11, base_family = "", ...) {
  ggchord_local_console_errors()
  ggchord_check_theme_arguments(...)
  ggchord_make_theme(base_size, base_family, "dark")
}

#' Publication sequence-plot appearance
#' @inheritParams theme_ggchord
#' @return An ordinary ggplot2 theme.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' p_basic + theme_ggchord_publication(base_size = 10)
theme_ggchord_publication <- function(base_size = 10, base_family = "", ...) {
  ggchord_local_console_errors()
  ggchord_check_theme_arguments(...)
  ggchord_make_theme(base_size, base_family, "publication")
}

#' Plasmid sequence-plot appearance
#'
#' The warm ivory background and dark sequence axis distinguish annotations
#' from major and minor ticks. Use `scale_seq_continuous(guide = "none")`
#' to omit the sequence guide.
#' @inheritParams theme_ggchord
#' @return An ordinary ggplot2 theme.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' ggchord(data = one_seq) + geom_seq() +
#'   geom_gene(data = subset(genes, accver == "A")) +
#'   coord_circular() + theme_ggchord_plasmid()
theme_ggchord_plasmid <- function(base_size = 11, base_family = "", ...) {
  ggchord_local_console_errors()
  ggchord_check_theme_arguments(...)
  ggchord_make_theme(base_size, base_family, "plasmid")
}

ggchord_check_theme_arguments <- function(...) {
  if (...length())
    ggchord_stop("Add + theme(...) to set theme elements. Presets accept base_size and base_family.")
}

ggchord_make_theme <- function(base_size, base_family, variant) {
  if (!is.numeric(base_size) || length(base_size) != 1L ||
      !is.finite(base_size) || base_size <= 0)
    ggchord_stop("base_size must be a positive font size in points")
  if (!is.character(base_family) || length(base_family) != 1L || is.na(base_family))
    ggchord_stop("base_family must be one font family")
  blank <- element_blank()
  dark <- variant == "dark"
  publication <- variant == "publication"
  plasmid <- variant == "plasmid"
  fg <- if (dark) "#F1F3F5" else if (publication) "#111111" else "#30353A"
  bg <- if (dark) "#17191C" else if (variant %in% c("default", "plasmid")) "#FFFCF5" else "white"
  axis_col <- if (dark) "#ADB5BD" else if (publication) "#333333" else if (variant == "default") "#60696F" else "#737A80"
  line <- element_line(colour = axis_col, linewidth = .3)
  axis_text <- element_text(colour = fg, size = rel(.8))
  title <- element_text(colour = fg, size = rel(.9))
  sequence_line <- if (plasmid) element_line(colour = "#252525", linewidth = .55) else line
  background <- element_rect(fill = bg, colour = NA)
  out <- ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(colour = fg),
      plot.title = element_text(hjust = .5, face = "bold", size = rel(1.27),
        margin = margin(b = 4)),
      plot.subtitle = element_text(size = rel(.9)),
      plot.caption = element_text(size = rel(.72)),
      panel.grid = blank,
      panel.background = background,
      plot.background = background,
      panel.border = blank,
      legend.position = if (variant %in% c("publication", "plasmid")) "bottom" else "right",
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.text = element_text(colour = fg, size = rel(.8)),
      legend.title = element_text(colour = fg, size = rel(.9), face = "bold"),
      axis.line = line, axis.ticks = line, axis.text = axis_text,
      axis.title = title,
      ggchord.axis.line = if (variant == "minimal") blank else sequence_line,
      ggchord.axis.ticks = if (plasmid) element_line(colour = "#252525", linewidth = .45) else line,
      ggchord.axis.minor.ticks = if (plasmid) element_line(colour = "#85898B", linewidth = .2) else line,
      ggchord.axis.ticks.length = unit(if (plasmid) 2.2 else 1.2, "mm"),
      ggchord.axis.text = if (plasmid) element_text(colour = "#252525", size = rel(.78)) else axis_text,
      ggchord.axis.title = if (variant == "plasmid") blank else title,
      ggchord.track.axis.line = if (variant == "minimal") blank else line,
      ggchord.track.axis.ticks = line,
      ggchord.track.axis.ticks.length = unit(1.2, "mm"),
      ggchord.track.axis.text = axis_text,
      ggchord.track.axis.title = title)
  if (variant == "default") out <- out + theme(
    legend.text = element_text(colour = fg, size = rel(.86)),
    legend.key.spacing.y = unit(1.5, "mm"))
  if (plasmid) out <- out + theme(ggchord.axis.gap = unit(0, "mm"),
    ggchord.axis.minor.ticks.length = unit(.8, "mm"), plot.margin = margin(6, 6, 6, 6))
  if (publication) out <- out + theme(plot.margin = margin(4, 4, 4, 4, unit = "mm"))
  attr(out, "ggchord.explicit_theme") <- list()
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
    if (!is.null(x) && (!is.unit(x) || length(x) != 1L || !is.finite(as.numeric(x))))
      ggchord_stop("theme_ggchord(): ", nm, " must be one grid::unit()")
    if (!is.null(x) && !unitType(x) %in%
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
  overlap <- values$axis.text.check.overlap %||% TRUE
  if (!is.logical(overlap) || length(overlap) != 1L || is.na(overlap))
    ggchord_stop("theme_ggchord(): axis.text.check.overlap must be TRUE or FALSE")
  fields <- c("position", "direction", "background", "key", "key.size",
    "key.width", "key.height", "text", "text.position", "title",
    "title.position", "margin", "spacing")
  roles <- c("seq", "link", "ribbon", "gene", "feature", "region")
  role_settings <- setNames(vector("list", length(roles)), roles)
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
  common <- setNames(lapply(common_names, function(x) {
    values[[paste0("legend.", x)]]
  }), common_names)
  list(axis = list(
    hidden = inherits(values$axis, "element_blank"),
    gap = values$axis.gap %||% unit(.8, "mm"),
    gap.explicit = !is.null(values$axis.gap),
    ticks.length = values$axis.ticks.length %||% unit(.85, "mm"),
    ticks.length.explicit = !is.null(values$axis.ticks.length),
    minor.ticks.length = values$axis.minor.ticks.length %||% unit(.4, "mm"),
    minor.ticks.length.explicit = !is.null(values$axis.minor.ticks.length),
    text.offset = values$axis.text.offset %||% unit(.8, "mm"),
    text.offset.explicit = !is.null(values$axis.text.offset),
    text.orientation = orientation, text.check.overlap = overlap),
    legend = common, legends = role_settings)
}

#' Convert an absolute unit to inches without opening a graphics device
#' @noRd
ggchord_unit_inches <- function(x) {
  value <- as.numeric(x)
  switch(unitType(x), inches = value, cm = value / 2.54,
    mm = value / 25.4, points = value / 72.27,
    bigpts = value / 72, picas = value / 6,
    ggchord_stop("ggchord axis units must be absolute physical units"))
}

#' @noRd
ggchord_plot_settings <- function(plot) {
  found <- attr(plot$theme, "ggchord.settings")
  if (is.null(found)) {
    found <- ggchord_theme_settings(list())
  }
  physical <- c("axis", "axis.gap", "axis.ticks.length", "axis.minor.ticks.length",
    "axis.text.offset", "axis.text.orientation", "axis.text.check.overlap")
  values <- setNames(lapply(physical, function(n) plot$theme[[paste0("ggchord.", n)]]), physical)
  if (any(!vapply(values, is.null, logical(1)))) {
    checked <- ggchord_theme_settings(values)$axis
    for (name in physical[!vapply(values, is.null, logical(1))]) {
      key <- if (name == "axis") "hidden" else sub("^axis[.]", "", name)
      found$axis[[key]] <- checked[[key]]
      if (key == "text.offset" && !isTRUE(attr(values[[name]],
          "ggchord.automatic_axis_offset"))) found$axis$text.offset.explicit <- TRUE
      if (key %in% c("gap", "ticks.length", "minor.ticks.length")) found$axis[[paste0(key, ".explicit")]] <- TRUE
    }
  }
  if (isTRUE(plot$coordinates$ggchord_circular)) {
    if (!isTRUE(found$axis$gap.explicit)) found$axis$gap <- unit(0, "mm")
    if (!isTRUE(found$axis$ticks.length.explicit))
      found$axis$ticks.length <- unit(1.5, "mm")
    if (!isTRUE(found$axis$minor.ticks.length.explicit))
      found$axis$minor.ticks.length <- unit(.7, "mm")
  }
  found
}

#' @noRd
ggchord_theme_element <- function(plot, name) {
  theme <- plot$theme
  if (name %in% c("ggchord.axis.ticks", "ggchord.axis.minor.ticks")) {
    element <- theme[[name]]
    if (inherits(element, "element_line") && is.null(element@linewidth)) {
      axis <- ggchord_plot_settings(plot)$axis
      widths <- if (isTRUE(plot$coordinates$ggchord_circular)) {
        axis$circular.tick.linewidth %||% c(major = .45, minor = .3)
      } else axis$tick.linewidth %||% c(major = .25, minor = .2)
      element@linewidth <- unname(widths[if (name == "ggchord.axis.ticks") "major" else "minor"])
      theme[[name]] <- element
    }
  }
  tryCatch(ggplot2::calc_element(name, theme), error = function(e) NULL)
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

ggchord_resolve_theme_component <- function(plot, lyr, param_name, element_name) {
  el <- ggchord_theme_element(plot, element_name)
  current <- lyr$geom_params[[param_name]] %||% list()
  if (identical(element_name, "ggchord.feature.label.segment")) {
    raw <- plot$theme[[element_name]]
    current$.automatic_feature_segment <- is.null(raw) ||
      isTRUE(attr(raw, "ggchord.automatic_feature_segment"))
  }
  if (is.null(el)) return(current)
  if (inherits(el, "element_blank")) current$alpha <- 0
  else if (inherits(el, "element_line")) {
    vals <- list(colour = el@colour, linewidth = el@linewidth,
      linetype = el@linetype, lineend = el@lineend, linejoin = el@linejoin)
    # A composite layer may map text colour while its segment component
    # still needs an independent theme colour. Component line styles
    # are applied after scale mapping, so a layer-level mapping must not
    # suppress them.
    for (nm in names(vals)) if (is.null(current[[nm]]) &&
        !is.null(vals[[nm]])) {
      current[[nm]] <- vals[[nm]]
    }
  } else if (inherits(el, "element_text")) {
    vals <- list(colour = el@colour, size = el@size / ggplot2::.pt,
      family = el@family,
      fontface = el@face, lineheight = el@lineheight)
    for (nm in names(vals)) if (is.null(current[[nm]]) &&
        is.null(lyr$mapping[[nm]]) && !is.null(vals[[nm]])) {
      current[[nm]] <- vals[[nm]]
    }
  }
  current
}

#' @noRd
ggchord_apply_theme_styles <- function(plot) {
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    if (identical(lyr$ggchord_type, "track")) {
      style <- ggchord_track_axis_style(plot, lyr$ggchord_params$track_id)
      for (name in c("axis_line_params", "axis_tick_params", "axis_text_params", "axis_title_params", "axis_check_overlap"))
        if (is.null(lyr$geom_params[[name]])) lyr$geom_params[[name]] <- style[[name]]
      plot$layers[[i]] <- lyr
    }
    components <- lyr$ggchord_theme_components
    if (!is.null(components)) {
      for (param_name in names(components)) {
        lyr$geom_params[[param_name]] <- ggchord_resolve_theme_component(
          plot, lyr, param_name, unname(components[[param_name]]))
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
      vals <- list(colour = el@colour, size = el@size / ggplot2::.pt,
        family = el@family,
        fontface = el@face, lineheight = el@lineheight)
      geom_vals <- list()
    } else next
    for (nm in names(vals)) {
      if (is.null(lyr$aes_params[[nm]]) && is.null(lyr$mapping[[nm]]) &&
          !is.null(vals[[nm]])) lyr$aes_params[[nm]] <- vals[[nm]]
    }
    for (nm in names(geom_vals)) if (is.null(lyr$geom_params[[nm]]) && !is.null(geom_vals[[nm]])) lyr$geom_params[[nm]] <- geom_vals[[nm]]
    plot$layers[[i]] <- lyr
  }
  plot
}
