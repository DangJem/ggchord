#' Common appearance for sequence visualizations
#'
#' A warm ivory (#FFFCF5), borderless panel and centred bold plot title keep the focus on
#' the sequence geometry. Standard ggplot2 theme elements can override them.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @details Add `theme()` to configure standard ggplot2 elements or the
#'   registered sequence and track axis decorations. Set reference axis
#'   position, visibility and label orientation with [scale_seq_continuous()]
#'   and [guide_seq_axis()]. Set data label and leader appearance on the label
#'   geometry. To enable the ggchord palettes, add [scales_ggchord()].
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
theme_ggchord <- function(base_size = 11, base_family = "") {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
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
theme_ggchord_minimal <- function(base_size = 11, base_family = "") {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
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
theme_ggchord_dark <- function(base_size = 11, base_family = "") {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
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
theme_ggchord_publication <- function(base_size = 10, base_family = "") {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
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
theme_ggchord_plasmid <- function(base_size = 11, base_family = "") {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  ggchord_make_theme(base_size, base_family, "plasmid")
}

ggchord_check_theme_elements <- function(theme) {
  keys <- names(theme) %||% character()
  keys <- keys[startsWith(keys, "ggchord.")]
  if (!length(keys)) return(invisible(NULL))
  migration <- c(
    ggchord.axis = paste0('scale_seq_continuous(guide = "none")',
      ' or use guide = "none" on the native x and y scales with coord_dotplot()'),
    ggchord.axis.gap = "scale_seq_continuous(guide = guide_seq_axis(gap = ...))",
    ggchord.axis.text.orientation =
      "scale_seq_continuous(guide = guide_seq_axis(orientation = ...)) or guide_seq_axis(angle = ...) for degrees",
    ggchord.axis.text.check.overlap = "scale_seq_continuous(guide = guide_seq_axis(check.overlap = ...))",
    ggchord.seq.label = "geom_seq_label(size = ..., colour = ..., family = ...)",
    ggchord.gene.label = "geom_gene_label(size = ..., colour = ..., family = ...)",
    ggchord.feature.label = "geom_feature_label(size = ..., colour = ..., family = ...)",
    ggchord.restriction.label = "geom_site_label(size = ..., colour = ..., family = ...)",
    ggchord.gene.label.segment = "geom_gene_label(leader = ..., segment_colour = ..., segment_linewidth = ...)",
    ggchord.feature.label.segment = "geom_feature_label(leader = ..., segment_colour = ..., segment_linewidth = ...)",
    ggchord.restriction.label.segment = "geom_site_label(leader = ..., segment_colour = ..., segment_linewidth = ...)",
    ggchord.seq.center.label = "geom_seq_label_center(size = ..., colour = ..., family = ...)",
    ggchord.seq.center.name = "geom_seq_label_center(aes(label = ...), size = ...)",
    ggchord.seq.center.length = "geom_seq_label_center(aes(label = ...), size = ...)"
  )
  retired <- intersect(keys, names(migration))
  if (length(retired)) {
    name <- retired[1L]
    replacement <- migration[[name]]
    if (endsWith(name, ".segment") && inherits(theme[[name]], "element_blank"))
      replacement <- sub("[(].*$", "(leader = 'none')", replacement)
    ggchord_stop("theme(): '", name, "' was removed; use ", replacement)
  }
  unknown <- setdiff(keys, names(ggplot2::get_element_tree()))
  if (length(unknown)) ggchord_stop("theme(): unknown ggchord element '", unknown[1L], "'")
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
  if (plasmid) out <- out + theme(ggchord.axis.minor.ticks.length = unit(.8, "mm"), plot.margin = margin(6, 6, 6, 6))
  if (publication) out <- out + theme(plot.margin = margin(4, 4, 4, 4, unit = "mm"))
  attr(out, "ggchord.explicit_theme") <- list()
  out
}

ggchord_check_axis_unit <- function(x, field, nonnegative = FALSE) {
  if (!is.unit(x) || length(x) != 1L || !is.finite(as.numeric(x)) ||
      !unitType(x) %in% c("inches", "cm", "mm", "points", "bigpts", "picas") ||
      (nonnegative && as.numeric(x) < 0)) {
    requirement <- if (nonnegative) "a nonnegative absolute grid::unit()" else "an absolute grid::unit()"
    ggchord_stop(field, " must be ", requirement)
  }
  invisible(x)
}

# Convert an absolute unit to inches without opening a graphics device.
ggchord_unit_inches <- function(x) {
  value <- as.numeric(x)
  switch(unitType(x), inches = value, cm = value / 2.54,
    mm = value / 25.4, points = value / 72.27,
    bigpts = value / 72, picas = value / 6,
    ggchord_stop("ggchord axis units must be absolute physical units"))
}

# Merge only the local sequence guide's decorations into a plotting copy.
ggchord_seq_axis_plot <- function(plot, guide) {
  if (!is.null(guide$theme)) {
    overrides <- setNames(as.list(guide$theme), paste0("ggchord.", names(guide$theme)))
    plot$theme <- plot$theme + do.call(theme, overrides)
  }
  plot
}

ggchord_seq_axis_settings <- function(plot, guide) {
  circular <- isTRUE(plot$coordinates$ggchord_circular)
  distance <- function(field, default) {
    key <- paste0("ggchord.axis.", field)
    value <- plot$theme[[key]]
    if (is.null(value)) return(unit(default, "mm"))
    ggchord_check_axis_unit(value, key)
    value
  }
  orientation <- guide$orientation %||% "auto"
  list(
    hidden = identical(guide$sequences, character()),
    gap = guide$gap %||% unit(if (circular) 0 else .8, "mm"),
    ticks.length = distance("ticks.length", if (circular) 1.5 else .85),
    minor.ticks.length = distance("minor.ticks.length", if (circular) .7 else .4),
    text.offset = distance("text.offset", .8),
    text.offset.explicit = !is.null(plot$theme[["ggchord.axis.text.offset"]]),
    text.orientation = if (is.numeric(guide$angle)) guide$angle else
      if (orientation == "auto") "parallel" else orientation,
    text.check.overlap = guide$check.overlap %||% TRUE
  )
}

#' @noRd
ggchord_theme_element <- function(plot, name) {
  theme <- plot$theme
  if (name %in% c("ggchord.axis.ticks", "ggchord.axis.minor.ticks")) {
    element <- theme[[name]]
    if (inherits(element, "element_line") && is.null(element@linewidth)) {
      widths <- if (isTRUE(plot$coordinates$ggchord_circular)) {
        c(major = .45, minor = .3)
      } else c(major = .25, minor = .2)
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
  guide <- plot$scales$get_scales("seq_position")$ggchord_axis_guide
  axis_plot <- ggchord_seq_axis_plot(plot, guide)
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
    if (identical(lyr$ggchord_type, "track")) {
      style <- ggchord_track_axis_style(plot, lyr$ggchord_params$track_id)
      fields <- c("axis_line_params", "axis_tick_params", "axis_text_params",
        "axis_title_params", "axis_check_overlap")
      for (name in fields)
        if (is.null(lyr$geom_params[[name]])) lyr$geom_params[[name]] <- style[[name]]
    }
    components <- lyr$ggchord_theme_components
    if (!is.null(components)) {
      for (param_name in names(components)) {
        lyr$geom_params[[param_name]] <- ggchord_resolve_theme_component(
          axis_plot, lyr, param_name, unname(components[[param_name]]))
      }
    }
    plot$layers[[i]] <- lyr
  }
  plot
}
