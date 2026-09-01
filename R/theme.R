# theme.R - compact themes for chord diagrams

#' ggchord themes
#'
#' Theme helpers remove Cartesian axes and provide consistent typography and
#' legend spacing for chord diagrams. `theme_ggchord()` uses a restrained
#' white publication canvas; the other variants provide a still quieter
#' canvas, a dark canvas, or tighter print typography.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @param axis_line,axis_ticks Optional line or blank elements for genomic
#'   axes. `NULL` keeps the preset value.
#' @param axis_text,seq_label,gene_label Optional text or blank elements.
#'   `NULL` keeps the preset value.
#' @param gene_label_segment Optional line or blank element for gene leader
#'   lines. `NULL` keeps the preset value.
#'
#' @return A ggplot2 theme object.
#' @export
theme_ggchord <- function(
    base_size = 11,
    base_family = "",
    axis_line = NULL,
    axis_ticks = NULL,
    axis_text = NULL,
    seq_label = NULL,
    gene_label = NULL,
    gene_label_segment = NULL) {
  out <- ggchord_theme_base(base_size, base_family)
  ggchord_theme_override(
    out, axis_line, axis_ticks, axis_text, seq_label, gene_label,
    gene_label_segment
  )
}

#' Construct the default complete theme before user overrides
#' @noRd
ggchord_theme_base <- function(base_size = 11, base_family = "") {
  scale <- base_size / 11
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, size = 14 * scale, face = "bold", colour = "#202428",
        margin = ggplot2::margin(b = 4 * scale)
      ),
      text = ggplot2::element_text(colour = "#30353A"),
      # A transparent legend box avoids masking labels that legitimately use
      # the plot margin when a compact export leaves little panel space.
      legend.background = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key = ggplot2::element_rect(fill = NA, colour = NA),
      legend.box.spacing = grid::unit(4 * scale, "mm"),
      legend.spacing = grid::unit(2 * scale, "mm"),
      legend.spacing.x = grid::unit(2 * scale, "mm"),
      legend.spacing.y = grid::unit(4 * scale, "mm"),
      legend.text = ggplot2::element_text(size = 8 * scale, colour = "#343A40"),
      legend.title = ggplot2::element_text(
        size = 9 * scale, face = "bold", colour = "#252A2E"
      ),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(5.5, 5.5, 5.5, 5.5),
      ggchord.axis.line = ggplot2::element_line(
        colour = "#737A80", linewidth = 0.3
      ),
      ggchord.axis.ticks = ggplot2::element_line(
        colour = "#8A9096", linewidth = 0.25
      ),
      ggchord.axis.text = ggplot2::element_text(
        colour = "#3D4348", size = 2.8 * ggplot2::.pt * scale
      ),
      ggchord.seq.label = ggplot2::element_text(
        colour = "#252A2E", size = 3 * ggplot2::.pt * scale
      ),
      ggchord.gene.label = ggplot2::element_text(
        colour = "#2B3035", size = 2.4 * ggplot2::.pt * scale
      ),
      ggchord.gene.label.segment = ggplot2::element_line(
        colour = "#858B91", linewidth = 0.25
      )
    )
}

#' @rdname theme_ggchord
#' @export
theme_ggchord_minimal <- function(
    base_size = 11,
    base_family = "",
    axis_line = NULL,
    axis_ticks = NULL,
    axis_text = NULL,
    seq_label = NULL,
    gene_label = NULL,
    gene_label_segment = NULL) {
  out <- ggchord_theme_base(base_size, base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size * 1.25),
      legend.box.spacing = grid::unit(3, "mm"),
      ggchord.axis.line = ggplot2::element_line(
        colour = "#90969B", linewidth = 0.25
      ),
      ggchord.axis.ticks = ggplot2::element_line(
        colour = "#A0A5AA", linewidth = 0.2
      ),
      ggchord.gene.label.segment = ggplot2::element_line(
        colour = "#A1A6AA", linewidth = 0.2
      )
    )
  ggchord_theme_override(
    out, axis_line, axis_ticks, axis_text, seq_label, gene_label,
    gene_label_segment
  )
}

#' @rdname theme_ggchord
#' @export
theme_ggchord_dark <- function(
    base_size = 11,
    base_family = "",
    axis_line = NULL,
    axis_ticks = NULL,
    axis_text = NULL,
    seq_label = NULL,
    gene_label = NULL,
    gene_label_segment = NULL) {
  light <- "#F1F3F5"
  out <- ggchord_theme_base(base_size, base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = light),
      plot.background = ggplot2::element_rect(fill = "#17191C", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#17191C", colour = NA),
      legend.background = ggplot2::element_rect(fill = "#17191C", colour = NA),
      plot.title = ggplot2::element_text(colour = light),
      legend.text = ggplot2::element_text(colour = light),
      legend.title = ggplot2::element_text(colour = light),
      ggchord.axis.line = ggplot2::element_line(colour = "#ADB5BD", linewidth = 0.35),
      ggchord.axis.ticks = ggplot2::element_line(colour = "#ADB5BD", linewidth = 0.3),
      ggchord.axis.text = ggplot2::element_text(colour = light),
      ggchord.seq.label = ggplot2::element_text(colour = light),
      ggchord.gene.label = ggplot2::element_text(colour = light),
      ggchord.gene.label.segment = ggplot2::element_line(
        colour = "#CED4DA", linewidth = 0.3
      )
    )
  ggchord_theme_override(
    out, axis_line, axis_ticks, axis_text, seq_label, gene_label,
    gene_label_segment
  )
}

#' @rdname theme_ggchord
#' @export
theme_ggchord_publication <- function(
    base_size = 9,
    base_family = "",
    axis_line = NULL,
    axis_ticks = NULL,
    axis_text = NULL,
    seq_label = NULL,
    gene_label = NULL,
    gene_label_segment = NULL) {
  out <- ggchord_theme_base(base_size, base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_size * 1.25),
      legend.box.spacing = grid::unit(4, "mm"),
      legend.spacing = grid::unit(2, "mm"),
      legend.key.height = grid::unit(4, "mm"),
      legend.key.width = grid::unit(5, "mm")
    )
  ggchord_theme_override(
    out, axis_line, axis_ticks, axis_text, seq_label, gene_label,
    gene_label_segment
  )
}

#' Apply user-supplied ggchord elements to a complete preset
#' @noRd
ggchord_theme_override <- function(
    theme,
    axis_line = NULL,
    axis_ticks = NULL,
    axis_text = NULL,
    seq_label = NULL,
    gene_label = NULL,
    gene_label_segment = NULL) {
  values <- list(
    ggchord.axis.line = axis_line,
    ggchord.axis.ticks = axis_ticks,
    ggchord.axis.text = axis_text,
    ggchord.seq.label = seq_label,
    ggchord.gene.label = gene_label,
    ggchord.gene.label.segment = gene_label_segment
  )
  values <- values[!vapply(values, is.null, logical(1))]
  if (length(values) == 0L) return(theme)
  theme + do.call(ggplot2::theme, values)
}

#' Resolve a ggchord theme element without leaking theme internals elsewhere
#' @noRd
ggchord_theme_element <- function(plot, name) {
  tryCatch(
    ggplot2::calc_element(name, plot$theme),
    error = function(e) NULL
  )
}

#' Resolve a theme text size (points) to a geom text size (millimetres)
#' @noRd
ggchord_theme_text_size <- function(plot, name, fallback) {
  el <- ggchord_theme_element(plot, name)
  if (is.null(el) || inherits(el, "element_blank") ||
      is.null(el@size) || !is.finite(el@size)) {
    return(fallback)
  }
  el@size / ggplot2::.pt
}

#' Resolve a theme text size in points for responsive guide construction
#' @noRd
ggchord_theme_point_size <- function(plot, name, fallback) {
  el <- ggchord_theme_element(plot, name)
  if (is.null(el) || inherits(el, "element_blank") ||
      is.null(el@size) || length(el@size) != 1L ||
      !is.finite(el@size)) {
    return(fallback)
  }
  el@size
}

#' Apply registered theme elements to ggchord annotation layers
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
        if (inherits(el, "element_blank")) {
          current$alpha <- 0
        } else if (inherits(el, "element_line")) {
          values <- list(
            colour = el@colour, linewidth = el@linewidth,
            linetype = el@linetype, lineend = el@lineend,
            linejoin = el@linejoin
          )
          for (nm in names(values)) {
            if (is.null(current[[nm]]) && !is.null(values[[nm]])) {
              current[[nm]] <- values[[nm]]
            }
          }
        } else if (inherits(el, "element_text")) {
          values <- list(
            colour = el@colour, family = el@family, fontface = el@face,
            lineheight = el@lineheight
          )
          for (nm in names(values)) {
            if (is.null(current[[nm]]) && !is.null(values[[nm]])) {
              current[[nm]] <- values[[nm]]
            }
          }
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
      next
    }

    if (inherits(el, "element_line")) {
      values <- list(
        colour = el@colour,
        linewidth = el@linewidth,
        linetype = el@linetype
      )
      geom_values <- list(lineend = el@lineend, linejoin = el@linejoin)
    } else if (inherits(el, "element_text")) {
      values <- list(
        colour = el@colour,
        family = el@family,
        fontface = el@face,
        lineheight = el@lineheight
      )
      geom_values <- list()
    } else {
      next
    }
    for (nm in names(values)) {
      if (is.null(lyr$aes_params[[nm]]) && !is.null(values[[nm]])) {
        lyr$aes_params[[nm]] <- values[[nm]]
      }
    }
    for (nm in names(geom_values)) {
      if (is.null(lyr$geom_params[[nm]]) && !is.null(geom_values[[nm]])) {
        lyr$geom_params[[nm]] <- geom_values[[nm]]
      }
    }
    plot$layers[[i]] <- lyr
  }
  plot
}
