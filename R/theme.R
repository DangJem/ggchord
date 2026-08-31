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
#'
#' @return A ggplot2 theme object.
#' @export
theme_ggchord <- function(base_size = 11, base_family = "") {
  scale <- base_size / 11
  theme_void(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        hjust = 0.5, size = 14 * scale, face = "bold", colour = "#202428",
        margin = margin(b = 4 * scale)
      ),
      text = element_text(colour = "#30353A"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.box.spacing = unit(4 * scale, "mm"),
      legend.spacing = unit(2 * scale, "mm"),
      legend.text = element_text(size = 8 * scale, colour = "#343A40"),
      legend.title = element_text(
        size = 9 * scale, face = "bold", colour = "#252A2E"
      ),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5),
      ggchord.axis.line = element_line(colour = "#737A80", linewidth = 0.3),
      ggchord.axis.ticks = element_line(colour = "#8A9096", linewidth = 0.25),
      ggchord.axis.text = element_text(
        colour = "#3D4348", size = 2.8 * ggplot2::.pt * scale
      ),
      ggchord.seq.label = element_text(
        colour = "#252A2E", size = 3 * ggplot2::.pt * scale
      ),
      ggchord.group.label = element_text(
        colour = "#202020", size = 3.5 * ggplot2::.pt * scale, face = "bold"
      ),
      ggchord.gene.label = element_text(
        colour = "#2B3035", size = 2.4 * ggplot2::.pt * scale
      ),
      ggchord.gene.label.segment = element_line(
        colour = "#858B91", linewidth = 0.25
      )
    )
}

#' @rdname theme_ggchord
#' @export
theme_ggchord_minimal <- function(base_size = 11, base_family = "") {
  theme_ggchord(base_size, base_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = base_size * 1.25),
      legend.box.spacing = unit(3, "mm")
    )
}

#' @rdname theme_ggchord
#' @export
theme_ggchord_dark <- function(base_size = 11, base_family = "") {
  light <- "#F1F3F5"
  theme_ggchord(base_size, base_family) +
    theme(
      text = element_text(colour = light),
      plot.background = element_rect(fill = "#17191C", colour = NA),
      panel.background = element_rect(fill = "#17191C", colour = NA),
      legend.background = element_rect(fill = "#17191C", colour = NA),
      plot.title = element_text(colour = light),
      legend.text = element_text(colour = light),
      legend.title = element_text(colour = light),
      ggchord.axis.line = element_line(colour = "#ADB5BD", linewidth = 0.35),
      ggchord.axis.ticks = element_line(colour = "#ADB5BD", linewidth = 0.3),
      ggchord.axis.text = element_text(colour = light),
      ggchord.seq.label = element_text(colour = light),
      ggchord.group.label = element_text(colour = light, face = "bold"),
      ggchord.gene.label = element_text(colour = light),
      ggchord.gene.label.segment = element_line(
        colour = "#CED4DA", linewidth = 0.3
      )
    )
}

#' @rdname theme_ggchord
#' @export
theme_ggchord_publication <- function(base_size = 9, base_family = "") {
  theme_ggchord(base_size, base_family) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = base_size * 1.25),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5),
      legend.box.spacing = unit(4, "mm"),
      legend.spacing = unit(2, "mm"),
      legend.key.height = unit(4, "mm"),
      legend.key.width = unit(5, "mm")
    )
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

#' Apply registered theme elements to ggchord annotation layers
#' @noRd
ggchord_apply_theme_styles <- function(plot) {
  for (i in seq_along(plot$layers)) {
    lyr <- plot$layers[[i]]
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
    if (identical(name, "ggchord.group.label")) values$colour <- NULL
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
