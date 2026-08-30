# theme.R - compact themes for chord diagrams

#' ggchord themes
#'
#' Theme helpers remove Cartesian axes and provide consistent typography and
#' legend spacing for chord diagrams. `theme_ggchord()` preserves the
#' package's established default appearance; the other variants provide a
#' quieter canvas, a dark canvas, or a print-oriented white canvas.
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
        hjust = 0.5, size = 20 * scale, face = "bold",
        margin = margin(b = 5.5 * scale)
      ),
      legend.background = element_blank(),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.box.spacing = unit(10 * scale, "mm"),
      legend.spacing = unit(5 * scale, "mm"),
      legend.text = element_text(size = 8 * scale),
      legend.title = element_text(size = 10 * scale, face = "bold"),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      ggchord.axis.line = element_line(colour = "#4D4D4D", linewidth = 0.35),
      ggchord.axis.ticks = element_line(colour = "#4D4D4D", linewidth = 0.3),
      ggchord.axis.text = element_text(colour = "#2F2F2F", size = 7 * scale),
      ggchord.seq.label = element_text(colour = "#202020", size = 9 * scale),
      ggchord.group.label = element_text(
        colour = "#202020", size = 9 * scale, face = "bold"
      ),
      ggchord.gene.label = element_text(colour = "#202020", size = 7 * scale),
      ggchord.gene.label.segment = element_line(
        colour = "#6B6B6B", linewidth = 0.3
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
      plot.title = element_text(size = base_size * 1.35),
      legend.box.spacing = unit(6, "mm")
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
