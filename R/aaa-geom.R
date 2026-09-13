# Shared geom helper loaded before the public layer files.

#' Clone a geom and expose selected standard aesthetics under role names
#' @noRd
rename_geom_aes <- function(geom = ggplot2::GeomPolygon, renames) {
  new_geom <- ggplot2::ggproto(
    paste0("GeomChord", sub("^Geom", "", class(geom)[1])), geom
  )

  aes_names <- names(new_geom$default_aes)
  for (old in names(renames)) {
    aes_names[aes_names == old] <- renames[[old]]
  }
  names(new_geom$default_aes) <- aes_names

  old_handle_na <- geom$handle_na
  new_geom$handle_na <- function(self, data, params) {
    for (old in names(renames)) {
      colnames(data)[colnames(data) == renames[[old]]] <- old
    }
    old_handle_na(data, params)
  }
  old_draw_key <- geom$draw_key
  new_geom$draw_key <- function(data, params, size) {
    for (old in names(renames)) {
      colnames(data)[colnames(data) == renames[[old]]] <- old
    }
    old_draw_key(data, params, size)
  }
  new_geom
}

#' Apply component-specific fixed aesthetics before delegating to ggplot2
#' @noRd
ggchord_component_style <- function(data, params, aesthetics) {
  if (!nrow(data)) return(data)
  for (nm in intersect(names(params), aesthetics)) {
    value <- params[[nm]]
    if (!is.null(value) && length(value) == 1L) data[[nm]] <- value
  }
  data
}

#' Draw genomic-axis components as one standard ggplot2 layer
#' @noRd
GeomChordAxis <- ggplot2::ggproto(
  "GeomChordAxis", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(
    xend = NA_real_, yend = NA_real_, label = NA_character_,
    .component = NA_character_, .draw_as_arc = FALSE,
    colour = "#737A80", alpha = NA, linewidth = 0.3,
    linetype = 1, size = 3, angle = 0, hjust = 0.5, vjust = 0.5,
    family = "", fontface = 1, lineheight = 1.2
  ),
  extra_params = c(
    "na.rm", "line_params", "tick_params", "minor_tick_params", "text_params"
  ),
  draw_key = function(data, params, size) grid::nullGrob(),
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        line_params = list(), tick_params = list(),
                        minor_tick_params = list(),
                        text_params = list()) {
    line <- data[data$.component %in% "line", , drop = FALSE]
    tick <- data[data$.component %in% "major_tick", , drop = FALSE]
    minor_tick <- data[data$.component %in% "minor_tick", , drop = FALSE]
    text <- data[data$.component %in% "text", , drop = FALSE]
    line <- ggchord_component_style(
      line, line_params, c("colour", "alpha", "linewidth", "linetype")
    )
    tick <- ggchord_component_style(
      tick, tick_params, c("colour", "alpha", "linewidth", "linetype")
    )
    if (nrow(tick) && "is_origin" %in% names(tick)) {
      origin <- !is.na(tick$is_origin) & tick$is_origin
      tick$linewidth[origin] <- tick$linewidth[origin] * 1.8
    }
    minor_tick <- ggchord_component_style(
      minor_tick, minor_tick_params,
      c("colour", "alpha", "linewidth", "linetype")
    )
    text <- ggchord_component_style(
      text, text_params,
      c("colour", "alpha", "size", "family", "fontface", "lineheight")
    )
    grobs <- list()
    if (nrow(line)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomPath$draw_panel(
        line, panel_params, coord,
        lineend = line_params$lineend %||% "butt",
        linejoin = line_params$linejoin %||% "round",
        linemitre = 10, na.rm = na.rm
      )
    }
    if (nrow(tick)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomSegment$draw_panel(
        tick, panel_params, coord,
        arrow = tick_params$arrow %||% NULL,
        arrow.fill = tick_params$arrow.fill %||% NULL,
        lineend = tick_params$lineend %||% "butt",
        linejoin = tick_params$linejoin %||% "round",
        na.rm = na.rm
      )
    }
    if (nrow(minor_tick)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomSegment$draw_panel(
        minor_tick, panel_params, coord,
        arrow = minor_tick_params$arrow %||% NULL,
        arrow.fill = minor_tick_params$arrow.fill %||% NULL,
        lineend = minor_tick_params$lineend %||% "butt",
        linejoin = minor_tick_params$linejoin %||% "round",
        na.rm = na.rm
      )
    }
    if (nrow(text)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
        text, panel_params, coord,
        parse = text_params$parse %||% FALSE,
        check_overlap = text_params$check_overlap %||% FALSE,
        na.rm = na.rm
      )
    }
    do.call(grid::grobTree, grobs)
  }
)

#' Draw automatic label text and leaders as one standard ggplot2 layer
#' @noRd
GeomChordGeneLabelRepel <- ggplot2::ggproto(
  "GeomChordGeneLabelRepel", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(
    xend = NA_real_, yend = NA_real_, label = NA_character_,
    .component = NA_character_,
    colour = "#2B3035", alpha = NA, linewidth = 0.25, linetype = 1,
    size = 2.5, angle = 0, hjust = 0.5, vjust = 0.5,
    family = "", fontface = 1, lineheight = 1.2
  ),
  extra_params = c("na.rm", "segment_params", "text_params"),
  draw_key = function(data, params, size) grid::nullGrob(),
  draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                        segment_params = list(), text_params = list()) {
    segment <- data[data$.component %in% "segment", , drop = FALSE]
    text <- data[data$.component %in% "text", , drop = FALSE]
    arc_text <- data[data$.component %in% "arc_text", , drop = FALSE]
    segment <- ggchord_component_style(
      segment, segment_params,
      c("colour", "alpha", "linewidth", "linetype")
    )
    text <- ggchord_component_style(
      text, text_params,
      c("colour", "alpha", "size", "family", "fontface", "lineheight")
    )
    arc_text <- ggchord_component_style(
      arc_text, text_params,
      c("colour", "alpha", "size", "family", "fontface", "lineheight")
    )
    grobs <- list()
    if (nrow(segment)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomSegment$draw_panel(
        segment, panel_params, coord,
        arrow = segment_params$arrow %||% NULL,
        arrow.fill = segment_params$arrow.fill %||% NULL,
        lineend = segment_params$lineend %||% "butt",
        linejoin = segment_params$linejoin %||% "round",
        na.rm = na.rm
      )
    }
    if (nrow(text)) {
      external <- if ("feature_label_mode" %in% names(text)) {
        text$feature_label_mode %in% "external"
      } else rep(FALSE, nrow(text))
      curved <- if (".draw_as_arc" %in% names(text)) {
        text$.draw_as_arc %in% TRUE
      } else rep(FALSE, nrow(text))
      plain <- text[!external & !curved, , drop = FALSE]
      callout <- text[external, , drop = FALSE]
      if (nrow(callout)) {
        callout$fill <- if ("feature_label_fill" %in% names(callout)) {
          callout$feature_label_fill
        } else "#F2F3F4"
        callout$colour <- "#353A3E"
        grobs[[length(grobs) + 1L]] <- ggplot2::GeomLabel$draw_panel(
          callout, panel_params, coord,
          label.padding = grid::unit(.12, "lines"),
          label.r = grid::unit(.16, "lines"), na.rm = na.rm
        )
      }
      if (nrow(plain)) {
        grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
          plain, panel_params, coord,
          parse = text_params$parse %||% FALSE,
          check_overlap = text_params$check_overlap %||% FALSE,
          na.rm = na.rm
        )
      }
    }
    if (nrow(arc_text)) {
      grobs[[length(grobs) + 1L]] <- ggplot2::GeomText$draw_panel(
        arc_text, panel_params, coord,
        parse = FALSE, check_overlap = FALSE, na.rm = na.rm
      )
    }
    do.call(grid::grobTree, grobs)
  }
)
