# geom-ribbon.R - alignment ribbon layer
# Fetches pre-computed ribbon polygon data from the package environment and renders it with geom_polygon
# Ribbon parameters are specified in this layer and stored for use at print time

# ---------------------------------------------------------------------------
# Internal helper: clone a polygon geom and expose selected aesthetics under
# internal names.  The ribbon layer uses this so its fill / outline / linetype
# scales are independent from the gene layer's plain "fill" scale and from the
# sequence layer's plain "colour" scale.
# ---------------------------------------------------------------------------
rename_geom_aes <- function(geom = GeomPolygon,
                            renames = c(fill = "zfill",
                                        colour = "zoutline",
                                        linetype = "zlinetype")) {
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

# The geom used by the ribbon layer; fill is exposed as "zfill" so that the
# ribbon and gene layers keep independent fill scales.
ribbon_geom <- rename_geom_aes(GeomPolygon, renames = c(fill = "zfill"))

make_ribbon_geom <- function(outline = FALSE, linetype = FALSE) {
  renames <- c(fill = "zfill")
  if (isTRUE(outline)) renames <- c(renames, colour = "zoutline")
  if (isTRUE(linetype)) renames <- c(renames, linetype = "zlinetype")
  rename_geom_aes(GeomPolygon, renames = renames)
}

#' Add an alignment ribbon layer
#'
#' Draws colored ribbons corresponding to alignment results. Color scheme and spacing parameters are specified here.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param ribbon_color_scheme Character. Color scheme \code{"pident"},
#'   \code{"query"}, \code{"subject"} or \code{"single"}, default \code{"pident"}
#' @param ribbon_colors Optional color vector. Ribbon color parameters
#' @param ribbon_color_by Optional character column name. When set, ribbon fill
#'   is mapped to a continuous colourbar for that numeric column instead of
#'   \code{pident} (e.g. \code{"bitscore"}).
#' @param ribbon_color_limits Optional numeric length-2 limits for
#'   \code{ribbon_color_by}.
#' @param ribbon_color_breaks Optional numeric breaks for the \code{ribbon_color_by}
#'   colourbar.
#' @param ribbon_color_name Optional legend title for the \code{ribbon_color_by}
#'   colourbar (defaults to the column name).
#' @param ribbon_alpha Numeric (0-1). Ribbon transparency, default 0.35
#' @param ribbon_alpha_by Optional character column name. When set, alpha is
#'   scaled continuously from that numeric column.
#' @param ribbon_alpha_range Numeric length-2. Alpha range used by
#'   \code{ribbon_alpha_by}, default \code{c(0.15, 0.9)}.
#' @param ribbon_ctrl_point Optional vector/list. Bezier control points, default c(0,0)
#' @param ribbon_gap Optional numeric/vector. Spacing between sequences and ribbons, default 0.15
#' @param alpha Ribbon transparency (overrides ribbon_alpha), defaults to the value used in the layout
#' @param ribbon_outline_color Character. Color of the ribbon outline (border), default "black"
#' @param ribbon_outline_width Numeric. Line width of the ribbon outline, default 0.05
#' @param ribbon_outline_linetype Numeric or character. Line type of the ribbon outline, default 1 (solid); see \code{linetype} in ggplot2 for options
#' @param ribbon_outline_by Optional discrete column name. When set, outline
#'   colour is mapped by that column and \code{ribbon_outline_colors} controls
#'   the palette.
#' @param ribbon_outline_colors Optional named color vector for
#'   \code{ribbon_outline_by}; unnamed vectors are recycled positionally.
#' @param ribbon_linetype_by Optional discrete column name. When set, outline
#'   linetype is mapped by that column and \code{ribbon_linetypes} controls the
#'   values.
#' @param ribbon_linetypes Optional named linetype vector for
#'   \code{ribbon_linetype_by}.
#' @param ribbon_direction Character. How to visually distinguish same- vs
#'   reverse-orientation alignments: \code{"none"}, \code{"alpha"},
#'   \code{"outline"} or \code{"linetype"}.
#' @param ribbon_direction_colors Named color vector with \code{same} and
#'   \code{reverse} entries, used when \code{ribbon_direction = "outline"}.
#' @param ribbon_direction_linetypes Named linetype vector with \code{same} and
#'   \code{reverse} entries, used when \code{ribbon_direction = "linetype"}.
#' @param ribbon_direction_alpha Named numeric vector with \code{same} and
#'   \code{reverse} entries, used when \code{ribbon_direction = "alpha"}.
#' @param show_legend Whether to show the legend, default TRUE
#' @param legend_position Position of this layer's legend (the Identity(%)
#'   colourbar): one of "left", "right", "top", "bottom" or "inside", default
#'   "left". Pass NULL to let the legend follow
#'   \code{theme(legend.position = ...)} together with the other legends.
#' @param legend_key_length Optional length of the Identity(%) colourbar (the
#'   long dimension: its height when placed on the left/right, its width when
#'   placed on the top/bottom). Accepts a grid unit, e.g.
#'   \code{unit(5, "cm")}, or a number interpreted as centimetres. Default
#'   NULL lets the vertical bar fill the available height (and the horizontal
#'   bar default to 4 cm).
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A list of ggplot2 layers
#' @export
#'
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' p <- ggchord(seq_data_example, ribbon_data_example) +
#'   geom_seq() + geom_ribbon()
#' ggplot2::ggplot_build(p)
geom_ribbon <- function(mapping = NULL, data = NULL,
                        ribbon_color_scheme = NULL,
                        ribbon_colors = NULL,
                        ribbon_color_by = NULL,
                        ribbon_color_limits = NULL,
                        ribbon_color_breaks = NULL,
                        ribbon_color_name = NULL,
                        ribbon_alpha = NULL,
                        ribbon_alpha_by = NULL,
                        ribbon_alpha_range = c(0.15, 0.9),
                        ribbon_ctrl_point = NULL,
                        ribbon_gap = NULL,
                        alpha = NULL,
                        ribbon_outline_color = "black",
                        ribbon_outline_width = 0.05,
                        ribbon_outline_linetype = 1,
                        ribbon_outline_by = NULL,
                        ribbon_outline_colors = NULL,
                        ribbon_linetype_by = NULL,
                        ribbon_linetypes = NULL,
                        ribbon_direction = c("none", "alpha", "outline", "linetype"),
                        ribbon_direction_colors = c(same = "black", reverse = "grey50"),
                        ribbon_direction_linetypes = c(same = "solid", reverse = "dashed"),
                        ribbon_direction_alpha = c(same = 1, reverse = 0.45),
                        show_legend = TRUE,
                        legend_position = "left",
                        legend_key_length = NULL,
                        ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  ribbon_alpha <- alpha %||% ribbon_alpha
  ribbon_direction <- match.arg(ribbon_direction)

  # Determine the fill aesthetic early: the real geometry is injected at build
  # time, but the mapping must already use the correct data column.
  scheme <- if (!is.null(ribbon_color_by)) "value" else ribbon_color_scheme %||% "pident"
  outline_mapped <- !is.null(ribbon_outline_by) || identical(ribbon_direction, "outline")
  linetype_mapped <- !is.null(ribbon_linetype_by) || identical(ribbon_direction, "linetype")

  empty_polys <- data.frame(
    x = numeric(0), y = numeric(0),
    group = integer(0),
    pident = numeric(0),
    value = numeric(0),
    fill = character(0),
    alpha = numeric(0),
    outline_col = character(0),
    linetype_val = character(0),
    stringsAsFactors = FALSE
  )

  if (scheme == "pident") {
    fill_mapping <- aes(fill = pident)
  } else if (scheme == "value") {
    fill_mapping <- aes(fill = value)
  } else {
    fill_mapping <- aes(fill = fill)
  }

  mapping_base <- aes(x = x, y = y, group = group, alpha = alpha)
  mapping_base[["fill"]] <- fill_mapping$fill
  if (isTRUE(outline_mapped)) mapping_base[["zoutline"]] <- as.name("outline_col")
  if (isTRUE(linetype_mapped)) mapping_base[["zlinetype"]] <- as.name("linetype_val")

  outline_params <- list()
  if (!outline_mapped) {
    if (!is.null(ribbon_outline_color)) outline_params$colour <- ribbon_outline_color
    if (!is.null(ribbon_outline_linetype)) outline_params$linetype <- ribbon_outline_linetype
  }
  if (!is.null(ribbon_outline_width)) outline_params$linewidth <- ribbon_outline_width

  lyr <- ggplot2::layer(
    data        = empty_polys,
    mapping     = mapping_base,
    stat        = "identity",
    geom        = make_ribbon_geom(outline = outline_mapped, linetype = linetype_mapped),
    position    = "identity",
    show.legend = if (identical(show_legend, TRUE)) {
                    c(fill = TRUE, colour = FALSE, linetype = FALSE)
                  } else show_legend,
    inherit.aes = FALSE,
    check.aes   = FALSE,
    check.param = FALSE,
    key_glyph   = key_glyph_ribbon,
    params      = c(list(...), outline_params)
  )
  lyr$ggchord_type <- "ribbon"
  lyr$ggchord_params <- list(
    type                      = "ribbon",
    ribbon_color_scheme       = ribbon_color_scheme,
    ribbon_colors             = ribbon_colors,
    ribbon_color_by           = ribbon_color_by,
    ribbon_color_limits       = ribbon_color_limits,
    ribbon_color_breaks       = ribbon_color_breaks,
    ribbon_color_name         = ribbon_color_name,
    ribbon_alpha              = ribbon_alpha,
    ribbon_alpha_by           = ribbon_alpha_by,
    ribbon_alpha_range        = ribbon_alpha_range,
    ribbon_ctrl_point         = ribbon_ctrl_point,
    ribbon_gap                = ribbon_gap,
    ribbon_outline_color      = ribbon_outline_color,
    ribbon_outline_width      = ribbon_outline_width,
    ribbon_outline_linetype   = ribbon_outline_linetype,
    ribbon_outline_by         = ribbon_outline_by,
    ribbon_outline_colors     = ribbon_outline_colors,
    ribbon_linetype_by        = ribbon_linetype_by,
    ribbon_linetypes          = ribbon_linetypes,
    ribbon_direction          = ribbon_direction,
    ribbon_direction_colors   = ribbon_direction_colors,
    ribbon_direction_linetypes = ribbon_direction_linetypes,
    ribbon_direction_alpha    = ribbon_direction_alpha,
    legend_position           = legend_position,
    legend_key_length         = legend_key_length
  )
  list(lyr)
}
