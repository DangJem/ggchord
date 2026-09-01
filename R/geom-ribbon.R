# geom-ribbon.R - alignment ribbon layer
# Fetches pre-computed ribbon polygon data from the package environment and renders it with geom_polygon
# Ribbon parameters are specified in this layer and stored for use at print time

# ---------------------------------------------------------------------------
# Internal helper: clone a polygon geom and expose selected aesthetics under
# internal names.  The ribbon layer uses this so its fill / outline / linetype
# scales are independent from the gene layer's plain "fill" scale and from the
# sequence layer's plain "colour" scale.
# ---------------------------------------------------------------------------
# The geom used by the ribbon layer; fill is exposed as "zfill" so that the
# ribbon and gene layers keep independent fill scales.
ribbon_geom <- rename_geom_aes(
  ggplot2::GeomPolygon, renames = c(fill = "ribbon_fill", alpha = "ribbon_alpha")
)

make_ribbon_geom <- function(outline = FALSE, linetype = FALSE) {
  renames <- c(fill = "ribbon_fill", alpha = "ribbon_alpha")
  if (isTRUE(outline)) renames <- c(renames, colour = "ribbon_colour")
  if (isTRUE(linetype)) renames <- c(renames, linetype = "ribbon_linetype")
  rename_geom_aes(ggplot2::GeomPolygon, renames = renames)
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
#' @param ribbon_alpha Numeric (0-1). Ribbon transparency, default 0.42
#' @param ribbon_alpha_by Optional character column name. When set, alpha is
#'   scaled continuously from that numeric column.
#' @param ribbon_alpha_range Numeric length-2. Alpha range used by
#'   \code{ribbon_alpha_by}, default \code{c(0.15, 0.9)}.
#' @param ribbon_ctrl_point Optional vector/list. Bezier control points, default c(0,0)
#' @param ribbon_gap Optional numeric/vector controlling spacing between
#'   sequences and ribbon endpoints. The default \code{NULL} uses local
#'   obstacle-aware spacing: endpoints move closer to \code{geom_seq()} where
#'   no \code{geom_gene()} or \code{geom_feature()} polygon overlaps that
#'   genomic interval, and retain enough clearance where one does. Text and
#'   leader lines are ignored. Supplying a number disables the automatic rule
#'   and uses that exact spacing.
#' @param alpha Ribbon transparency (overrides ribbon_alpha), defaults to the value used in the layout
#' @param ribbon_outline_color Character. Colour of the ribbon outline,
#'   default \code{"#59636D"}, a restrained dark neutral that remains visible
#'   after ribbon transparency is applied.
#' @param ribbon_outline_width Numeric. Line width of the ribbon outline, default 0.08
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
#' @param legend_position Position of this layer's legend (the Identity (%)
#'   colourbar): one of "left", "right", "top", "bottom" or "inside", default
#'   "left". Pass NULL to let the legend follow
#'   \code{theme(legend.position = ...)} together with the other legends.
#' @param legend_key_width Optional width of the Identity (%) colourbar key.
#'   Accepts a grid unit, e.g. \code{unit(1, "cm")}, or a number interpreted
#'   as centimetres. Default NULL uses the package's 3.6 mm vertical-bar width
#'   or 50 mm horizontal-bar length.
#' @param legend_key_height Optional height of the Identity (%) colourbar key.
#'   Accepts a grid unit, e.g. \code{unit(5, "cm")}, or a number interpreted
#'   as centimetres. Default NULL uses the package's 50 mm vertical-bar length
#'   or 3.6 mm horizontal-bar height so the guide remains stable across
#'   output-device sizes.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A list of ggplot2 layers
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' p <- ggchord(seq_data_example, ribbon_data_example) +
#'   geom_seq() + geom_ribbon()
#' p
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
                        ribbon_outline_color = "#59636D",
                        ribbon_outline_width = 0.08,
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
                        legend_key_width = NULL,
                        legend_key_height = NULL,
                        ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!missing(legend_position) || !missing(legend_key_width) ||
      !missing(legend_key_height)) {
    ggchord_deprecate_once(
      "geom_ribbon(legend_position/legend_key_width/legend_key_height)",
      "guides(ribbon_fill = guide_ggchord_colourbar(...))"
    )
  }

  ribbon_alpha <- alpha %||% ribbon_alpha
  ribbon_direction <- match.arg(ribbon_direction)

  # Determine the fill aesthetic early: the real geometry is injected at build
  # time, but the mapping must already use the correct data column.
  scheme <- if (!is.null(ribbon_color_by)) "value" else ribbon_color_scheme %||% "pident"
  outline_mapped <- !is.null(ribbon_outline_by) ||
    identical(ribbon_direction, "outline") ||
    (!is.null(mapping) && "ribbon_colour" %in% names(mapping))
  linetype_mapped <- !is.null(ribbon_linetype_by) ||
    identical(ribbon_direction, "linetype") ||
    (!is.null(mapping) && "ribbon_linetype" %in% names(mapping))

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
    fill_mapping <- ggplot2::aes(fill = pident)
  } else if (scheme == "value") {
    fill_mapping <- ggplot2::aes(fill = value)
  } else {
    fill_mapping <- ggplot2::aes(fill = fill)
  }

  mapping_base <- ggplot2::aes(x = x, y = y, group = group, ribbon_alpha = alpha)
  mapping_base[["ribbon_fill"]] <- fill_mapping$fill
  if (isTRUE(outline_mapped)) mapping_base[["ribbon_colour"]] <- as.name("outline_col")
  if (isTRUE(linetype_mapped)) mapping_base[["ribbon_linetype"]] <- as.name("linetype_val")

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
                    c(ribbon_fill = TRUE, ribbon_colour = FALSE,
                      ribbon_linetype = FALSE)
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
    legend_key_width          = legend_key_width,
    legend_key_height         = legend_key_height
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping,
    c("qaccver", "saccver", "length", "pident", "qstart", "qend",
      "sstart", "send")
  )
  legacy <- list(
    c("ribbon_color_scheme", "ribbon_fill", "ggplot2::aes(ribbon_fill = ...)",
      !missing(ribbon_color_scheme) && !is.null(ribbon_color_scheme)),
    c("ribbon_colors", "ribbon_fill", "scale_ribbon_fill_*()",
      !missing(ribbon_colors) && !is.null(ribbon_colors)),
    c("ribbon_color_by", "ribbon_fill", "ggplot2::aes(ribbon_fill = ...)",
      !missing(ribbon_color_by) && !is.null(ribbon_color_by)),
    c("ribbon_color_limits", "ribbon_fill", "scale_ribbon_fill_*(limits = ...)",
      !missing(ribbon_color_limits) && !is.null(ribbon_color_limits)),
    c("ribbon_color_breaks", "ribbon_fill", "scale_ribbon_fill_*(breaks = ...)",
      !missing(ribbon_color_breaks) && !is.null(ribbon_color_breaks)),
    c("ribbon_color_name", "ribbon_fill", "scale_ribbon_fill_*(name = ...)",
      !missing(ribbon_color_name) && !is.null(ribbon_color_name)),
    c("ribbon_alpha_by", "ribbon_alpha", "ggplot2::aes(ribbon_alpha = ...)",
      !missing(ribbon_alpha_by) && !is.null(ribbon_alpha_by)),
    c("ribbon_alpha_range", "ribbon_alpha", "scale_ribbon_alpha_continuous(range = ...)",
      !missing(ribbon_alpha_range)),
    c("ribbon_outline_by", "ribbon_colour", "ggplot2::aes(ribbon_colour = ...)",
      !missing(ribbon_outline_by) && !is.null(ribbon_outline_by)),
    c("ribbon_outline_colors", "ribbon_colour", "scale_ribbon_colour_manual(values = ...)",
      !missing(ribbon_outline_colors) && !is.null(ribbon_outline_colors)),
    c("ribbon_linetype_by", "ribbon_linetype", "ggplot2::aes(ribbon_linetype = ...)",
      !missing(ribbon_linetype_by) && !is.null(ribbon_linetype_by)),
    c("ribbon_linetypes", "ribbon_linetype", "scale_ribbon_linetype_manual(values = ...)",
      !missing(ribbon_linetypes) && !is.null(ribbon_linetypes))
  )
  for (spec in legacy) {
    lyr <- ggchord_add_legacy_scale(
      lyr, identical(spec[[4]], "TRUE"), spec[[1]], spec[[2]], spec[[3]]
    )
  }
  if (!missing(ribbon_direction) && !identical(ribbon_direction, "none")) {
    direction_aesthetic <- switch(
      ribbon_direction, alpha = "ribbon_alpha", outline = "ribbon_colour",
      linetype = "ribbon_linetype", "ribbon_fill"
    )
    lyr <- ggchord_add_legacy_scale(
      lyr, TRUE, "ribbon_direction", direction_aesthetic,
      paste0("ggplot2::aes(", direction_aesthetic, " = after_stat(direction))")
    )
  }
  lyr <- ggchord_add_legacy_scale(
    lyr, !missing(ribbon_direction_colors), "ribbon_direction_colors",
    "ribbon_colour", "scale_ribbon_colour_manual(values = ...)"
  )
  lyr <- ggchord_add_legacy_scale(
    lyr, !missing(ribbon_direction_linetypes), "ribbon_direction_linetypes",
    "ribbon_linetype", "scale_ribbon_linetype_manual(values = ...)"
  )
  lyr <- ggchord_add_legacy_scale(
    lyr, !missing(ribbon_direction_alpha), "ribbon_direction_alpha",
    "ribbon_alpha", "scale_ribbon_alpha_manual(values = ...)"
  )
  list(lyr)
}
