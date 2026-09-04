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
#' Draws alignment ribbons. Identity is mapped to `ribbon_fill` by default;
#' alternative mappings and legends are controlled with `aes()`, `scale_*()`
#' and `guides()`.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param ribbon_ctrl_point Optional vector/list. Bezier control points, default c(0,0)
#' @param ribbon_gap Optional numeric/vector controlling spacing between
#'   sequences and ribbon endpoints. The default \code{NULL} uses local
#'   obstacle-aware spacing: endpoints move closer to \code{geom_seq()} where
#'   no \code{geom_gene()} or \code{geom_feature()} polygon overlaps that
#'   genomic interval, and retain enough clearance where one does. Text and
#'   leader lines are ignored. Supplying a number disables the automatic rule
#'   and uses that exact spacing.
#' @param fill Optional fixed ribbon fill. `NULL` keeps the default Identity
#'   mapping.
#' @param alpha,colour,linewidth,linetype Standard fixed ribbon styles.
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to \code{geom_polygon()}
#'
#' @return A ggplot2 layer.
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
                        ribbon_ctrl_point = NULL,
                        ribbon_gap = NULL,
                        fill = NULL,
                        alpha = 0.42,
                        colour = "#59636D",
                        linewidth = 0.08,
                        linetype = 1,
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE,
                        ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  mapping <- ggchord_normalize_mapping(mapping)
  dots <- list(...)
  ggchord_reject_retired(dots, "geom_ribbon()", c(
    ribbon_color_scheme = "aes(ribbon_fill = ...) and scale_ribbon_fill_*()",
    ribbon_colors = "scale_ribbon_fill_*()",
    ribbon_color_by = "aes(ribbon_fill = ...)",
    ribbon_color_limits = "scale_ribbon_fill_*(limits = ...)",
    ribbon_color_breaks = "scale_ribbon_fill_*(breaks = ...)",
    ribbon_color_name = "scale_ribbon_fill_*(name = ...)",
    ribbon_alpha = "alpha",
    ribbon_alpha_by = "aes(ribbon_alpha = ...)",
    ribbon_alpha_range = "scale_ribbon_alpha_continuous(range = ...)",
    ribbon_outline_color = "colour",
    ribbon_outline_width = "linewidth",
    ribbon_outline_linetype = "linetype",
    ribbon_outline_by = "aes(ribbon_colour = ...)",
    ribbon_outline_colors = "scale_ribbon_colour_manual(values = ...)",
    ribbon_linetype_by = "aes(ribbon_linetype = ...)",
    ribbon_linetypes = "scale_ribbon_linetype_manual(values = ...)",
    ribbon_direction = "aes() with a direction column and a matching scale",
    ribbon_direction_colors = "scale_ribbon_colour_manual(values = ...)",
    ribbon_direction_linetypes = "scale_ribbon_linetype_manual(values = ...)",
    ribbon_direction_alpha = "scale_ribbon_alpha_manual(values = ...)",
    show_legend = "show.legend",
    legend_position = "guides(ribbon_fill = guide_ggchord_colourbar(position = ...))",
    legend_key_width = "guide_ggchord_colourbar(barwidth = ...)",
    legend_key_height = "guide_ggchord_colourbar(barheight = ...)"
  ))
  alias <- ggchord_colour_alias(colour, dots, "geom_ribbon()", sys.call())
  colour <- alias$colour
  dots <- alias$dots
  outline_mapped <- !is.null(mapping) && "ribbon_colour" %in% names(mapping)
  linetype_mapped <- !is.null(mapping) && "ribbon_linetype" %in% names(mapping)

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

  mapping_base <- ggplot2::aes(x = x, y = y, group = group, ribbon_alpha = alpha)
  if (is.null(fill)) mapping_base[["ribbon_fill"]] <- as.name("pident")
  if (isTRUE(outline_mapped)) mapping_base[["ribbon_colour"]] <- as.name("outline_col")
  if (isTRUE(linetype_mapped)) mapping_base[["ribbon_linetype"]] <- as.name("linetype_val")

  fixed <- list(linewidth = linewidth)
  if (!is.null(fill)) fixed$ribbon_fill <- fill
  if (!outline_mapped) fixed$colour <- colour
  if (!linetype_mapped) fixed$linetype <- linetype

  lyr <- ggplot2::layer(
    data        = empty_polys,
    mapping     = mapping_base,
    stat        = "identity",
    geom        = make_ribbon_geom(outline = outline_mapped, linetype = linetype_mapped),
    position    = position,
    show.legend = if (identical(show.legend, TRUE)) {
                    c(ribbon_fill = TRUE, ribbon_colour = FALSE,
                      ribbon_linetype = FALSE)
                  } else show.legend,
    inherit.aes = inherit.aes,
    check.aes   = FALSE,
    check.param = FALSE,
    key_glyph   = key_glyph_ribbon,
    params      = c(dots, fixed)
  )
  lyr$ggchord_type <- "ribbon"
  lyr$ggchord_params <- list(
    type                      = "ribbon",
    ribbon_color_scheme       = "pident",
    ribbon_colors             = NULL,
    ribbon_color_by           = NULL,
    ribbon_color_limits       = NULL,
    ribbon_color_breaks       = NULL,
    ribbon_color_name         = NULL,
    ribbon_alpha              = alpha,
    ribbon_alpha_by           = NULL,
    ribbon_alpha_range        = c(0.15, 0.9),
    ribbon_ctrl_point         = ribbon_ctrl_point,
    ribbon_gap                = ribbon_gap,
    ribbon_outline_color      = colour,
    ribbon_outline_width      = linewidth,
    ribbon_outline_linetype   = linetype,
    ribbon_outline_by         = NULL,
    ribbon_outline_colors     = NULL,
    ribbon_linetype_by        = NULL,
    ribbon_linetypes          = NULL,
    ribbon_direction          = "none",
    ribbon_direction_colors   = NULL,
    ribbon_direction_linetypes = NULL,
    ribbon_direction_alpha    = NULL,
    legend_position           = NULL,
    legend_key_width          = NULL,
    legend_key_height         = NULL
  )
  lyr <- ggchord_capture_layer_input(
    lyr, data, mapping,
    c("qaccver", "saccver", "length", "pident", "qstart", "qend",
      "sstart", "send")
  )
  lyr
}
