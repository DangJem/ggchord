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
#' When pident is absent, the default is a fixed grey fill with no Identity
#' guide. Alignment length and identity are never inferred from endpoints.
#' Alternative mappings and legends are controlled with `aes()`, `scale_*()`
#' and `guides()`.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param ribbon_ctrl_point Optional vector/list of Bezier edge controls.
#'   Two numbers specify one shared control; four numbers specify controls for
#'   the two ribbon edges. A per-row list may contain pairs of coordinate vectors.
#'   NULL derives controls from the sequence endpoints. With shared branches,
#'   controls shape the far ends while the junction follows the shared trunk.
#' @param ribbon_gap Optional numeric/vector controlling spacing between
#'   sequences and ribbon endpoints. With the default
#'   \code{link_avoid = "uniform"}, \code{NULL} uses one obstacle-aware
#'   clearance across each endpoint. Use \code{"smooth"} for clearance that
#'   changes only next to \code{geom_gene()} and \code{geom_feature()} polygons,
#'   or \code{"none"} for fixed spacing 0.035. Text and leader lines are ignored.
#'   Supplying a number disables the automatic rule and uses that exact spacing.
#'   Linear layouts reduce automatic baseline clearance to fit narrow channels;
#'   an explicit gap or feature clearance that cannot fit raises an error.
#' @param link_avoid One of \code{"uniform"} (default), \code{"none"}, or
#'   \code{"smooth"}. Smooth avoidance changes clearance locally; uniform
#'   avoidance uses one clearance across an endpoint.
#' @inheritParams geom_link_line
#' @param fill Optional fixed ribbon fill. `NULL` keeps the default Identity
#'   mapping.
#' @param alpha,colour,linewidth Standard fixed ribbon styles.
#' @param linetype Fixed outline line type. The default `NULL` uses solid
#'   outlines for adjacent displayed sequences and dashed outlines for
#'   non-adjacent sequences. Chord adjacency includes the first/last pair.
#'   This follows the final display order, independently of alignment direction.
#'   Explicit fixed styles or `aes(ribbon_linetype = ...)` (also
#'   `aes(linetype = ...)`) replace the automatic styles. Use a matching scale
#'   for data mappings; automatic styles do not train a linetype scale.
#'   Global line-type mappings are used when `inherit.aes = TRUE`.
#' @param adjacent_only Whether only neighbouring sequences are connected.
#'   `NULL` uses TRUE in [coord_collinear()] and FALSE in [coord_chord()].
#'   Set FALSE to retain cross-track ribbons. Chord neighbours include the
#'   first/last sequence pair.
#' @param pairs Optional data frame with `qaccver` and `saccver` selecting
#'   unordered comparison pairs. With `adjacent_only = NULL`, explicit pairs
#'   take precedence over the coordinate default. Omitted rows are reported
#'   by [export_ggchord_layout()].
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
#'   geom_seq() + geom_link_ribbon()
#' p
#'
#' # Cross-track connections remain available explicitly in linear layouts.
#' ggchord(seq_data_example, ribbon_data_example) + geom_seq() +
#'   geom_link_ribbon(adjacent_only = FALSE) +
#'   coord_collinear()
#' ggchord(seq_data_example, ribbon_data_example) + geom_seq() +
#'   geom_link_ribbon(pairs = data.frame(
#'     qaccver = seq_data_example$accver[1],
#'     saccver = seq_data_example$accver[3])) + coord_collinear()
geom_link_ribbon <- function(mapping = NULL, data = NULL,
                        ribbon_ctrl_point = NULL,
                        ribbon_gap = NULL,
                        link_avoid = c("uniform", "none", "smooth"),
                        link_branch = c("none", "query", "subject"),
                        link_branch_fraction = 0.15,
                        fill = NULL,
                        alpha = 0.42,
                        colour = "#59636D",
                        linewidth = 0.08,
                        linetype = NULL,
                        adjacent_only = NULL,
                        pairs = NULL,
                        position = "identity",
                        show.legend = TRUE,
                        inherit.aes = FALSE,
                        ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  link_avoid <- match.arg(link_avoid)
  link_branch <- match.arg(link_branch)
  if (!is.null(adjacent_only) && (!is.logical(adjacent_only) ||
      length(adjacent_only) != 1L || is.na(adjacent_only))) {
    ggchord_stop("geom_link_ribbon(): adjacent_only must be NULL, TRUE or FALSE")
  }
  if (!is.null(pairs)) ggchord_require_columns(pairs,
    c("qaccver", "saccver"), "geom_link_ribbon(pairs)")
  ggchord_check_branch_fraction(link_branch_fraction)
  mapping <- ggchord_normalize_mapping(mapping)
  dots <- list(...)
  ggchord_reject_retired(dots, "geom_link_ribbon()", c(
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
  alias <- ggchord_colour_alias(colour, dots, "geom_link_ribbon()", sys.call())
  colour <- alias$colour
  dots <- alias$dots
  outline_mapped <- !is.null(mapping) && "ribbon_colour" %in% names(mapping)
  linetype_mapped <- !is.null(mapping) && "ribbon_linetype" %in% names(mapping)
  plain_linetype_mapped <- !is.null(mapping) && "linetype" %in% names(mapping)
  if (linetype_mapped && plain_linetype_mapped) ggchord_stop(
    "geom_link_ribbon(): map only one of linetype and ribbon_linetype"
  )
  automatic_linetype <- is.null(linetype) &&
    !linetype_mapped && !plain_linetype_mapped

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
  if (automatic_linetype) {
    # A per-ribbon geom default must not add categories to another layer's
    # mapped scale. Explicit mappings replace this ordinary AsIs aesthetic.
    mapping_base[["linetype"]] <- ggplot2::aes(linetype = I(linetype_val))[["linetype"]]
  }

  fixed <- list(linewidth = linewidth)
  if (!is.null(fill)) fixed$ribbon_fill <- fill
  if (!outline_mapped) fixed$colour <- colour
  if (!linetype_mapped && !plain_linetype_mapped && !is.null(linetype)) {
    fixed$linetype <- linetype
  }

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
    adjacent_only = adjacent_only, pairs = pairs,
    link_avoid = link_avoid, link_branch = link_branch, link_branch_fraction = link_branch_fraction,
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

# Plot mappings are not available when a layer is constructed. For an
# explicitly inheriting ribbon, resolve its line type before binding data so
# the automatic AsIs default does not mask an inherited ggplot2 mapping.
ggchord_inherit_ribbon_linetype <- function(layer, mapping) {
  if (!identical(layer$ggchord_type, "ribbon") ||
      !isTRUE(layer$inherit.aes) ||
      !is.null(layer$ggchord_params$ribbon_outline_linetype) ||
      any(c("linetype", "ribbon_linetype") %in%
        names(layer$ggchord_input_mapping))) return(layer)
  inherited <- intersect(c("linetype", "ribbon_linetype"), names(mapping))
  if (!length(inherited)) return(layer)
  if (length(inherited) > 1L) ggchord_stop(
    "geom_link_ribbon(): map only one of linetype and ribbon_linetype"
  )
  # Clone the layer so preparing a plot never stores an inherited mapping on
  # the original reusable layer or a later plot with another global mapping.
  layer <- reconstruct_layer(layer, layer$data)
  layer$mapping$linetype <- NULL
  input_mapping <- layer$ggchord_input_mapping %||% ggplot2::aes()
  input_mapping[[inherited]] <- mapping[[inherited]]
  layer$ggchord_input_mapping <- input_mapping
  if (identical(inherited, "ribbon_linetype")) {
    draw_key <- layer$geom$draw_key
    layer$geom <- make_ribbon_geom(
      outline = "ribbon_colour" %in% names(layer$geom$default_aes),
      linetype = TRUE
    )
    layer$geom$draw_key <- draw_key
  }
  layer
}
