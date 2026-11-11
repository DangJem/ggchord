#' Draw HSP alignment segments
#'
#' Projects each standard ribbon-data row from `(qstart, sstart)` to
#' `(qend, send)` in [coord_dotplot()]. All four endpoints are retained;
#' reverse alignments must have normalized descending endpoints, as returned
#' by [read_paf()]. Explicit unknown `direction` or `strand` values are omitted
#' with a warning and recorded in the layout connection report.
#' This is an HSP-level view, not a reconstruction of mismatches or gaps.
#'
#' @param mapping,data Standard ggplot2 mapping and layer data. Data defaults
#'   to the constructor's `ribbon_data`. Geometry roles are `qaccver`,
#'   `saccver`, `qstart`, `qend`, `sstart`, and `send`.
#' @param colour,alpha,linewidth,linetype Fixed styles. Mapped styles use
#'   `link_colour`, `link_alpha`, `link_linewidth` and `link_linetype`, with
#'   the corresponding [scale_link] functions. Explicit fixed styles win.
#' @param allow_repeated Allow explicitly identified `source_hsp_id` records to
#'   be displayed more than once. All participating layers must opt in. Rows
#'   without declared IDs are never deduplicated by their values.
#' @param position Only `"identity"` is supported: data adjustment would
#'   change the alignment endpoints.
#' @param show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional arguments passed to ggplot2's segment geometry.
#' @return A ggplot2 layer for a ggchord Dotplot.
#' @export
#' @examples
#' data(seq_data_example)
#' data(ribbon_data_example)
#' ggchord(seq_data_example, ribbon_data_example) +
#'   geom_alignment_segment(aes(link_colour = pident)) +
#'   coord_dotplot()
geom_alignment_segment <- function(mapping = NULL, data = NULL,
    colour = "#59636D", alpha = 0.6, linewidth = 0.4, linetype = 1,
    position = "identity", show.legend = NA, inherit.aes = FALSE,
    allow_repeated = FALSE, ...) {
  if (!is.logical(allow_repeated) || length(allow_repeated) != 1L || is.na(allow_repeated))
    ggchord_stop("allow_repeated must be TRUE or FALSE")
  if (!identical(position, "identity") && !inherits(position, "PositionIdentity"))
    ggchord_stop("geom_alignment_segment(): only position = 'identity' is supported")
  defaults <- c(link_colour = missing(colour), link_alpha = missing(alpha),
    link_linewidth = missing(linewidth), link_linetype = missing(linetype))
  alias <- ggchord_colour_alias(colour, list(...), "geom_alignment_segment()", sys.call())
  if ("color" %in% names(list(...))) defaults[["link_colour"]] <- FALSE
  layer <- ggplot2::layer(
    data = data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric()),
    mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    geom = GeomAlignmentSegment, stat = "identity", position = position,
    show.legend = show.legend, inherit.aes = inherit.aes, check.aes = FALSE,
    params = c(list(link_colour = alias$colour, link_alpha = alpha,
      link_linewidth = linewidth, link_linetype = linetype), alias$dots))
  layer$ggchord_type <- "alignment"
  layer$ggchord_params <- list(type = "alignment", style_defaults = defaults, allow_repeated = allow_repeated)
  ggchord_capture_layer_input(layer, data, mapping,
    c("qaccver", "saccver", "qstart", "qend", "sstart", "send"))
}

GeomAlignmentSegment <- ggplot2::ggproto("GeomAlignmentSegment", ggplot2::GeomSegment,
  rename_size = FALSE,
  default_aes = ggplot2::aes(link_colour = "#59636D", link_alpha = 0.6,
    link_linewidth = 0.4, link_linetype = 1),
  draw_panel = function(data, panel_params, coord, arrow = NULL,
      arrow.fill = NULL, lineend = "butt", linejoin = "round", na.rm = FALSE) {
    for (name in c("colour", "alpha", "linewidth", "linetype"))
      data[[name]] <- data[[paste0("link_", name)]]
    ggplot2::GeomSegment$draw_panel(data, panel_params, coord,
      arrow = arrow, arrow.fill = arrow.fill, lineend = lineend,
      linejoin = linejoin, na.rm = na.rm)
  },
  draw_key = function(data, params, size) {
    for (name in c("colour", "alpha", "linewidth", "linetype"))
      data[[name]] <- data[[paste0("link_", name)]]
    ggplot2::draw_key_path(data, params, size)
  })
