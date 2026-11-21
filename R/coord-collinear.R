# Linear multi-sequence coordinates own their physical layout strategy.
CoordCollinear <- ggproto("CoordCollinear", CoordGgchord)

#' Project references onto parallel linear paths
#'
#' @inheritParams coord_dotplot
#' @param axis Direction of the reference paths, either `"x"` or `"y"`.
#' @return A coordinate system for a ggchord plot.
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
#' links <- data.frame(
#'   qaccver = "A", saccver = "B", qstart = c(100, 400),
#'   qend = c(240, 620), sstart = c(150, 450), send = c(290, 650),
#'   qstrand = "+", sstrand = c("+", "-"),
#'   pident = c(98, 90), link_id = c("l1", "l2")
#' )
#'
#' ggchord(data = seqs) +
#'   geom_seq() +
#'   geom_gene(data = genes) +
#'   geom_link_ribbon(data = links) +
#'   layout_sequences(orientation = c(A = 1, B = -1), row_gap = 8) +
#'   coord_collinear(axis = "x")
coord_collinear <- function(axis = "x", ratio = NULL, xlim = NULL, ylim = NULL,
    expand = FALSE, clip = "on", fit = "geometry") {
  ggchord_check_argument_names(sys.call(), sys.function())
  axis <- ggchord_match_choice(axis, c("x", "y"), "axis")
  settings <- ggchord_coord_settings(ratio, xlim, ylim, expand, clip, fit)
  ggproto(NULL, CoordCollinear, !!!settings,
    ggchord_collinear = TRUE, rotation = 0,
    collinear_align = "left", collinear_orientation = "manual",
    collinear_seq_axis = if (axis == "x") "horizontal" else "vertical")
}

ggchord_validate_coord_contract <- function(plot) {
  coord <- plot$coordinates
  if (!isTRUE(coord$ggchord_coord) && !isTRUE(coord$ggchord_dotplot))
    ggchord_stop("ggchord(): choose a sequence coordinate system")
  if (any(unlist(plot$facet$params$free)) &&
      (!isTRUE(coord$ggchord_dotplot) || !is.null(coord$ratio)))
    ggchord_stop("Free facet axes require coord_dotplot(ratio = NULL)")
  if (isTRUE(coord$ggchord_dotplot) && !is.null(plot$scales$get_scales("seq_position")))
    ggchord_stop("scale_seq_position() requires a single reference path")
  value_scale <- plot$scales$get_scales("y")
  has_track <- any(vapply(plot$layers, function(layer)
    identical(layer$ggchord_type, "track") && !identical(layer$ggchord_params$track_kind, "tile"), logical(1)))
  if (has_track && !is.null(value_scale) && !identical(value_scale$get_transformation()$name, "identity"))
    ggchord_stop("Track heights use scale_track_value(transform = ...)")
  invisible(plot)
}
