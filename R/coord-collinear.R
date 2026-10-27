# Linear multi-sequence coordinates own their physical layout strategy.
CoordCollinear <- ggplot2::ggproto("CoordCollinear", CoordGgchord)

#' Linear multi-sequence collinearity coordinates
#'
#' Draws sequences as left-aligned horizontal backbones, ordered from top to
#' bottom by `geom_seq(seq_order = ...)`. All sequences share one linear bp
#' scale. `seq_orientation = -1` reverses a sequence together with its features,
#' links and axis without modifying the input data.
#'
#' In this coordinate system, `geom_seq(seq_gap = ...)` specifies the distance
#' between consecutive backbones as a proportion of the longest backbone's
#' width (default 0.3). Positive Position offsets point upward and negative
#' offsets downward, independently of sequence orientation. Sequence labels
#' appear to the left of their backbone. Automatic feature labels use local
#' horizontal rows; fixed radial labels are vertical and tangent labels are
#' horizontal. Label side `"outside"` means above and `"inside"` below.
#'
#' Explicit sequence radius, curvature, ring mappings/scales and sequence-label
#' radius cannot be converted and produce an error. Remove these arguments
#' when switching from [coord_chord()]. Circular-only primer, restriction-site
#' and centre-label layers require [coord_circular()]; use [geom_feature()] for
#' general linear annotations. Cross-track links may cross intervening tracks.
#' Explicit link control points use the exported linear layout x/y units;
#' the longest backbone has width 2.
#' The circular `external = FALSE` feature-label constraint cannot be converted;
#' remove it to use linear label rows. Gene/feature repel layout modes share
#' the linear row strategy rather than selecting circular label tracks.
#'
#' @inheritParams coord_chord
#' @return A ggplot2 Coord object for composition with a ggchord plot.
#' @md
#' @export
#' @examples
#' seqs <- data.frame(accver = c("A", "B"), length = c(1000, 800))
#' links <- data.frame(qaccver = "A", saccver = "B",
#'                     qstart = 100, qend = 400, sstart = 200, send = 500)
#' ggchord(seqs, links) + geom_seq() + geom_link_ribbon() +
#'   geom_seq_label() + coord_collinear()
coord_collinear <- function(
    ratio = 1, xlim = NULL, ylim = NULL,
    expand = FALSE, clip = "off",
    fit = c("labels", "geometry", "manual")) {
  fit <- match.arg(fit)
  validated <- coord_chord(
    rotation = 0, ratio = ratio, xlim = xlim, ylim = ylim,
    expand = expand, clip = clip, fit = fit
  )
  coord <- new_ggchord_coord(
    CoordCollinear, 0, ratio, xlim, ylim, expand, clip, fit,
    validated$user_xlim, validated$user_ylim
  )
  coord$ggchord_collinear <- TRUE
  coord
}

ggchord_validate_coord_contract <- function(plot) {
  if (!isTRUE(plot$coordinates$ggchord_collinear)) return(invisible(plot))
  for (layer in plot$layers) {
    params <- layer$ggchord_params
    if (is.null(params)) next
    if (identical(params$type, "seq")) {
      for (name in c("seq_radius", "seq_curvature")) {
        if (!is.null(params[[name]])) ggchord_stop(
          "coord_collinear(): explicit `", name, "` cannot be converted; ",
          "remove it from geom_seq() and use seq_gap for track spacing"
        )
      }
      if ("seq_ring" %in% names(layer$ggchord_input_mapping)) ggchord_stop(
        "coord_collinear(): `seq_ring` is circular; remove the mapping and ",
        "scale_seq_ring_manual(), and use seq_order/seq_gap for linear tracks"
      )
    }
    if (isTRUE(params$seq_label_radius_supplied)) ggchord_stop(
      "coord_collinear(): explicit `seq_label_radius` cannot be converted; ",
      "remove it from geom_seq_label() to use labels beside each track"
    )
    if (isTRUE(params$is_feature_label) &&
        identical(params$feature_label_external, FALSE)) ggchord_stop(
      "coord_collinear(): `external = FALSE` constrains circular feature labels; ",
      "remove it to use linear label rows, or use coord_circular()"
    )
    if ((params$type %||% "") %in% c("restriction_site", "seq_center_label") ||
        isTRUE(params$is_primer) || isTRUE(params$is_primer_label)) ggchord_stop(
      "coord_collinear(): circular-only primer, restriction-site and centre-label ",
      "layers require coord_circular(); use geom_feature() and ",
      "geom_feature_label() for linear annotations"
    )
  }
  if (!is.null(plot$scales$get_scales("seq_ring"))) ggchord_stop(
    "coord_collinear(): scale_seq_ring_manual() controls circular radii; ",
    "remove it and use geom_seq(seq_gap = ...) for linear track spacing"
  )
  invisible(plot)
}
