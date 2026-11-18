# Feature stacking is solved before polygon generation. The Position object is
# identity after that point; its fields tell the shared layout how to allocate
# the minimum number of radial lanes.

PositionFeatureStack <- ggplot2::ggproto(
  "PositionFeatureStack", ggplot2::PositionIdentity,
  ggchord_feature_stack = TRUE,
  spacing = 0.10,
  side = "strand",
  base_position = NULL,
  direction = "auto"
)

#' Stack overlapping genes or features on local lanes
#'
#' Assigns overlapping intervals on the same sequence and side to radial lanes.
#' Explicit display priority is ordered first, then intervals by decreasing
#' length. Every interval takes the nearest collision-free lane; feature type,
#' name and colour never choose a lane.
#' The resulting anchors are shared with gene label layers because stacking is
#' solved before polygons and labels are generated.
#' In [coord_collinear()], these local lanes become parallel horizontal bands;
#' they do not represent circular physical tracks or radii.
#'
#' The default object can also be selected on a feature layer with
#' `position = "feature_stack"`. Use the constructor to configure the lanes.
#' @param spacing Positive minimum radial distance between adjacent feature
#'   lanes. Circular plasmid layout enlarges the first transition into a main
#'   label corridor and reserves at least one label track between deeper
#'   overlapping feature lanes.
#' @param side Legacy track-side selector retained for compatibility when
#'   used alone. Historical `outside` starts at -0.1 and grows inward;
#'   `inside` starts at +0.1 and grows outward. The default `strand` starts
#'   at +0.1 on plus and -0.1 on minus, growing away from the backbone.
#' @param base_position Optional identity, strand, or plasmid Position applied
#'   before interval-overlap detection and lane allocation. Use identity to
#'   start both strands on the backbone. The default starting bands follow strand.
#' @param direction Growth of successive lanes within each starting band.
#'   `"outward"` follows the positive local normal (up on a horizontal
#'   backbone), `"inward"` the negative normal. `"auto"` follows the sign
#'   of the starting offset, growing outward at zero. Lane zero stays at
#'   its starting offset. Supply legacy side on its own.
#' @examples
#' position_feature_stack(base_position = "identity", direction = "inward")
#' position_feature_stack(base_position = position_plasmid(0.1), direction = "outward")
#' @return A ggplot2 Position object for `geom_gene()` or `geom_feature()`.
#' @export
position_feature_stack <- function(
    spacing = 0.10,
    side = c("strand", "outside", "inside"),
    base_position = NULL, direction = c("auto", "outward", "inward")) {
  if (!missing(side) && !missing(direction))
    ggchord_stop("position_feature_stack(): use direction with base_position, or supply legacy side alone")
  direction <- match.arg(direction)
  side_missing <- missing(side)
  side <- match.arg(side)
  if (!is.numeric(spacing) || length(spacing) != 1L ||
      !is.finite(spacing) || spacing <= 0) {
    ggchord_stop("position_feature_stack(): spacing must be one positive number")
  }
  if (!is.null(base_position)) {
    if (!side_missing) {
      ggchord_stop("position_feature_stack(): `side` cannot be combined with `base_position`")
    }
    base_position <- ggchord_as_feature_position(
      base_position, "position_feature_stack()"
    )
    if (isTRUE(base_position$ggchord_feature_stack)) {
      ggchord_stop("position_feature_stack(): base_position cannot itself be a stack")
    }
  }
  ggplot2::ggproto(
    NULL, PositionFeatureStack,
    spacing = as.numeric(spacing), side = side, base_position = base_position, direction = direction
  )
}
