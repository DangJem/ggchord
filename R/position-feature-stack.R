# Feature stacking is solved before polygon generation. The Position object is
# identity after that point; its fields tell the shared layout how to allocate
# the minimum number of radial lanes.

PositionFeatureStack <- ggplot2::ggproto(
  "PositionFeatureStack", ggplot2::PositionIdentity,
  ggchord_feature_stack = TRUE,
  spacing = 0.10,
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
#' @param base_position Position applied before overlap detection and lane
#'   allocation. Use `"identity"` to start both strands on the backbone,
#'   [position_offset()] for one shared offset, or [position_strand()] for
#'   strand-specific offsets. The default starts at +0.1 on the plus strand
#'   and -0.1 on the minus strand.
#' @param direction Growth of successive lanes within each starting band.
#'   `"outward"` follows the positive local normal (up on a horizontal
#'   backbone), `"inward"` the negative normal. `"auto"` follows the sign
#'   of the starting offset, growing outward at zero. Lane zero stays at
#'   its starting offset.
#' @examples
#' position_feature_stack(base_position = "identity", direction = "inward")
#' position_feature_stack(base_position = position_offset(0.1), direction = "outward")
#' @return A ggplot2 Position object for `geom_feature(preset = "gene")` or `geom_feature()`.
#' @export
position_feature_stack <- function(
    spacing = 0.10,
    base_position = position_strand(), direction = c("auto", "outward", "inward")) {
  direction <- match.arg(direction)
  if (!is.numeric(spacing) || length(spacing) != 1L ||
      !is.finite(spacing) || spacing <= 0) {
    ggchord_stop("position_feature_stack(): spacing must be one positive number")
  }
  base_position <- ggchord_as_feature_position(
    base_position %||% position_strand(), "position_feature_stack()"
  )
  if (isTRUE(base_position$ggchord_feature_stack)) {
    ggchord_stop("position_feature_stack(): base_position cannot itself be a stack")
  }
  ggplot2::ggproto(
    NULL, PositionFeatureStack,
    spacing = as.numeric(spacing), base_position = base_position, direction = direction
  )
}
