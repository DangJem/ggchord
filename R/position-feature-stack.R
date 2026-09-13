# Feature stacking is solved before polygon generation. The Position object is
# identity after that point; its fields tell the shared layout how to allocate
# the minimum number of radial lanes.

PositionFeatureStack <- ggplot2::ggproto(
  "PositionFeatureStack", ggplot2::PositionIdentity,
  ggchord_feature_stack = TRUE,
  spacing = 0.10,
  side = "strand",
  base_position = NULL
)

#' Stack overlapping genes or features on radial tracks
#'
#' Assigns overlapping intervals on the same sequence and side to radial lanes.
#' Explicit display priority is ordered first, then intervals by decreasing
#' length. Every interval takes the nearest collision-free lane; feature type,
#' name and colour never choose a lane.
#' The resulting anchors are shared with gene label layers because stacking is
#' solved before polygons and labels are generated.
#'
#' @param spacing Positive minimum radial distance between adjacent feature
#'   lanes. Circular plasmid layout enlarges the first transition into a main
#'   label corridor and reserves at least one label track between deeper
#'   overlapping feature lanes.
#' @param side Legacy track-side selector retained for compatibility when
#'   `base_position` is `NULL`. Its historical names retain their existing
#'   geometry.
#' @param base_position Optional identity, strand, or plasmid Position applied
#'   before interval-overlap detection and lane allocation.
#' @return A ggplot2 Position object for `geom_gene()` or `geom_feature()`.
#' @export
position_feature_stack <- function(
    spacing = 0.10,
    side = c("strand", "outside", "inside"),
    base_position = NULL) {
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
    spacing = as.numeric(spacing), side = side, base_position = base_position
  )
}
