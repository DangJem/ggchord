# Feature stacking is solved before polygon generation. The Position object is
# identity after that point; its fields tell the shared layout how to allocate
# the minimum number of radial lanes.

PositionFeatureStack <- ggplot2::ggproto(
  "PositionFeatureStack", ggplot2::PositionIdentity,
  ggchord_feature_stack = TRUE,
  spacing = 0.08,
  side = "strand"
)

#' Stack overlapping genes or features on radial tracks
#'
#' Assigns overlapping intervals on the same sequence and side to the minimum
#' number of radial lanes. Non-overlapping intervals reuse the nearest lane.
#' The resulting anchors are shared with gene label layers because stacking is
#' solved before polygons and labels are generated.
#'
#' @param spacing Positive radial distance between adjacent lanes.
#' @param side Track side. `"strand"` preserves the normal gene convention
#'   (`+` inside and `-` outside); `"outside"` or `"inside"` places both
#'   strands on that side while retaining arrow direction.
#' @return A ggplot2 Position object for `geom_gene()` or `geom_feature()`.
#' @export
position_feature_stack <- function(
    spacing = 0.08,
    side = c("strand", "outside", "inside")) {
  side <- match.arg(side)
  if (!is.numeric(spacing) || length(spacing) != 1L ||
      !is.finite(spacing) || spacing <= 0) {
    ggchord_stop("position_feature_stack(): spacing must be one positive number")
  }
  ggplot2::ggproto(
    NULL, PositionFeatureStack,
    spacing = as.numeric(spacing), side = side
  )
}

ggchord_stack_feature_tracks <- function(data, position) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0L) return(data)
  if (!all(c("accver", "start", "end", "strand") %in% names(data))) {
    ggchord_stop(
      "position_feature_stack() requires accver, start, end, and strand"
    )
  }
  side <- position$side %||% "strand"
  effective_side <- if (side == "strand") {
    ifelse(as.character(data$strand) == "+", "inside", "outside")
  } else {
    rep(side, nrow(data))
  }
  lane <- integer(nrow(data))
  keys <- paste(as.character(data$accver), effective_side, sep = "\r")
  groups <- split(seq_len(nrow(data)), factor(keys, levels = unique(keys)))
  for (idx in groups) {
    lo <- pmin(data$start[idx], data$end[idx])
    hi <- pmax(data$start[idx], data$end[idx])
    ord <- order(lo, hi, idx)
    lane_end <- numeric(0)
    for (local in ord) {
      available <- which(lane_end < lo[local])
      chosen <- if (length(available)) available[1L] else length(lane_end) + 1L
      if (chosen > length(lane_end)) lane_end <- c(lane_end, -Inf)
      lane_end[chosen] <- hi[local]
      lane[idx[local]] <- chosen - 1L
    }
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  data$.feature_stack_lane <- lane
  data$.feature_stack_side <- effective_side
  data$.feature_stack_spacing <- position$spacing
  data
}
