#' Stack overlapping features in physical lanes
#' @export
#' @param spacing Minimum physical gap between adjacent lanes, in millimetres.
#'   An inward circular stack can widen this gap to fit repelled labels.
#' @param base_position Placement applied before lane assignment. When omitted
#'   for an inward circular stack, the first lane clears the backbone and its
#'   visible inward ticks and labels.
#' @param direction Placement direction, including strand-dependent lanes.
#' @param order Place features by source position (`"source"`) or prefer
#'   higher `display_priority` and then longer occupied intervals (`"priority"`).
#' @return A placement declaration for compatible ggchord layers.
#' @examples
#' library(ggplot2)
#'
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' ggchord(data = one_seq) + geom_seq() +
#'   geom_feature(data = subset(genes, accver == "A"), shape = "arrow",
#'     id = "features",
#'     position = position_feature_stack(
#'       base_position = position_offset(-2), spacing = 1,
#'       direction = "negative", order = "priority")) +
#'   geom_feature_label_repel(aes(label = label), anchor_layer = "features") +
#'   coord_circular()
position_feature_stack <- function(
    spacing = 1, base_position = position_identity(), direction = "positive",
    order = "source") {
  ggchord_check_argument_names(sys.call(), sys.function())
  automatic_base <- missing(base_position)
  direction <- match.arg(direction, c("positive", "negative", "strand"))
  order <- match.arg(order, c("source", "priority"))
  spacing <- ggchord_physical_mm(spacing, "position_feature_stack()", nonnegative = TRUE)
  base_position <- ggchord_as_feature_position(base_position, "position_feature_stack()")
  if (isTRUE(base_position$ggchord_feature_stack))
    ggchord_stop("position_feature_stack(): base_position cannot itself be a stack")
  ggproto(NULL, ggplot2::PositionIdentity, ggchord_feature_stack = TRUE,
    spacing = spacing, base_position = base_position, direction = direction,
    lane_order = order, automatic_base = automatic_base)
}
