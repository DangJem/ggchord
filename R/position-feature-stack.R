#' Stack overlapping features in physical lanes
#' @export
#' @param spacing Nonnegative physical gap between adjacent lanes, in millimetres.
#' @param base_position Placement applied before lane assignment.
#' @param direction Placement direction, including strand-dependent lanes.
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
#'       direction = "negative")) +
#'   geom_feature_label_repel(aes(label = label), anchor_layer = "features") +
#'   coord_circular()
position_feature_stack <- function(
    spacing = 1, base_position = position_identity(), direction = "positive") {
  ggchord_check_argument_names(sys.call(), sys.function())
  direction <- match.arg(direction, c("positive", "negative", "strand"))
  spacing <- ggchord_physical_mm(spacing, "position_feature_stack()", nonnegative = TRUE)
  base_position <- ggchord_as_feature_position(base_position, "position_feature_stack()")
  if (isTRUE(base_position$ggchord_feature_stack))
    ggchord_stop("position_feature_stack(): base_position cannot itself be a stack")
  ggproto(NULL, ggplot2::PositionIdentity, ggchord_feature_stack = TRUE,
    spacing = spacing, base_position = base_position, direction = direction)
}
