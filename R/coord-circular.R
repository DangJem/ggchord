# Dedicated coordinate contract for one circular sequence.

CoordCircular <- ggproto(
  "CoordCircular", CoordGgchord,
  resolve_annotation_registry = function(self, layer_geometry, layout) {
    ggchord_resolve_circular_annotation_registry(
      layer_geometry, layout, coord = self
    )
  }
)

#' Project a reference sequence onto a circle
#'
#' @inheritParams coord_chord
#' @return A coordinate system for a ggchord plot.
#' @export
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
#' ggchord(data = one_seq) +
#'   geom_seq() +
#'   geom_gene(data = subset(genes, accver == "A"), aes(fill = type)) +
#'   layout_sequences(angular_gap = 12) +
#'   coord_circular(start_angle = 90)
coord_circular <- function(start_angle = 0, direction = "clockwise",
    xlim = NULL, ylim = NULL, expand = FALSE, clip = "on", fit = "geometry") {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  direction <- ggchord_coord_direction(start_angle, direction)
  settings <- ggchord_coord_settings(1, xlim, ylim, expand, clip, fit)
  ggproto(NULL, CoordCircular, !!!settings,
    rotation = start_angle %% 360, reference_direction = direction,
    ggchord_circular = TRUE)
}
