#' Selected ggplot2 helpers for ggchord plots
#'
#' ggchord re-exports a deliberately small set of ggplot2 helpers that are
#' directly useful when constructing, annotating, styling, or saving a chord
#' plot. Cartesian geoms, coordinates, facets, generic scales, and complete
#' theme presets remain in ggplot2 so IDE completion stays focused.
#'
#' @name ggplot2-helpers
#' @keywords internal
NULL

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 aes
#' @export
ggplot2::aes

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 after_stat
#' @export
ggplot2::after_stat

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 after_scale
#' @export
ggplot2::after_scale

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 annotate
#' @export
ggplot2::annotate

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 labs
#' @export
ggplot2::labs

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 ggtitle
#' @export
ggplot2::ggtitle

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 guides
#' @export
ggplot2::guides

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 theme
#' @export
ggplot2::theme

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 element_blank
#' @export
ggplot2::element_blank

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 element_line
#' @export
ggplot2::element_line

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 element_rect
#' @export
ggplot2::element_rect

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 element_text
#' @export
ggplot2::element_text

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 margin
#' @export
ggplot2::margin

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 rel
#' @export
ggplot2::rel

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 ggsave
#' @export
ggplot2::ggsave

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 last_plot
#' @export
ggplot2::last_plot

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 waiver
#' @export
ggplot2::waiver

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 expansion
#' @export
ggplot2::expansion

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 position_identity
#' @export
ggplot2::position_identity

#' @rdname ggplot2-helpers
#' @importFrom ggplot2 PositionIdentity
#' @export
ggplot2::PositionIdentity

#' Physical units and direction arrows
#'
#' These are the original grid functions, available after loading ggchord.
#' Use physical units for theme spacing and arrows for point links.
#' @param x Numeric unit values.
#' @param units Unit names, such as mm, inches or npc.
#' @param data Optional supplementary data for special grid units.
#' @param angle Arrow-head angle in degrees.
#' @param length Arrow-head length as a unit object.
#' @param ends Which ends receive arrows: first, last or both.
#' @param type Arrow-head style: open or closed.
#' @return A grid unit or arrow object.
#' @name grid-helpers
#' @examples
#' unit(2, "mm")
#' arrow(length = unit(2, "mm"), type = "closed")
NULL

#' @rdname grid-helpers
#' @importFrom grid unit
#' @export
grid::unit

#' @rdname grid-helpers
#' @importFrom grid arrow
#' @export
grid::arrow
