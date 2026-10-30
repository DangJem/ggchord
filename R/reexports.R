#' Selected ggplot2 helpers for ggchord plots
#'
#' ggchord re-exports a deliberately small set of ggplot2 helpers that are
#' directly useful when constructing, annotating, styling, or saving a chord
#' plot. Cartesian geoms, coordinates, facets, generic scales, and complete
#' theme presets remain in ggplot2 so IDE completion stays focused.
#'
#' Mapping helpers are `aes()`, `after_stat()` and `after_scale()`;
#' `ggplot2::stage()` combines stages of one aesthetic. Use `labs()` and
#' `ggtitle()` for titles, `guides()` for guide selection, and `theme()` with
#' the `element_*()` helpers for appearance. `annotate()` adds fixed annotations.
#' `margin()` and `rel()` specify theme sizes; `waiver()` and `expansion()`
#' specify automatic or expanded scale settings. `position_identity()` and
#' `PositionIdentity` are the upstream identity placement helpers.
#' `ggsave()` saves a supplied plot; `last_plot()` retrieves ggplot2's last plot.
#' Each symbol is the original ggplot2 function or object, and can also be
#' called with its `ggplot2::` namespace. `unit()` and `arrow()` come from grid.
#'
#' @name ggplot2-helpers
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
