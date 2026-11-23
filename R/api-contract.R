#' Plotting, data and output contracts
#'
#' [ggchord()] and [geom_seq()] offer two ways to declare the same
#' authoritative reference. Build geometry against that reference even when
#' a plot draws no backbone or uses several backbone layers. Pass feature,
#' link and track observations to their respective layers.
#'
#' @section Plot components:
#' Use [layout_sequences()] for order, display orientation, physical distances
#' and source windows. Use a coordinate system for projection and viewport,
#' a Position for local placement and a label layout for text arrangement.
#' For Dotplot, [layout_dotplot()] defines the query and subject axes.
#' The native ggplot2 facet functions divide observations among panels.
#'
#' Use the ordinary ggplot2 aesthetics and visual scales. [scale_for()]
#' assigns a separate visual scale to explicitly named layers.
#' [layout_track()] declares a physical container, [`scale_y_continuous()`][ggplot2::scale_y_continuous]
#' controls its numerical height, and [guide_track_axis()] explains the
#' resulting scale. [theme_ggchord()] controls appearance.
#'
#' @section Coordinates and identity:
#' Source intervals are one-based and closed. Keep biological strand
#' separate from display orientation. Interbase cut coordinates count
#' from the boundary before the first residue. Use bp or aa according to
#' the original sequence. Position offsets and track thickness are in mm.
#' Give each referenced layer an explicit identifier. Source entity,
#' display instance and Fusion panel path have separate identities.
#'
#' @section Reproducible output:
#' [export_ggchord_layout()] reports source, mapped, statistical and
#' geometric stages with their availability and measurement conditions.
#' [map_ggchord_position()] queries source, local and display coordinates.
#' Preview a composition at saved dimensions with [view_ggchord()] and
#' save it with [ggplot2::ggsave()].
#'
#' @seealso [ggchord_data_contract], [ggchord_layout_contract]
#' @name ggchord_api_contract
NULL
