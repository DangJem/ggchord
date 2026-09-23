# Primer binding-site geometry and labels.

#' Draw primer binding ranges on a circular sequence
#'
#' Primers are slender directional arcs drawn between the two lines of a
#' double sequence backbone by default, not arrows in the feature stack. Add
#' [geom_primer_label_repel()] for the matching unboxed label and leader.
#'
#' @param mapping,data Standard layer inputs.
#' @param side Backbone lane used by the primer arc. `"between"` uses the
#'   backbone centre line; `"outside"` and `"inside"` use `arc_offset`.
#' @param arc_offset Distance from the backbone centre when `side` is
#'   `"outside"` or `"inside"`.
#' @param arc_linewidth Radial width of the primer arc.
#' @param colour Primer arc colour.
#' @param position,show.legend,inherit.aes Standard layer arguments.
#' @param ... Additional parameters passed to [geom_feature()].
#' @return A ggplot2 layer.
#' @export
#' @examples
#' primer <- data.frame(
#'   accver = "circle", start = 100, end = 120,
#'   strand = "+", name = "sequencing primer"
#' )
#' ggchord(data.frame(accver = "circle", length = 1000)) +
#'   geom_seq() + geom_primer(data = primer) +
#'   geom_primer_label_repel(data = primer) + coord_circular()
geom_primer <- function(mapping = NULL, data = NULL,
                        side = c("between", "outside", "inside"),
                        arc_offset = .045,
                        arc_linewidth = .026,
                        colour = "#A020F0",
                        position = NULL,
                        show.legend = FALSE,
                        inherit.aes = FALSE, ...) {
  side <- match.arg(side)
  values <- c(arc_offset, arc_linewidth)
  if (!is.numeric(values) || any(!is.finite(values)) ||
      arc_offset < 0 || arc_linewidth <= 0) {
    ggchord_stop("geom_primer(): arc_offset and arc_linewidth must be finite non-negative/positive values")
  }
  if (is.null(position)) {
    offset <- switch(side,
      between = 0,
      outside = arc_offset,
      inside = -arc_offset
    )
    position <- position_plasmid(offset)
  }
  layer <- geom_feature(
    mapping = mapping, data = data, feature_shape = "primer_arc",
    feature_width = arc_linewidth, arrow_head_length = arc_linewidth,
    arrow_head_width = 1, short_feature = "auto", position = position,
    show.legend = show.legend, inherit.aes = inherit.aes,
    feature_fill = colour, colour = colour, ...
  )
  inherited_transform <- layer$ggchord_input_transform
  layer$ggchord_input_transform <- function(x) {
    if (!"feature_type" %in% names(x)) x$feature_type <- "primer"
    if (!"feature_label" %in% names(x) && "name" %in% names(x)) {
      x$feature_label <- x$name
    }
    x$annotation_class <- "primer"
    inherited_transform(x)
  }
  layer$ggchord_params$feature_role <- "primer"
  layer$ggchord_params$is_primer <- TRUE
  layer$ggchord_params$primer_side <- side
  layer
}

#' Arrange primer labels and leaders
#'
#' Primer labels are unboxed, use the primer colour, and show one-based
#' inclusive binding coordinates by default. Their leaders start at the
#' directional primer arc's arrow tip (the primer 3-prime end).
#'
#' @inheritParams geom_primer
#' @param side Side on which the external primer label is arranged.
#' @param show_location Include the binding range in the label.
#' @param max_overlaps Maximum unresolved overlaps.
#' @return A composite text and leader-line layer.
#' @export
geom_primer_label_repel <- function(
    mapping = NULL, data = NULL,
    side = c("outside", "inside"),
    show_location = TRUE,
    colour = "#A020F0",
    max_overlaps = Inf,
    position = NULL,
    show.legend = FALSE, inherit.aes = FALSE, ...) {
  side <- match.arg(side)
  if (!is.logical(show_location) || length(show_location) != 1L ||
      is.na(show_location)) {
    ggchord_stop("geom_primer_label_repel(): show_location must be TRUE or FALSE")
  }
  if (is.null(position)) {
    # `side` selects the exterior label field, not a second biological anchor.
    # Keep the label's source on the primer arc between the backbone lines.
    position <- position_plasmid(0)
  }
  layer <- geom_feature_label_repel(
    mapping = mapping, data = data, label_layout = "callout",
    label_side = side, max_overlaps = max_overlaps, external = TRUE,
    position = position, show.legend = show.legend,
    inherit.aes = inherit.aes, colour = colour, ...
  )
  layer$ggchord_input_transform <- function(x) {
    ggchord_require_columns(
      x, c("start", "end"), "geom_primer_label_repel()"
    )
    name <- if ("label" %in% names(x)) x$label else if (
      "name" %in% names(x)) x$name else rep("primer", nrow(x))
    range <- paste0("(", as.integer(x$start), " .. ",
      as.integer(x$end), ")")
    x$feature_label <- if (show_location) {
      paste(as.character(name), range)
    } else as.character(name)
    x$annotation_class <- "primer"
    x$primer_name <- as.character(name)
    x$primer_start <- as.integer(x$start)
    x$primer_end <- as.integer(x$end)
    x$primer_show_location <- show_location
    x$.feature_label_anchor <- "head"
    x$feature_label_mode <- "external"
    x$feature_label_fill <- NA_character_
    ggchord_feature_label_data(x)
  }
  layer$mapping[["annotation_class"]] <- ggplot2::aes(
    annotation_class = I(annotation_class)
  )$annotation_class
  layer$mapping[["primer_name"]] <- ggplot2::aes(
    primer_name = I(primer_name)
  )$primer_name
  layer$mapping[["primer_start"]] <- ggplot2::aes(
    primer_start = I(primer_start)
  )$primer_start
  layer$mapping[["primer_end"]] <- ggplot2::aes(
    primer_end = I(primer_end)
  )$primer_end
  layer$mapping[["primer_show_location"]] <- ggplot2::aes(
    primer_show_location = I(primer_show_location)
  )$primer_show_location
  layer$ggchord_params$is_primer_label <- TRUE
  layer$ggchord_params$feature_label_external <- TRUE
  # The composite label geom styles its segment separately from its text.
  # A text colour alone would leave the primer leader in the feature-theme
  # grey, even though both are one purple annotation in the reference.
  layer$geom_params$segment_params$colour <- colour
  layer
}
