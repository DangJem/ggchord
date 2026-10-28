# Decorative endpoint marks for cropped linear sequence views.

#' Mark cropped ends of linear sequences
#'
#' Adds two short diagonal strokes at each truncated sequence end in
#' \code{coord_collinear()}. Metadata returned by \code{focus_ggchord_data()}
#' determines the truncated genomic ends; reversed sequence orientation moves
#' each mark to the corresponding displayed end. These marks are decorations
#' and do not change the source sequence or feature coordinates.
#'
#' @param mapping,data Standard layer mapping and sequence data. By default,
#'   use the plot's sequence table. Truncation is read from
#'   \code{.truncated_left}/\code{.truncated_right}, with source interval and
#'   length metadata as fallback when available.
#' @param break_width,break_height Width and height of the complete break mark
#'   in linear layout units. These are layout distances, not millimetres.
#' @param linewidth,colour Fixed line appearance. Ordinary visual aesthetics
#'   can also be mapped and controlled by their ggplot2 scales.
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional parameters passed to \code{ggplot2::geom_path()}.
#' @return A ggplot2 layer. Other ggchord coordinate systems are unsupported.
#' @export
#' @examples
#' seqs <- data.frame(accver = "A", length = 1000)
#' focused <- focus_ggchord_data(seqs,
#'   loci = data.frame(accver = "A", start = 200, end = 700))
#' ggchord(focused$seq_data) + geom_seq() + geom_seq_break() +
#'   coord_collinear()
geom_seq_break <- function(mapping = NULL, data = NULL,
                           break_width = .02, break_height = .04,
                           linewidth = .4, colour = "#3D4348",
                           position = "identity", show.legend = FALSE,
                           inherit.aes = FALSE, ...) {
  linewidth_supplied <- !missing(linewidth)
  colour_supplied <- !missing(colour)
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)
  dots <- list(...)
  colour_supplied <- colour_supplied || "color" %in% names(dots)
  alias <- ggchord_colour_alias(colour, dots, "geom_seq_break()", sys.call())
  mapping <- ggchord_normalize_mapping(mapping)
  for (name in c("break_width", "break_height")) {
    value <- if (name == "break_width") break_width else break_height
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value <= 0) {
      ggchord_stop("geom_seq_break(): ", name, " must be one finite positive number")
    }
  }
  styles <- list(linewidth = linewidth, colour = alias$colour)
  if (!linewidth_supplied && "linewidth" %in% names(mapping)) styles$linewidth <- NULL
  if (!colour_supplied && "colour" %in% names(mapping)) styles$colour <- NULL
  empty <- ggchord_empty_seq_break_geometry()
  layer <- do.call(ggplot2::geom_path, c(list(
    data = empty, mapping = ggplot2::aes(x = x, y = y, group = group),
    position = position, show.legend = show.legend, inherit.aes = inherit.aes
  ), styles, alias$dots))
  layer$ggchord_type <- "seq_break"
  layer$ggchord_params <- list(type = "seq_break",
    break_width = as.numeric(break_width), break_height = as.numeric(break_height))
  ggchord_capture_layer_input(layer, data, mapping,
    c("accver", "length", ".truncated_left", ".truncated_right"))
}

ggchord_empty_seq_break_geometry <- function() data.frame(
  x = numeric(), y = numeric(), group = integer(), source_row = integer(),
  accver = character(), source_end = character(), display_end = character(),
  .component = character(), stringsAsFactors = FALSE)

ggchord_seq_break_flags <- function(data) {
  n <- nrow(data)
  number <- function(name) {
    value <- data[[name]]
    if (is.null(value) || !is.numeric(value)) rep(NA_real_, n) else value
  }
  first <- number(".source_start")
  origin <- number(".source_origin")
  first[!is.finite(first)] <- origin[!is.finite(first)]
  focused_start <- number(".focus_start")
  first[!is.finite(first)] <- focused_start[!is.finite(first)]
  last <- number(".source_end")
  focused_end <- number(".focus_end")
  last[!is.finite(last)] <- focused_end[!is.finite(last)]
  full_length <- number(".source_length")
  left <- is.finite(first) & first > 1
  right <- is.finite(last) & is.finite(full_length) & last < full_length
  for (side in c("left", "right")) {
    value <- data[[paste0(".truncated_", side)]]
    if (is.null(value)) next
    if (!is.logical(value)) ggchord_stop(
      "geom_seq_break(): .truncated_left/right must contain logical values")
    known <- !is.na(value)
    if (side == "left") left[known] <- value[known] else right[known] <- value[known]
  }
  list(left = left, right = right)
}

ggchord_seq_break_geometry <- function(data, params, layout) {
  if (!isTRUE(layout$collinear)) ggchord_stop(
    "geom_seq_break() requires coord_collinear(); remove the layer for other coordinates")
  if (is.null(data) || !nrow(data)) return(ggchord_empty_seq_break_geometry())
  ggchord_require_columns(data, "accver", "geom_seq_break()")
  ref <- layout$sequence_reference
  ids <- as.character(data$accver)
  if (anyNA(ids) || any(!nzchar(ids))) ggchord_stop(
    "geom_seq_break(): accessions must be non-missing and non-empty")
  flags <- ggchord_seq_break_flags(data)
  width <- params$break_width %||% .02
  height <- params$break_height %||% .04
  out <- list()
  group <- 0L
  for (i in seq_len(nrow(data))) {
    id <- ids[i]
    if (!id %in% names(ref$lens)) next
    orientation <- unname(ref$orientation[id])
    shift <- if (is.null(ref$x_offset)) 0 else unname(ref$x_offset[id])
    for (side in c("left", "right")) {
      if (!flags[[side]][i]) next
      position <- if (side == "left") 0 else unname(ref$lens[id])
      x <- shift + (if (orientation == 1) position else
        unname(ref$lens[id]) - position) * ref$x_scale
      y <- unname(ref$y[id])
      displayed <- if ((side == "left") == (orientation == 1)) "left" else "right"
      for (stroke in c(-1, 1)) {
        group <- group + 1L
        out[[length(out) + 1L]] <- data.frame(
          x = x + stroke * width / 4 + c(-1, 1) * width / 4,
          y = y + c(-1, 1) * height / 2,
          group = group, source_row = i, accver = id,
          source_end = side, display_end = displayed,
          .component = "seq_break", stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out)) do.call(rbind, out) else ggchord_empty_seq_break_geometry()
}
