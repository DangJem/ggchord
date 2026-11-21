# Decorative endpoint marks for cropped linear sequence views.

#' Mark display breaks on reference paths
#' @export
#' @param mapping Aesthetic mapping from `aes()`. A local mapping overrides the matching inherited channel.
#' @param data Optional reference table or function selecting from the registered reference.
#' @param mark Symbol used at a display discontinuity.
#' @param size Physical size of the symbol in millimetres.
#' @param position A compatible Position for placement of the layer geometry.
#' @param ... Additional named visual styles and parameters accepted by this layer.
#' @param na.rm Suppress warnings for missing observations that the layer can omit.
#' @param show.legend Logical or named logical indicating whether the layer contributes legend keys.
#' @param inherit.aes Whether to inherit plot mappings and evaluate them on the layer data.
#' @param id Optional plot-local layer identifier. Name a layer explicitly to scope a scale, anchor a label or export it.
#' @return A ggplot2 layer for addition to a plot.
#' @details The source is the complete reference table with `accver` and `length`. A backbone can select an authoritative subset. An explicit complete table in geom_seq() can register a reference when the plot has no global data.
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' windows <- data.frame(accver = "A", start = 1, end = 500, window_id = "w1")
#'
#' display_windows <- data.frame(
#'   instance_id = c("A_left", "A_right"), accver = "A",
#'   start = c(1, 600), end = c(300, 1000), row_id = "A"
#' )
#' ggchord(one_seq) +
#'   layout_sequences(windows = display_windows, segment_gap = 3) +
#'   geom_seq() +
#'   geom_seq_break(mark = "double_slash") +
#'   coord_collinear()
geom_seq_break <- function(mapping = NULL, data = NULL, mark = "slash", size = 2,
    position = "identity", ..., na.rm = FALSE, show.legend = NA,
    inherit.aes = TRUE, id = NULL) {
  ggchord_check_argument_names(sys.call(), sys.function())
  mark <- match.arg(mark, c("slash", "double_slash", "zigzag"))
  size <- ggchord_physical_mm(size, "geom_seq_break(size)", nonnegative = TRUE)
  if (length(size) != 1L || size <= 0)
    ggchord_stop("geom_seq_break(): size must be a positive scalar")
  position <- ggchord_as_feature_position(position, "geom_seq_break()")
  if (ggchord_local_position(position) && !identical(position$placement, "offset"))
    ggchord_stop("geom_seq_break(): use identity or position_offset()")
  dots <- ggchord_colour_dots(list(...), "geom_seq_break()")
  ggchord_check_geom_params(dots, ggplot2::GeomPath, "geom_seq_break()")
  layer <- ggplot2::layer(
    data = ggchord_empty_seq_break_geometry(),
    mapping = aes(x = x, y = y, group = group),
    geom = ggplot2::GeomPath, stat = "identity", position = "identity",
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = c(list(na.rm = na.rm), dots)
  )
  layer$ggchord_type <- "seq_break"
  layer$ggchord_params <- list(id = id, type = "seq_break", mark = mark,
    size_mm = size, feature_position = position, physical_layout = TRUE)
  ggchord_capture_layer_input(layer, data, mapping, character())
}

ggchord_empty_seq_break_geometry <- function() data.frame(
  x = numeric(), y = numeric(), group = integer(), source_row = integer(),
  accver = character(), source_end = character(), display_end = character(),
  .component = character(), stringsAsFactors = FALSE)

# Automatic structural discontinuities are part of the sequence path style,
# independent of an optional user geom_seq_break() decoration layer.
ggchord_segment_break_geometry <- function(layout) {
  data <- layout$sequence_layout
  if (!"row_id" %in% names(data)) return(NULL)
  result <- list()
  for (i in seq_len(nrow(data))) for (side in c("left", "right")) {
    if (!isTRUE(data[[paste0("break_", side)]][i])) next
    x <- data[[if (side == "left") "xmin" else "xmax"]][i]
    width <- min(.02, data$gap_bp[i] * data$units_per_bp[i] * .4)
    for (stroke in c(-1, 1)) {
      result[[length(result) + 1L]] <- data.frame(
        x = x + stroke * width / 4 + c(-1, 1) * width / 4,
        y = data$y[i] + c(-1, 1) * .02,
        accver = data$accver[i],
        group = paste(data$accver[i], "segment-break", side, stroke, sep = "\r"),
        .component = "break", .decoration = "segment_break", stringsAsFactors = FALSE)
    }
  }
  if (length(result)) do.call(rbind, result) else NULL
}

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
  if (is.null(data) || !nrow(data)) return(ggchord_empty_seq_break_geometry())
  ggchord_require_columns(data, "accver", "geom_seq_break()")
  ids <- as.character(data$accver)
  ref <- layout$sequence_reference
  flags <- ggchord_seq_break_flags(data)
  row_breaks <- layout$segment_positions
  size <- params$size_mm * ggchord_units_per_mm()
  offsets <- ggchord_resolve_position_offset(params$feature_position$offset %||% 0,
    names(ref$lens), "offset", "geom_seq_break()")
  result <- list()
  add <- function(id, row, side, position) {
    basis <- ggchord_track_frame(position, id, layout)
    normal <- drop(basis$normal); tangent <- c(normal[2L], -normal[1L])
    center <- drop(basis$base) + normal * offsets[[id]][[1L]] * ggchord_units_per_mm()
    pattern <- switch(params$mark,
      slash = list(rbind(c(-.25, -.5), c(.25, .5))),
      double_slash = list(rbind(c(-.50, -.5), c(0, .5)),
        rbind(c(0, -.5), c(.50, .5))),
      zigzag = list(rbind(c(-.5, -.5), c(-.1, .5), c(.1, -.5), c(.5, .5))))
    for (points in pattern) {
      xy <- outer(points[, 1L] * size, tangent) + outer(points[, 2L] * size, normal)
      result[[length(result) + 1L]] <<- data.frame(
        x = center[1L] + xy[, 1L], y = center[2L] + xy[, 2L],
        group = length(result) + 1L, source_row = row, accver = id,
        source_end = side, display_end = side, .component = "display_break")
    }
  }
  for (i in seq_along(ids)) {
    id <- ids[i]; len <- ref$lens[[id]]
    if (is.null(len)) next
    if (isTRUE(flags$left[i])) add(id, i, "left", 0)
    if (isTRUE(flags$right[i])) add(id, i, "right", len)
  }
  ggchord_rbind_fill(result) %||% ggchord_empty_seq_break_geometry()
}
