# A free-standing bp reference uses the coordinate's actual shared bp scale.

ggchord_seq_scale_distance <- function(value, axis, layout, name) {
  context <- if (exists("ggchord_render_context", mode = "function"))
    ggchord_render_context() else NULL
  if (grid::is.unit(value)) {
    if (length(value) != 1L || any(!grid::unitType(value) %in%
        c("mm", "cm", "inches", "points", "bigpts", "picas"))) ggchord_stop(
      "geom_seq_scale(): ", name, " must use one absolute physical unit or one layout-unit number")
    # Absolute units have a fixed conversion and must not open a graphics
    # device merely to build a device-free layout.
    factor <- switch(grid::unitType(value), mm = 1, cm = 10, inches = 25.4,
      points = 25.4 / 72.27, bigpts = 25.4 / 72, picas = 12 * 25.4 / 72.27)
    mm <- as.numeric(value) * factor
    units <- context$units_per_mm[[axis]]
    if (is.null(units) || length(units) != 1L || !is.finite(units) || units <= 0)
      units <- (layout$text_units_per_inch %||% .35) / 25.4
    value <- mm * units
  }
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0)
    ggchord_stop("geom_seq_scale(): ", name, " must be finite and non-negative")
  as.numeric(value)
}

ggchord_empty_seq_scale_geometry <- function() data.frame(
  x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
  label = character(), group = integer(), source_row = integer(),
  bp_length = numeric(), measurement_status = character(),
  .component = character(), stringsAsFactors = FALSE)

#' Add a shared base-pair scale bar to a linear display
#'
#' Adds a free-standing reference bar in the empty space outside the sequence
#' envelope. Its length uses the same bp-to-x conversion as every backbone in
#' [coord_collinear()], including independent windows on shared display rows.
#' Display gaps are never interpreted as genomic length. Local reversals and
#' original-source axis labels do not change the bar's bp scale. Chord and
#' circular coordinates have local sequence axes instead and are rejected.
#'
#' The bar does not claim to span one source interval. An explicit length can
#' therefore exceed a single window and still be a valid distance reference
#' in the shared linear metric. With no length, a small whole-number 1/2/5
#' reference is chosen from displayed window lengths. No source coordinates
#' are modified and no missing genomic region is manufactured.
#'
#' @param mapping,data Standard ggplot2 layer inputs. Optional `scale_length`
#'   maps one positive integer bp length per bar; `label` supplies its text.
#'   Multiple bars require explicit `x` and `y` layout positions in `data`
#'   (or their role mappings), to avoid accidental overprinting.
#' @param length Optional fixed positive integer length in bp. Cannot be
#'   combined with a `scale_length` mapping. `NULL` uses input lengths or an
#'   automatic reference length.
#' @param label Optional fixed character label, or a function receiving the
#'   chosen bp lengths. `NULL` uses input labels or formatted `"... bp"` text.
#' @param location Corner of the sequence envelope used for automatic placement.
#' @param offset Two non-negative distances, horizontal then vertical. Absolute
#'   grid units are supported; plain numeric values are layout-unit distances.
#' @param tick_height,label_gap Non-negative absolute grid units or layout-unit
#'   numbers for the end ticks and space between the bar and caption. Physical
#'   conversion uses the current render context when available; device-free
#'   builds retain an explicitly estimated measurement status.
#' @param linewidth,colour,size Fixed line colour/width and text size. Text size
#'   follows ordinary ggplot2 millimetres. These may also use normal mappings
#'   and scales when not supplied as fixed arguments.
#' @param position,show.legend,inherit.aes Standard ggplot2 layer arguments.
#' @param ... Additional fixed aesthetics passed to [ggplot2::layer()].
#' @return One independent ggplot2 layer. Exported axis geometry retains the
#'   bp length, input row and measurement status for each bar.
#' @md
#' @export
#' @examples
#' sequences <- data.frame(accver = c("A", "B"), length = c(1000, 800))
#' ggchord(sequences) + geom_seq(seq_orientation = c(A = 1, B = -1)) +
#'   scale_seq_position_continuous(sequences = "B") +
#'   geom_seq_scale(length = 200) + coord_collinear()
geom_seq_scale <- function(mapping = NULL, data = NULL, length = NULL, label = NULL,
    location = c("bottomleft", "bottomright", "topleft", "topright"),
    offset = grid::unit(c(2, 6), "mm"), tick_height = grid::unit(1.5, "mm"),
    label_gap = grid::unit(1.5, "mm"), linewidth = .4, colour = "#3D4348", size = 3,
    position = "identity", show.legend = FALSE, inherit.aes = FALSE, ...) {
  fixed <- c(linewidth = !missing(linewidth), colour = !missing(colour), size = !missing(size))
  location <- match.arg(location)
  mapping <- ggchord_normalize_mapping(mapping)
  if (!is.null(length) && (!is.numeric(length) || length(length) != 1L ||
      !is.finite(length) || length < 1 || length != floor(length))) ggchord_stop(
    "geom_seq_scale(): length must be one positive integer number of base pairs")
  if (!is.null(length) && "scale_length" %in% names(mapping)) ggchord_stop(
    "geom_seq_scale(): use a fixed length or a scale_length mapping, not both")
  if (!is.null(label) && !is.function(label) &&
      (!is.character(label) || length(label) != 1L || is.na(label))) ggchord_stop(
    "geom_seq_scale(): label must be NULL, one character string, or a labelling function")
  if (length(offset) != 2L) ggchord_stop("geom_seq_scale(): offset must have horizontal and vertical distances")
  dots <- list(...)
  fixed[["colour"]] <- fixed[["colour"]] || "color" %in% names(dots)
  alias <- ggchord_colour_alias(colour, dots, "geom_seq_scale()", sys.call())
  styles <- list(linewidth = linewidth, colour = alias$colour, size = size)
  for (name in names(styles)) if (!fixed[[name]] && name %in% names(mapping)) styles[[name]] <- NULL
  if (is.null(data)) data <- data.frame(.scale_row = 1L)
  layer <- ggplot2::layer(geom = GeomChordAxis, stat = "identity", data = ggchord_empty_seq_scale_geometry(),
    mapping = ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
      label = label, group = group, .component = .component),
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = c(styles, alias$dots))
  layer$ggchord_type <- "seq_scale"
  layer$ggchord_params <- list(type = "seq_scale", physical_layout = TRUE, length = length, label = label,
    location = location, offset = offset, tick_height = tick_height, label_gap = label_gap)
  ggchord_capture_layer_input(layer, data, mapping, c("scale_length", "label", "x", "y"))
}

ggchord_seq_scale_geometry <- function(data, params, layout) {
  caller <- "geom_seq_scale()"
  if (!isTRUE(layout$collinear)) ggchord_stop(caller,
    ": a shared straight bp scale requires coord_collinear(); use local sequence axes for other coordinates")
  if (is.null(data) || !nrow(data)) return(ggchord_empty_seq_scale_geometry())
  ref <- layout$sequence_reference
  scale <- ref$x_scale
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0)
    ggchord_stop(caller, ": the layout does not expose one shared positive bp scale")
  supplied <- params$length %||% data$scale_length
  if (is.null(supplied)) {
    target <- max(unname(ref$lens)) / 4
    power <- 10^floor(log10(max(target, 1)))
    candidates <- c(1, 2, 5, 10) * power
    supplied <- max(1, candidates[candidates <= max(target, 1)])
  }
  if (length(supplied) == 1L) supplied <- rep(supplied, nrow(data))
  if (!is.numeric(supplied) || length(supplied) != nrow(data) ||
      any(!is.finite(supplied)) || any(supplied < 1 | supplied != floor(supplied)))
    ggchord_stop(caller, ": each bar requires one positive integer scale_length")
  labels <- if (is.function(params$label)) params$label(supplied) else params$label %||% data$label
  automatic_labels <- paste0(format(supplied, big.mark = ",", scientific = FALSE, trim = TRUE), " bp")
  if (is.null(labels)) labels <- automatic_labels
  if (length(labels) == 1L) labels <- rep(labels, nrow(data))
  if (!is.character(labels) || length(labels) != nrow(data)) ggchord_stop(caller, ": labels must return one character string per bar")
  missing <- is.na(labels)
  labels[missing] <- automatic_labels[missing]
  has_position <- c("x", "y") %in% names(data)
  if (any(has_position) && !all(has_position)) ggchord_stop(caller, ": explicit bar placement requires both x and y")
  if (nrow(data) > 1L && !all(has_position)) ggchord_stop(caller, ": multiple bars require explicit x/y positions")
  if (all(has_position) && (!is.numeric(data$x) || !is.numeric(data$y) ||
      any(!is.finite(data$x)) || any(!is.finite(data$y)))) ggchord_stop(caller, ": x/y positions must be finite layout coordinates")
  offset <- vapply(seq_len(2L), function(i) ggchord_seq_scale_distance(params$offset[i],
    c("x", "y")[i], layout, "offset"), numeric(1))
  height <- ggchord_seq_scale_distance(params$tick_height, "y", layout, "tick_height")
  gap <- ggchord_seq_scale_distance(params$label_gap, "y", layout, "label_gap")
  arcs <- ggchord_rbind_fill(layout$seq_arcs)
  bounds <- c(xmin = min(arcs$x), xmax = max(arcs$x), ymin = min(arcs$y), ymax = max(arcs$y))
  if (is.data.frame(layout$track_registry) && nrow(layout$track_registry)) {
    bands <- layout$track_registry
    y <- unname(ref$y[bands$accver])
    bounds[c("ymin", "ymax")] <- range(c(bounds[c("ymin", "ymax")],
      y + bands$normal_min, y + bands$normal_max))
  }
  bottom <- grepl("bottom", params$location)
  left <- grepl("left", params$location)
  context <- if (exists("ggchord_render_context", mode = "function")) ggchord_render_context() else NULL
  measured <- !is.null(context) && identical(context$status, "measured") &&
    length(context$units_per_mm) == 2L && all(is.finite(context$units_per_mm)) &&
    all(context$units_per_mm > 0)
  status <- if (measured) "measured" else "estimated"
  output <- lapply(seq_len(nrow(data)), function(i) {
    width <- supplied[i] * scale
    x <- if (all(has_position)) data$x[i] else if (left) bounds[["xmin"]] + offset[1L] else bounds[["xmax"]] - offset[1L] - width
    y <- if (all(has_position)) data$y[i] else if (bottom) bounds[["ymin"]] - offset[2L] else bounds[["ymax"]] + offset[2L]
    text_y <- y + (if (bottom) -1 else 1) * (height / 2 + gap)
    group <- (i - 1L) * 4L
    data.frame(x = c(x, x + width, x, x + width, x + width / 2),
      y = c(y, y, y - height / 2, y - height / 2, text_y),
      xend = c(NA, NA, x, x + width, NA),
      yend = c(NA, NA, y + height / 2, y + height / 2, NA),
      label = c(NA_character_, NA_character_, NA_character_, NA_character_, labels[i]),
      group = group + c(1L, 1L, 2L, 3L, 4L), source_row = i,
      bp_length = supplied[i], measurement_status = status,
      .component = c("line", "line", "major_tick", "major_tick", "text"),
      stringsAsFactors = FALSE)
  })
  ggchord_rbind_fill(output)
}
