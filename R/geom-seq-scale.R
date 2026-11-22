# A free-standing bp reference uses the coordinate's actual shared bp scale.

ggchord_seq_scale_distance <- function(value, axis, layout, name) {
  context <- if (exists("ggchord_render_context", mode = "function"))
    ggchord_render_context() else NULL
  if (is.unit(value)) {
    if (length(value) != 1L || any(!unitType(value) %in%
        c("mm", "cm", "inches", "points", "bigpts", "picas"))) ggchord_stop(
      "geom_seq_scale(): ", name, " must use one absolute physical unit or one layout-unit number")
    # Absolute units have a fixed conversion and must not open a graphics
    # device merely to build a device-free layout.
    factor <- switch(unitType(value), mm = 1, cm = 10, inches = 25.4,
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

#' Draw a reference length annotation
#' @export
#' @param length Positive scale-bar length in the reference coordinate unit.
#' @param sequence Optional source reference accession for the scale bar.
#' @param instance Optional display instance for the scale bar.
#' @param location Placement of the scale bar in the panel.
#' @param label Optional text for the scale bar.
#' @param height Scale-bar height in millimetres.
#' @param text_size Scale-bar text size in ggplot2 text units.
#' @param colour Fixed colour.
#' @param id Layer identifier. [auto_id()] assigns a plot-local ID and independent default visual scales. Use a name for explicit references, or `NULL` or `NA` to share scales.
#' @return A ggplot2 layer for addition to a plot.
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' one_seq <- seqs[1, , drop = FALSE]
#'
#' ggchord(one_seq) +
#'   geom_seq() +
#'   coord_collinear() +
#'   annotation_seq_scalebar(length = 100, sequence = "A")
annotation_seq_scalebar <- function(length, sequence = NULL, instance = NULL,
    location = "bottomleft", label = NULL, height = 1, text_size = 3,
    colour = "black", id = auto_id()) {
  ggchord_local_console_errors()
  ggchord_check_argument_names(sys.call(), sys.function())
  if (missing(length) || !is.numeric(length) || base::length(length) != 1L ||
      !is.finite(length) || length <= 0)
    ggchord_stop("annotation_seq_scalebar(): length must be positive")
  for (key in c("sequence", "instance")) {
    value <- get(key)
    if (!is.null(value)) ggchord_check_layer_id(value)
  }
  height <- ggchord_physical_mm(height, "scalebar height", TRUE)
  text_size <- ggchord_physical_mm(text_size, "scalebar text size", TRUE)
  if (base::length(height) != 1L || height <= 0 || base::length(text_size) != 1L || text_size <= 0)
    ggchord_stop("annotation_seq_scalebar(): height and text_size must be positive scalars")
  layer <- ggchord_seq_scale_layer(length = length, label = label, location = location,
    tick_height = unit(height, "mm"), size = text_size, colour = colour,
    inherit.aes = FALSE, id = id)
  layer$ggchord_params$sequence <- sequence
  layer$ggchord_params$instance <- instance
  layer
}

ggchord_seq_scale_layer <- function(mapping = NULL, data = NULL, length = NULL, label = NULL,
    location = c("bottomleft", "bottomright", "topleft", "topright"),
    offset = grid::unit(c(0, 6), "mm"), tick_height = grid::unit(1.5, "mm"),
    label_gap = grid::unit(1.5, "mm"), linewidth = .4, colour = "#3D4348", size = 3,
    position = "identity", show.legend = FALSE, inherit.aes = TRUE, ..., id = auto_id()) {
  ggchord_check_argument_names(sys.call(), sys.function())
  fixed <- c(linewidth = !missing(linewidth), colour = !missing(colour), size = !missing(size))
  location <- ggchord_match_arg(location)
  mapping <- ggchord_normalize_mapping(mapping)
  if (!is.null(length) && (!is.numeric(length) || length(length) != 1L ||
      !is.finite(length) || length <= 0)) ggchord_stop(
    "annotation_seq_scalebar(): length must be positive")
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
  ggchord_check_geom_params(alias$dots, GeomChordAxis, "geom_seq_scale()", position)
  layer <- ggplot2::layer(geom = GeomChordAxis, stat = "identity", data = ggchord_empty_seq_scale_geometry(),
    mapping = aes(x = x, y = y, xend = xend, yend = yend,
      label = label, group = group, .component = .component),
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = c(styles, alias$dots))
  layer$ggchord_type <- "seq_scale"
  layer$ggchord_params <- list(id = id, type = "seq_scale", physical_layout = TRUE, length = length, label = label,
    location = location, offset = offset, tick_height = tick_height, label_gap = label_gap)
  ggchord_capture_layer_input(layer, data, mapping, c("scale_length", "label", "x", "y"))
}

ggchord_seq_scale_geometry <- function(data, params, layout) {
  caller <- "annotation_seq_scalebar()"
  if (is.null(data) || !nrow(data)) return(ggchord_empty_seq_scale_geometry())
  ref <- layout$sequence_reference
  ids <- layout$seqs
  instances <- layout$display_instances
  selected <- ids
  if (!is.null(params$sequence))
    selected <- if (is.null(instances)) intersect(selected, params$sequence) else
      intersect(selected, instances$key[instances$accver == params$sequence])
  if (!is.null(params$instance))
    selected <- if (is.null(instances)) intersect(selected, params$instance) else
      intersect(selected, instances$key[instances$instance_id == params$instance])
  if (!length(selected)) ggchord_stop(caller, ": unknown sequence or display instance")
  widths <- vapply(selected, function(id) {
    if (isTRUE(layout$collinear)) return(ggchord_collinear_scale(id, ref))
    arc <- layout$seq_arcs[[id]]
    measure <- ggchord_render_context()$units_per_mm
    if (length(measure) != 2L) measure <- rep(ggchord_units_per_mm(), 2L)
    mm <- sum(sqrt((diff(arc$x) / measure[1L])^2 + (diff(arc$y) / measure[2L])^2))
    mm * measure[1L] / ref$lens[[id]]
  }, numeric(1))
  scale <- widths[1L]
  if (!is.finite(scale) || scale <= 0 || any(abs(widths / scale - 1) > 1e-4))
    ggchord_stop(caller, ": choose a sequence or instance with a unique length scale")
  supplied <- params$length %||% data$scale_length
  if (is.null(supplied)) {
    target <- max(unname(ref$lens)) / 4
    power <- 10^floor(log10(max(target, 1)))
    candidates <- c(1, 2, 5, 10) * power
    supplied <- max(1, candidates[candidates <= max(target, 1)])
  }
  if (length(supplied) == 1L) supplied <- rep(supplied, nrow(data))
  if (!is.numeric(supplied) || length(supplied) != nrow(data) ||
      any(!is.finite(supplied)) || any(supplied <= 0))
    ggchord_stop(caller, ": each bar requires one positive source length")
  labels <- if (is.function(params$label)) params$label(supplied) else params$label %||% data$label
  unit <- layout$sequence_unit %||% "bp"
  automatic_labels <- paste(format(supplied, big.mark = ",", scientific = FALSE, trim = TRUE), unit)
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
      scale_length = supplied[i], sequence_unit = unit,
      bp_length = if (identical(unit, "bp")) supplied[i] else NA_real_,
      measurement_status = status, component_type = "scale_annotation",
      .component = c("line", "line", "major_tick", "major_tick", "text"),
      stringsAsFactors = FALSE)
  })
  ggchord_rbind_fill(output)
}
