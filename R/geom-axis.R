# geom-axis.R - axis layer
# Fetches pre-computed axis lines, tick marks, and label data from the package environment
# Axis parameters are specified in this layer and stored for use at print time

#' Add an axis layer
#'
#' Draws axes for each sequence in the chord diagram (including axis lines, major/minor ticks, and labels).
#' Axis parameters (spacing, tick count/length, label size/orientation, etc.) are specified here.
#'
#' @param mapping Default NULL (uses pre-computed data)
#' @param data Default NULL (retrieved automatically from the layout)
#' @param show_axis Logical. Whether to show the axis, default TRUE
#' @param axis_gap Optional numeric/vector. Spacing between sequence and axis, default 0.05
#' @param axis_tick_major_number Optional integer/vector. Number of major ticks, default 3
#' @param axis_tick_major_length Optional numeric/vector. Major tick length ratio, default 0.02
#' @param axis_tick_minor_number Optional integer/vector. Number of minor ticks, default 4
#' @param axis_tick_minor_length Optional numeric/vector. Minor tick length ratio, default 0.01
#' @param axis_label_size Optional numeric/vector. Tick label font size, default 3
#' @param axis_label_offset Optional numeric/vector. Label offset ratio, default 2
#' @param axis_label_orientation Optional character/numeric/vector. Label
#'   orientation, default "parallel". Accepted values: "horizontal" (text stays
#'   horizontal), "parallel" (text runs parallel to the axis, i.e. along the
#'   arc), "perpendicular" (text runs perpendicular to the axis, i.e. along the
#'   radial direction), or a numeric angle in degrees (ggplot2 convention:
#'   counter-clockwise from horizontal, in the final rendered plot space). A
#'   vector or named vector can be used to specify a different orientation per
#'   sequence.
#' @param axis_label_hide_overlaps Logical, default FALSE. When TRUE, axis
#'   labels whose boxes would overlap the plot content (sequence arcs, genes,
#'   ribbons) or other axis labels are automatically hidden.
#' @param line_params,tick_params,text_params Named lists of fixed style
#'   arguments for the axis path, tick segments and tick labels respectively.
#'   These override the corresponding `ggchord.axis.*` theme elements.
#' @param ... Shared fixed style arguments. `colour`, `alpha` and
#'   `na.rm` apply to all three components; line-specific and
#'   text-specific arguments are routed only to compatible geoms.
#'
#' @return A list of ggplot2 layers
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' p <- ggchord(seq_data_example) + geom_seq() + geom_axis()
#' p
geom_axis <- function(mapping = NULL, data = NULL,
                      show_axis = NULL,
                      axis_gap = NULL,
                      axis_tick_major_number = NULL,
                      axis_tick_major_length = NULL,
                      axis_tick_minor_number = NULL,
                      axis_tick_minor_length = NULL,
                      axis_label_size = NULL,
                      axis_label_offset = NULL,
                      axis_label_orientation = NULL,
                      axis_label_hide_overlaps = FALSE,
                      line_params = list(),
                      tick_params = list(),
                      text_params = list(),
                      ...) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  split_params <- ggchord_axis_params(
    list(...), line_params, tick_params, text_params
  )

  empty_id <- data.frame(x = numeric(0), y = numeric(0),
                         seq_id = character(0))
  empty_seg <- data.frame(x0 = numeric(0), y0 = numeric(0),
                          x1 = numeric(0), y1 = numeric(0),
                          label = character(0), label_x = numeric(0),
                          label_y = numeric(0), size = numeric(0),
                          label_hjust = numeric(0), label_vjust = numeric(0),
                          label_angle = numeric(0),
                          seq_id = character(0))

  path_layer <- do.call(ggplot2::geom_path, c(list(
    data = empty_id,
    mapping = ggplot2::aes(x = x, y = y, group = seq_id),
    inherit.aes = FALSE, show.legend = FALSE
  ), split_params$line))
  path_layer$ggchord_type <- "axis_line"
  path_layer$ggchord_theme_element <- "ggchord.axis.line"
  path_layer$ggchord_params <- list(
    type                    = "axis",
    show_axis               = show_axis,
    axis_gap                = axis_gap,
    axis_tick_major_number  = axis_tick_major_number,
    axis_tick_major_length  = axis_tick_major_length,
    axis_tick_minor_number  = axis_tick_minor_number,
    axis_tick_minor_length  = axis_tick_minor_length,
    axis_label_size         = axis_label_size,
    axis_label_offset       = axis_label_offset,
    axis_label_orientation  = axis_label_orientation,
    axis_label_hide_overlaps = axis_label_hide_overlaps
  )
  path_layer <- ggchord_capture_layer_input(
    path_layer, data, mapping, c("seq_id", "length")
  )
  path_layer <- ggchord_add_legacy_scale(
    path_layer,
    (!missing(axis_tick_major_number) && !is.null(axis_tick_major_number)) ||
      (!missing(axis_tick_minor_number) && !is.null(axis_tick_minor_number)),
    "axis_tick_major_number/axis_tick_minor_number", "seq_position",
    "scale_seq_position_continuous(breaks = ..., minor_breaks = ...)"
  )

  seg_layer <- do.call(ggplot2::geom_segment, c(list(
    data = empty_seg,
    mapping = ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
    inherit.aes = FALSE, show.legend = FALSE
  ), split_params$ticks))
  seg_layer$ggchord_type <- "axis_seg"
  seg_layer$ggchord_theme_element <- "ggchord.axis.ticks"
  seg_layer <- ggchord_capture_layer_input(
    seg_layer, data, mapping, c("seq_id", "length")
  )

  text_layer <- do.call(ggplot2::geom_text, c(list(
    data = empty_seg[integer(0), ],
    mapping = ggplot2::aes(x = label_x, y = label_y,
                  label = label, size = I(size),
                  hjust = label_hjust, vjust = label_vjust,
                  angle = label_angle),
    inherit.aes = FALSE, show.legend = FALSE
  ), split_params$text))
  text_layer$ggchord_type <- "axis_text"
  text_layer$ggchord_theme_element <- "ggchord.axis.text"
  text_layer <- ggchord_capture_layer_input(
    text_layer, data, mapping, c("seq_id", "length")
  )

  list(path_layer, seg_layer, text_layer)
}

#' Route geom_axis style arguments to compatible child geoms
#' @noRd
ggchord_axis_params <- function(dots, line, ticks, text) {
  for (item in list(line = line, ticks = ticks, text = text)) {
    if (!is.list(item) || (length(item) > 0 && is.null(names(item)))) {
      ggchord_stop("geom_axis(): *_params arguments must be named lists")
    }
  }
  normalize <- function(x) {
    if ("color" %in% names(x)) {
      if ("colour" %in% names(x)) {
        ggchord_stop("geom_axis(): use only one of colour and color")
      }
      names(x)[names(x) == "color"] <- "colour"
    }
    x
  }
  dots <- normalize(dots)
  line <- normalize(line)
  ticks <- normalize(ticks)
  text <- normalize(text)

  common <- c("colour", "alpha", "na.rm")
  line_names <- c(common, "linewidth", "linetype", "lineend", "linejoin")
  tick_names <- c(line_names, "arrow", "arrow.fill")
  text_names <- c(
    common, "family", "fontface", "lineheight", "parse", "check_overlap"
  )
  known <- union(line_names, union(tick_names, text_names))
  if (length(dots) > 0 && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    ggchord_stop("geom_axis(): all shared style arguments in ... must be named")
  }
  unknown <- setdiff(names(dots), known)
  if (length(unknown) > 0) {
    if ("show_legend" %in% unknown) {
      ggchord_stop("geom_axis(): show_legend was removed because axes never create a legend")
    }
    ggchord_stop(
      "geom_axis(): unsupported shared style argument(s): ",
      paste(unknown, collapse = ", "),
      ". Use line_params, tick_params or text_params for component styles."
    )
  }
  merge <- function(shared, specific) {
    shared <- shared[setdiff(names(shared), names(specific))]
    c(shared, specific)
  }
  allowed <- list(line = line_names, ticks = tick_names, text = text_names)
  supplied <- list(line = line, ticks = ticks, text = text)
  for (nm in names(supplied)) {
    bad <- setdiff(names(supplied[[nm]]), allowed[[nm]])
    if (length(bad) > 0) {
      ggchord_stop(
        "geom_axis(): unsupported ", nm, "_params argument(s): ",
        paste(bad, collapse = ", ")
      )
    }
  }
  list(
    line = merge(dots[intersect(names(dots), line_names)], line),
    ticks = merge(dots[intersect(names(dots), tick_names)], ticks),
    text = merge(dots[intersect(names(dots), text_names)], text)
  )
}
