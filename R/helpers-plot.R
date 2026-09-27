#' Colour-vision-friendly categorical palette
#' @noRd
chord_palette_discrete <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#56B4E9", "#E69F00", "#6F6F6F", "#F0E442"
)

#' Generate a default categorical palette
#'
#' Returns a colour-vision-friendly discrete palette. Larger palettes use a
#' qualitative HCL palette rather than interpolating through muddy midpoints.
#' @param n Number of colors requested
#' @return A character vector of \code{n} colors
#' @noRd
chord_default_palette <- function(n) {
  if (n <= 0) {
    character(0)
  } else if (n <= length(chord_palette_discrete)) {
    chord_palette_discrete[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, palette = "Dark 3")
  }
}

#' Missing value handling operator
#'
#' Used to safely handle NULL values: returns y if x is NULL, otherwise returns x
#'
#' @param x Any R object (may be NULL)
#' @param y Default value to return when x is NULL
#' @return x if x is not NULL, otherwise y
#' @noRd
if_null_else <- function (x, y)
{
  if (is.null(x)) y else x
}


#' Calculate plot extremes
#'
#' Extracts x/y coordinate extremes from all plot elements (sequence arcs, ribbons, gene arrows, etc.) for adjusting the plot range
#'
#' @param allRibbon data.frame, ribbon data (with x, y columns), default NULL
#' @param seqArcs List, sequence arc data (each element is a data frame with x, y, accver), default NULL
#' @param axisLines data.frame, axis line data (with x, y, accver columns), default NULL
#' @param axisTicks data.frame, tick mark data (with x0, y0, x1, y1, label_x, label_y columns), default NULL
#' @param gene_arrows data.frame, gene label data (with text_x, text_y columns), default NULL
#' @param gene_polys data.frame, gene arrow polygon data (with x, y columns), default NULL
#' @param show_axis Logical, whether to include extreme value calculation for axis-related elements, default FALSE
#' @return List containing x_min (minimum x), x_max (maximum x), y_min (minimum y), y_max (maximum y)
#' @noRd
get_plot_extremes <- function(allRibbon = NULL, seqArcs = NULL,
                                axisLines = NULL, axisTicks = NULL,
                                gene_arrows = NULL, gene_polys = NULL,
                                seq_labels = NULL, show_axis = FALSE) {
  x_min <- Inf
  x_max <- -Inf
  y_min <- Inf
  y_max <- -Inf

  include <- function(x, y) {
    if (length(x) == 0L || length(y) == 0L) return(invisible(NULL))
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
      x_min <<- min(x_min, min(x[ok]))
      x_max <<- max(x_max, max(x[ok]))
      y_min <<- min(y_min, min(y[ok]))
      y_max <<- max(y_max, max(y[ok]))
    }
    invisible(NULL)
  }

  if (!is.null(allRibbon) && nrow(allRibbon) > 0) {
    include(allRibbon$x, allRibbon$y)
  }
  if (!is.null(seqArcs) && length(seqArcs) > 0) {
    for (arc in seqArcs) {
      if (nrow(arc) > 0) include(arc$x, arc$y)
    }
  }
  if (!is.null(gene_arrows) && nrow(gene_arrows) > 0) {
    include(gene_arrows$text_x, gene_arrows$text_y)
  }
  if (show_axis && !is.null(axisLines) && nrow(axisLines) > 0) {
    include(axisLines$x, axisLines$y)
  }
  if (show_axis && !is.null(axisTicks) && nrow(axisTicks) > 0) {
    include(axisTicks$x0, axisTicks$y0)
    include(axisTicks$x1, axisTicks$y1)
    include(axisTicks$label_x, axisTicks$label_y)
  }
  if (!is.null(gene_polys) && nrow(gene_polys) > 0) {
    include(gene_polys$x, gene_polys$y)
  }
  if (!is.null(seq_labels) && nrow(seq_labels) > 0) {
    include(seq_labels$text_x, seq_labels$text_y)
  }

  list(
    x_min = if (is.finite(x_min)) x_min else NA_real_,
    x_max = if (is.finite(x_max)) x_max else NA_real_,
    y_min = if (is.finite(y_min)) y_min else NA_real_,
    y_max = if (is.finite(y_max)) y_max else NA_real_
  )
}

ggchord_device_units_per_inch <- function(x, y,
                                          fallback_inches = 6,
                                          margin_inches = 1.25,
                                          device_inches = NULL) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  x_span <- if (length(x) > 1) diff(range(x)) else 0
  y_span <- if (length(y) > 1) diff(range(y)) else 0
  geometry_span <- max(x_span, y_span, 1)

  if (is.null(device_inches)) {
    # Querying dev.size() on the null device opens R's default device. That
    # side effect could leave Rplots.pdf open after layout-only operations and
    # make a following ggsave() close the wrong device. Use the documented
    # fallback until an actual render device exists.
    device_inches <- if (grDevices::dev.cur() == 1L) {
      c(NA_real_, NA_real_)
    } else tryCatch(
      grDevices::dev.size("in"), error = function(e) c(NA_real_, NA_real_)
    )
  }
  device_short_side <- suppressWarnings(min(device_inches, na.rm = TRUE))
  if (!is.finite(device_short_side) || device_short_side <= 0) {
    usable_inches <- fallback_inches
  } else {
    # Keep a modest allowance for titles and legends, but never replace a
    # genuinely small device with the much larger fallback canvas. The old
    # `< 2` fallback made 4 x 3 inch exports behave as if they were six inches
    # wide and consequently underestimated every text box.
    reserve <- min(margin_inches, device_short_side * 0.4)
    usable_inches <- max(device_short_side - reserve, 0.75)
  }
  geometry_span / usable_inches
}

#' Convert text layers into fixed obstacle rectangles for label repulsion.
#' @noRd
