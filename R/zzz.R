# zzz.R - package environment and infrastructure
# The package keeps no global state that affects rendering: plot data and
# parameters are stored on the plot object itself. The package environment
# only holds a cache of the most recently computed layout so that the
# get_chord_layout() accessor can inspect it after rendering.

#' Package-level environment
#'
#' Internal environment that caches the most recently computed chord layout.
#'
#' @keywords internal
.chord_env <- new.env(parent = emptyenv())

# ====================================================================
# Error signalling
# ====================================================================

#' Signal an error without invoking a custom global error handler
#'
#' ggchord errors must always print a clear message and return control to the
#' user; they must never drop into an interactive debugger (for example when
#' `options(error = browser)` or RStudio's "Break in Code" error handler is
#' active). This helper temporarily restores the base error handler while the
#' error is signalled, then restores the user's handler afterwards.
#'
#' @param ... Message parts, passed through to [stop()].
#' @param call. Logical. Whether to include the call in the error message.
#' @noRd
ggchord_stop <- function(..., call. = FALSE) {
  old <- getOption("error")
  on.exit(options(error = old), add = TRUE)
  options(error = NULL)
  stop(..., call. = call.)
}

# Disable any custom global error handler (options(error = browser), RStudio's
# "Break in Code") for the duration of the calling function. Exported functions
# call this at the top of their body so that even R-generated errors (e.g. a
# missing required argument) print a plain message instead of dropping into an
# interactive debugger. The error itself still propagates normally.
ggchord_disable_debug <- function() {
  old <- getOption("error")
  options(error = NULL)
  old
}


# ====================================================================
# Layout cache (set at build time; used by the get_chord_layout() accessor)
# ====================================================================

#' Set the chord layout into the package environment
#' @keywords internal
set_chord_layout <- function(layout) {
  .chord_env$layout <- layout
}

#' Get the chord layout from the package environment
#'
#' Returns the most recently computed chord layout (after the plot was built,
#' e.g. via \code{print()} or \code{ggplot_build()}). This is useful for
#' building custom layers or annotations on top of the chord geometry.
#'
#' @return A chord layout list containing the computed geometry (sequence
#'   arcs, ribbon polygons, gene arrows, axis elements, extremes, colors, etc.)
#' @export
get_chord_layout <- function() {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  layout <- .chord_env$layout
  if (is.null(layout)) {
    ggchord_stop(
      "Chord layout data not found. Please render the plot first.",
      call. = FALSE
    )
  }
  layout
}

# ====================================================================
# Environment cleanup
# ====================================================================

#' Clear the package environment (used to reset state)
#' @keywords internal
clear_chord_env <- function() {
  rm(list = ls(.chord_env, all.names = TRUE), envir = .chord_env)
}

# ====================================================================
# Lazy layer data (plotly::ggplotly and other tools call layer$layer_data()
# directly, bypassing ggplot_build()).  Each ggchord layer is tagged with a
# ggchord_type and given a lazy data function that computes (on demand) and
# returns the geometry for that layer.
# ====================================================================

#' Attach the shared plot reference and a lazy data function to a ggchord layer
#' @keywords internal
wire_ggchord_layer <- function(lyr, plot) {
  if (!inherits(lyr, "LayerInstance") || is.null(lyr$ggchord_type)) return(lyr)
  if (is.null(lyr$ggchord_ref) && !is.null(plot$ggchord$ref)) {
    lyr$ggchord_placeholder <- lyr$data
    lyr$ggchord_ref <- plot$ggchord$ref
    lyr$data <- make_ggchord_lazy_data(lyr)
  }
  lyr
}

#' Build a lazy data function for a ggchord layer
#' @keywords internal
make_ggchord_lazy_data <- function(lyr) {
  force(lyr)
  function(plot_data) ggchord_layer_data(lyr)
}

#' Return the computed geometry for a ggchord layer, computing the layout on
#' demand if it has not been computed yet.
#' @keywords internal
ggchord_layer_data <- function(lyr) {
  ref <- lyr$ggchord_ref
  plot <- ref$plot
  layout <- ref$layout
  if (is.null(layout)) {
    # Compute the layout (and, for tools such as plotly::ggplotly() that read
    # layer data directly, attach scales and coordinates to the plot in place).
    plot <- prepare_ggchord_plot(plot)
    layout <- plot$ggchord$layout
  }
  extract_ggchord_layer_data(lyr, layout)
}

#' Extract the geometry for one layer from a computed layout
#' @keywords internal
extract_ggchord_layer_data <- function(lyr, layout) {
  fallback <- lyr$ggchord_placeholder
  switch(lyr$ggchord_type %||% "",
    seq       = if (length(layout$seq_arcs) > 0) do.call(rbind, layout$seq_arcs) else fallback,
    ribbon    = if (!is.null(layout$ribbon_polys)) layout$ribbon_polys else fallback,
    gene_poly = if (nrow(layout$gene_polys) > 0) layout$gene_polys else fallback,
    gene_text = if (nrow(layout$gene_labels) > 0) layout$gene_labels else fallback,
    gene_text_repel = if (nrow(layout$gene_labels) > 0) layout$gene_labels else fallback,
    gene_label_segment = if (nrow(layout$gene_label_segments) > 0) layout$gene_label_segments else fallback,
    seq_label = if (nrow(layout$seq_labels_df) > 0) layout$seq_labels_df else fallback,
    seq_group_label = if (nrow(layout$group_labels) > 0) layout$group_labels else fallback,
    axis_line = if (nrow(layout$axis_lines) > 0) layout$axis_lines else fallback,
    axis_seg  = if (nrow(layout$axis_ticks) > 0) layout$axis_ticks else fallback,
    axis_text = {
      d <- layout$axis_ticks
      if (nrow(d) > 0) d[!is.na(d$label), , drop = FALSE] else fallback
    },
    fallback
  )
}

# ====================================================================
# Helper operators
# ====================================================================

#' NULL coalescing operator
#'
#' Returns \code{y} if \code{x} is NULL, otherwise returns \code{x}.
#'
#' @param x Any R object (may be NULL)
#' @param y Default value returned when \code{x} is NULL
#' @name null-coalescing-operator
#' @rdname null-coalescing-operator
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
