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

# Register package-specific theme elements with ggplot2. These elements affect
# annotations drawn by ggchord; data-dependent fill and colour remain scales.
.onLoad <- function(libname, pkgname) {
  ggplot2::register_theme_elements(
    ggchord.axis.line = ggplot2::element_line(),
    ggchord.axis.ticks = ggplot2::element_line(),
    ggchord.axis.text = ggplot2::element_text(),
    ggchord.seq.label = ggplot2::element_text(),
    ggchord.group.label = ggplot2::element_text(),
    ggchord.gene.label = ggplot2::element_text(),
    ggchord.gene.label.segment = ggplot2::element_line(),
    element_tree = list(
      "ggchord.axis.line" = ggplot2::el_def("element_line", inherit = "line"),
      "ggchord.axis.ticks" = ggplot2::el_def("element_line", inherit = "line"),
      "ggchord.axis.text" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.seq.label" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.group.label" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.gene.label" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.gene.label.segment" = ggplot2::el_def(
        "element_line", inherit = "line"
      )
    )
  )
}

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

#' Emit one migration warning per old argument in a session
#' @noRd
ggchord_deprecate_once <- function(argument, replacement) {
  key <- paste0("deprecated:", argument)
  if (!isTRUE(.chord_env[[key]])) {
    warning(
      "`", argument, "` is deprecated for v0.9.0; use `", replacement,
      "` instead.", call. = FALSE
    )
    .chord_env[[key]] <- TRUE
  }
  invisible(NULL)
}

#' Record an old scale argument on a layer for build-time migration checks
#' @noRd
ggchord_add_legacy_scale <- function(
    lyr, supplied, argument, aesthetic, replacement) {
  if (!isTRUE(supplied)) return(lyr)
  spec <- data.frame(
    argument = argument, aesthetic = aesthetic, replacement = replacement,
    stringsAsFactors = FALSE
  )
  lyr$ggchord_legacy_scales <- rbind(lyr$ggchord_legacy_scales, spec)
  lyr
}

#' Check legacy scale arguments and emit their one-time migration warning
#' @noRd
ggchord_check_legacy_scales <- function(plot) {
  for (lyr in plot$layers) {
    specs <- lyr$ggchord_legacy_scales
    if (is.null(specs) || nrow(specs) == 0) next
    for (i in seq_len(nrow(specs))) {
      if (plot$scales$has_scale(specs$aesthetic[i])) {
        ggchord_stop(
          "`", specs$argument[i], "` conflicts with a user-supplied scale for `",
          specs$aesthetic[i], "`; remove the old argument and use `",
          specs$replacement[i], "`"
        )
      }
      ggchord_deprecate_once(specs$argument[i], specs$replacement[i])
    }
  }
  invisible(plot)
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
#' @param plot Optional ggchord plot. Supplying the plot is the reliable way to
#'   retrieve its own layout when several plots are built in one session.
#' @param build Logical. Build the layout when it is not cached, default TRUE.
#' @return A chord layout list containing the computed geometry (sequence
#'   arcs, ribbon polygons, gene arrows, axis elements, extremes, colors, etc.)
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' p <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
#' invisible(ggplot2::ggplot_build(p))
#' names(get_chord_layout(p)$seq_arcs)
get_chord_layout <- function(plot = NULL, build = TRUE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!is.logical(build) || length(build) != 1 || is.na(build)) {
    ggchord_stop("get_chord_layout(): build must be TRUE or FALSE")
  }
  if (!is.null(plot)) {
    if (!inherits(plot, "ggchord") || is.null(plot$ggchord)) {
      ggchord_stop("get_chord_layout(): plot must be a ggchord object")
    }
    layout <- plot$ggchord$ref$layout
    if (is.null(layout) && isTRUE(build)) layout <- compute_chord_geometry(plot)
  } else {
    if (!isTRUE(.chord_env$warned_get_layout)) {
      warning(
        "get_chord_layout() without a plot is deprecated; use get_chord_layout(plot)",
        call. = FALSE
      )
      .chord_env$warned_get_layout <- TRUE
    }
    layout <- .chord_env$layout
  }
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

#' Capture user data and mappings separately from the computed placeholder
#' @noRd
ggchord_capture_layer_input <- function(lyr, data, mapping, roles) {
  lyr$ggchord_input_data <- data
  lyr$ggchord_input_mapping <- mapping
  lyr$ggchord_role_aes <- roles
  lyr
}

#' Resolve a layer's input table and role mappings
#' @noRd
ggchord_resolve_layer_input <- function(lyr, fallback = NULL) {
  data <- lyr$ggchord_input_data %||% fallback
  if (is.null(data)) return(NULL)
  if (!is.data.frame(data)) {
    ggchord_stop("ggchord layer data must be a data.frame")
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  mapping <- lyr$ggchord_input_mapping
  roles <- intersect(lyr$ggchord_role_aes %||% character(0), names(mapping))
  for (role in roles) {
    value <- tryCatch(
      rlang::eval_tidy(mapping[[role]], data = data),
      error = function(e) ggchord_stop(
        "Cannot evaluate `", role, "` in layer mapping: ", conditionMessage(e)
      )
    )
    if (length(value) == 1 && nrow(data) != 1) value <- rep(value, nrow(data))
    if (length(value) != nrow(data)) {
      ggchord_stop("Mapped `", role, "` must return one value per input row")
    }
    data[[role]] <- value
  }
  data
}

#' Combine computed and user visual mappings
#' @noRd
ggchord_effective_mapping <- function(lyr) {
  out <- lyr$mapping
  user <- lyr$ggchord_input_mapping
  if (is.null(user)) return(out)
  visual <- setdiff(names(user), lyr$ggchord_role_aes %||% character(0))
  for (nm in visual) out[[nm]] <- user[[nm]]
  out
}

#' Add original input columns to expanded geometry
#' @noRd
ggchord_attach_input_columns <- function(geometry, input) {
  if (is.null(input) || !is.data.frame(input) || nrow(geometry) == 0) {
    return(geometry)
  }
  protected <- c(
    "x", "y", "x0", "y0", "x1", "y1", "xend", "yend", "group",
    "text_x", "text_y", "text", "label_x", "label_y", "label",
    "text_angle", "label_angle", "hjust", "vjust", "size", "alpha",
    "colour", "fill", "zfill", "zcolour", "zregionfill", "zoutline",
    "zlinetype", "outline_col", "linetype_val", "seq_colour",
    "group_colour", "ribbon_fill", "ribbon_alpha", "ribbon_colour",
    "ribbon_linetype", "gene_fill", "feature_fill", "region_fill"
  )
  if ("source_row" %in% names(geometry)) {
    idx <- geometry$source_row
  } else if ("seq_id" %in% names(geometry) && "seq_id" %in% names(input)) {
    idx <- match(as.character(geometry$seq_id), as.character(input$seq_id))
    geometry <- geometry[!is.na(idx), , drop = FALSE]
    idx <- idx[!is.na(idx)]
  } else {
    return(geometry)
  }
  valid <- !is.na(idx) & idx >= 1 & idx <= nrow(input)
  geometry <- geometry[valid, , drop = FALSE]
  idx <- idx[valid]
  cols <- setdiff(names(input), protected)
  for (nm in cols) geometry[[nm]] <- input[[nm]][idx]
  geometry
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
  registry <- layout$layer_geometry[[lyr$ggchord_layer_id %||% ""]]
  if (!is.null(registry) && !is.null(registry[[lyr$ggchord_type]])) {
    geometry <- registry[[lyr$ggchord_type]]
  } else {
    geometry <- switch(lyr$ggchord_type %||% "",
    seq       = if (length(layout$seq_arcs) > 0) do.call(rbind, layout$seq_arcs) else fallback,
    ribbon    = if (!is.null(layout$ribbon_polys)) layout$ribbon_polys else fallback,
    gene_poly = if (nrow(layout$gene_polys) > 0) layout$gene_polys else fallback,
    gene_text = if (nrow(layout$gene_labels) > 0) layout$gene_labels else fallback,
    gene_text_repel = if (nrow(layout$gene_labels) > 0) layout$gene_labels else fallback,
    gene_label_segment = if (nrow(layout$gene_label_segments) > 0) layout$gene_label_segments else fallback,
    seq_label = if (nrow(layout$seq_labels_df) > 0) layout$seq_labels_df else fallback,
    seq_group_label = if (nrow(layout$group_labels) > 0) layout$group_labels else fallback,
    seq_region = if (nrow(layout$region_polys) > 0) layout$region_polys else fallback,
    ribbon_highlight = if (nrow(layout$ribbon_highlight_polys) > 0) layout$ribbon_highlight_polys else fallback,
    axis_line = if (nrow(layout$axis_lines) > 0) layout$axis_lines else fallback,
    axis_seg  = if (nrow(layout$axis_ticks) > 0) layout$axis_ticks else fallback,
    axis_text = {
      d <- layout$axis_ticks
      if (nrow(d) > 0) d[!is.na(d$label), , drop = FALSE] else fallback
    },
    fallback)
  }
  input <- layout$layer_inputs[[lyr$ggchord_layer_id %||% ""]][[
    lyr$ggchord_type %||% ""
  ]]
  input <- input %||% ggchord_resolve_layer_input(lyr)
  ggchord_attach_input_columns(geometry, input)
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
