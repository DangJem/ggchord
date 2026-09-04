# zzz.R - package environment and infrastructure
# The package keeps no global rendering or layout state. Plot data, parameters
# and any computed layout cache belong to the plot object itself.

# Register package-specific theme elements with ggplot2. These elements affect
# annotations drawn by ggchord; data-dependent fill and colour remain scales.
.onLoad <- function(libname, pkgname) {
  ggplot2::register_theme_elements(
    ggchord.axis.line = ggplot2::element_line(),
    ggchord.axis.ticks = ggplot2::element_line(),
    ggchord.axis.minor.ticks = ggplot2::element_line(),
    ggchord.axis.text = ggplot2::element_text(),
    ggchord.seq.label = ggplot2::element_text(),
    ggchord.gene.label = ggplot2::element_text(),
    ggchord.gene.label.segment = ggplot2::element_line(),
    element_tree = list(
      "ggchord.axis.line" = ggplot2::el_def("element_line", inherit = "line"),
      "ggchord.axis.ticks" = ggplot2::el_def("element_line", inherit = "line"),
      "ggchord.axis.minor.ticks" = ggplot2::el_def(
        "element_line", inherit = "line"
      ),
      "ggchord.axis.text" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.seq.label" = ggplot2::el_def("element_text", inherit = "text"),
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

#' Reject removed public arguments with a concise migration message
#' @noRd
ggchord_reject_retired <- function(dots, caller, migrations,
                                   call = sys.call(-1L)) {
  raw <- names(as.list(call)[-1L]) %||% character()
  supplied <- unique(intersect(c(names(dots), raw), names(migrations)))
  if (!length(supplied)) return(invisible(NULL))
  advice <- unique(unname(migrations[supplied]))
  ggchord_stop(
    caller, ": removed argument(s): ", paste(supplied, collapse = ", "),
    ". Use ", paste(advice, collapse = "; "), "."
  )
}

# ====================================================================
# Plot-owned layout access
# ====================================================================

#' Get the chord layout from a ggchord plot
#'
#' Returns the layout owned by an explicit ggchord plot. This avoids ambiguous
#' cross-talk when several plots are built in one R session.
#'
#' @param plot A ggchord plot.
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
get_chord_layout <- function(plot, build = TRUE) {
  old_error <- ggchord_disable_debug()
  on.exit(options(error = old_error), add = TRUE)

  if (!is.logical(build) || length(build) != 1 || is.na(build)) {
    ggchord_stop("get_chord_layout(): build must be TRUE or FALSE")
  }
  if (missing(plot) || !inherits(plot, "ggchord") || is.null(plot$ggchord)) {
    ggchord_stop("get_chord_layout(): plot must be supplied as a ggchord object")
  }
  layout <- plot$ggchord$ref$layout
  if (is.null(layout) && isTRUE(build)) layout <- compute_chord_geometry(plot)
  if (is.null(layout)) {
    ggchord_stop(
      "Chord layout data not found. Please render the plot first.",
      call. = FALSE
    )
  }
  layout
}

#' Capture user data and mappings separately from the computed placeholder
#' @noRd
ggchord_capture_layer_input <- function(lyr, data, mapping, roles) {
  mapping <- ggchord_normalize_mapping(mapping)
  lyr$ggchord_input_data <- data
  lyr$ggchord_input_mapping <- mapping
  lyr$ggchord_role_aes <- roles
  lyr
}

#' Normalise public US spelling aliases before geometry or scale inference
#' @noRd
ggchord_normalize_mapping <- function(mapping) {
  if (is.null(mapping)) return(mapping)
  if (anyDuplicated(names(mapping))) {
    duplicated_names <- unique(names(mapping)[duplicated(names(mapping))])
    ggchord_stop(
      "Duplicated aesthetic after colour/color normalisation: ",
      paste(duplicated_names, collapse = ", ")
    )
  }
  aliases <- c(seq_color = "seq_colour", ribbon_color = "ribbon_colour")
  for (alias in names(aliases)) {
    canonical <- aliases[[alias]]
    if (alias %in% names(mapping) && canonical %in% names(mapping)) {
      ggchord_stop("Use only one of `", canonical, "` and `", alias, "`")
    }
    if (alias %in% names(mapping)) names(mapping)[names(mapping) == alias] <- canonical
  }
  mapping
}

#' Resolve a fixed colour/color alias without ambiguous precedence
#' @noRd
ggchord_colour_alias <- function(colour, dots, caller, call = sys.call(-1L)) {
  raw <- names(as.list(call)[-1L]) %||% character()
  has_color <- "color" %in% names(dots)
  has_colour <- "colour" %in% raw || "colour" %in% names(dots)
  if (has_color && has_colour) {
    ggchord_stop(caller, ": use only one of colour and color")
  }
  if (has_color) {
    colour <- dots$color
    dots$color <- NULL
  }
  if ("color" %in% names(dots)) names(dots)[names(dots) == "color"] <- "colour"
  list(colour = colour, dots = dots)
}

#' Normalise colour/color inside variadic fixed aesthetics
#' @noRd
ggchord_colour_dots <- function(dots, caller) {
  if ("color" %in% names(dots) && "colour" %in% names(dots))
    ggchord_stop(caller, ": use only one of colour and color")
  if ("color" %in% names(dots)) names(dots)[names(dots) == "color"] <- "colour"
  dots
}

#' Resolve a layer's input table and role mappings
#' @noRd
ggchord_resolve_layer_input <- function(lyr, fallback = NULL) {
  data <- lyr$ggchord_input_data %||% fallback
  if (is.null(data)) return(NULL)
  if (is.function(data)) {
    data_function <- data
    data <- tryCatch(
      data_function(fallback),
      error = function(e) ggchord_stop(
        "ggchord layer data function failed: ", conditionMessage(e)
      )
    )
  }
  if (!is.data.frame(data)) {
    ggchord_stop("ggchord layer data must be a data.frame or a function returning one")
  }
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  mapping <- lyr$ggchord_input_mapping
  roles <- intersect(lyr$ggchord_role_aes %||% character(0), names(mapping))
  for (role in roles) {
    expr <- rlang::quo_get_expr(mapping[[role]])
    contains_stage <- function(x) {
      if (!is.call(x)) return(FALSE)
      if (identical(as.character(x[[1L]]), "after_stat") ||
          identical(as.character(x[[1L]]), "after_scale")) return(TRUE)
      any(vapply(as.list(x)[-1L], contains_stage, logical(1)))
    }
    if (contains_stage(expr)) {
      ggchord_stop(
        "Role aesthetic `", role,
        "` controls geometry and cannot use after_stat() or after_scale()"
      )
    }
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
  transform <- lyr$ggchord_input_transform
  if (!is.null(transform)) data <- transform(data)
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
    "ribbon_fill", "ribbon_alpha", "ribbon_colour",
    "ribbon_linetype", "gene_fill", "feature_fill", "feature_shape",
    "region_fill"
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
    seq_region = if (nrow(layout$region_polys) > 0) layout$region_polys else fallback,
    ribbon_highlight = if (nrow(layout$ribbon_highlight_polys) > 0) layout$ribbon_highlight_polys else fallback,
    axis_line = if (nrow(layout$axis_lines) > 0) layout$axis_lines else fallback,
    axis_seg  = if (nrow(layout$axis_ticks) > 0) layout$axis_ticks else fallback,
    axis_text = {
      d <- layout$axis_ticks
      if (nrow(d) > 0) d[!is.na(d$label), , drop = FALSE] else fallback
    },
    axis = ggchord_axis_geometry(layout),
    gene_label_repel = ggchord_repel_geometry(layout),
    fallback)
  }
  input <- layout$layer_inputs[[lyr$ggchord_layer_id %||% ""]][[
    lyr$ggchord_type %||% ""
  ]]
  input <- input %||% ggchord_resolve_layer_input(lyr)
  ggchord_attach_input_columns(geometry, input)
}

#' Combine axis components for GeomChordAxis
#' @noRd
ggchord_axis_geometry <- function(layout) {
  line <- layout$axis_lines %||% data.frame()
  tick <- layout$axis_ticks %||% data.frame()
  label <- if (nrow(tick) && "label" %in% names(tick)) {
    tick[!is.na(tick$label), , drop = FALSE]
  } else {
    tick[integer(0), , drop = FALSE]
  }
  if (nrow(line)) {
    line$.component <- "line"
    line$group <- line$seq_id
  }
  if (nrow(tick)) {
    tick$.component <- ifelse(tick$is_major, "major_tick", "minor_tick")
    tick$x <- tick$x0
    tick$y <- tick$y0
    tick$xend <- tick$x1
    tick$yend <- tick$y1
  }
  if (nrow(label)) {
    label$.component <- "text"
    label$x <- label$label_x
    label$y <- label$label_y
    label$angle <- label$label_angle
    label$hjust <- label$label_hjust
    label$vjust <- label$label_vjust
  }
  values <- Filter(nrow, list(line, tick, label))
  if (!length(values)) {
    return(data.frame(x = numeric(), y = numeric(), .component = character()))
  }
  ggchord_rbind_fill(values)
}

#' Internal automatic sequence-axis layer
#' @noRd
ggchord_axis_layer <- function(layout) {
  data <- ggchord_axis_geometry(layout)
  lyr <- ggplot2::layer(
    data = data,
    mapping = ggplot2::aes(
      x = x, y = y, xend = xend, yend = yend, label = label,
      group = group, size = I(size), angle = angle,
      hjust = hjust, vjust = vjust, .component = I(.component)
    ),
    stat = "identity", geom = GeomChordAxis, position = "identity",
    show.legend = FALSE, inherit.aes = FALSE, check.aes = FALSE,
    check.param = FALSE,
    params = list(line_params = list(), tick_params = list(),
      minor_tick_params = list(), text_params = list(), na.rm = FALSE)
  )
  lyr$ggchord_theme_components <- c(
    line_params = "ggchord.axis.line",
    tick_params = "ggchord.axis.ticks",
    minor_tick_params = "ggchord.axis.minor.ticks",
    text_params = "ggchord.axis.text"
  )
  lyr
}

#' Combine automatic label components for GeomChordGeneLabelRepel
#' @noRd
ggchord_repel_geometry <- function(layout) {
  segment <- layout$gene_label_segments %||% data.frame()
  text <- layout$gene_labels %||% data.frame()
  if (nrow(segment)) {
    # Leader paths are stored as one or more segments per label and therefore
    # only carry the label group.  Restore the source-row identity before the
    # geometry is joined to user columns; otherwise the generic join drops all
    # segment rows while retaining the text rows.
    if (nrow(text) && "group" %in% names(segment) &&
        "group" %in% names(text)) {
      label_index <- match(segment$group, text$group)
      for (nm in intersect(c("source_row", "seq_id"), names(text))) {
        segment[[nm]] <- text[[nm]][label_index]
      }
    }
    segment$.component <- "segment"
    segment$x <- segment$x0
    segment$y <- segment$y0
    segment$xend <- segment$x1
    segment$yend <- segment$y1
  }
  if (nrow(text)) {
    text$.component <- "text"
    text$x <- text$text_x
    text$y <- text$text_y
    text$label <- text$text
    text$angle <- text$text_angle
  }
  values <- Filter(nrow, list(segment, text))
  if (!length(values)) {
    return(data.frame(x = numeric(), y = numeric(), .component = character()))
  }
  ggchord_rbind_fill(values)
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
