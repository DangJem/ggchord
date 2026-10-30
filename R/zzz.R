# zzz.R - package environment and infrastructure
# Plot data, parameters and persistent layout caches belong to each plot.
# Temporary render measurement context is scoped and restored in render-context.R.

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
    ggchord.feature.label = ggplot2::element_text(),
    ggchord.feature.label.segment = ggplot2::element_line(),
    ggchord.restriction.label = ggplot2::element_text(),
    ggchord.restriction.label.segment = ggplot2::element_line(),
    ggchord.seq.center.label = ggplot2::element_text(),
    ggchord.seq.center.name = ggplot2::element_text(),
    ggchord.seq.center.length = ggplot2::element_text(),
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
      ),
      "ggchord.feature.label" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.feature.label.segment" = ggplot2::el_def("element_line", inherit = "line"),
      "ggchord.restriction.label" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.restriction.label.segment" = ggplot2::el_def("element_line", inherit = "line"),
      "ggchord.seq.center.label" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.seq.center.name" = ggplot2::el_def("element_text", inherit = "text"),
      "ggchord.seq.center.length" = ggplot2::el_def("element_text", inherit = "text")
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
#' p <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_link_ribbon()
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
  if (ggchord_physical_requested(plot) &&
      !identical(plot$ggchord$ref$context_key, ggchord_render_cache_key()))
    plot$ggchord$ref$layout <- NULL
  layout <- plot$ggchord$ref$layout
  if (is.null(layout) && isTRUE(build)) {
    branching <- any(vapply(plot$layers, function(l) (l$ggchord_params$link_branch %||% "none") != "none", logical(1)))
    if (branching) {
      invisible(ggplot2::ggplot_build(plot))
      layout <- plot$ggchord$ref$layout
    } else layout <- compute_chord_geometry(plot)
  }
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
  lyr$ggchord_input_data <- ggchord_normalize_accver(data)
  lyr$ggchord_input_mapping <- mapping
  lyr$ggchord_role_aes <- roles
  lyr
}

#' Normalise public US spelling aliases before geometry or scale inference
#' @noRd
ggchord_normalize_mapping <- function(mapping) {
  if (is.null(mapping)) return(mapping)
  if ("seq_id" %in% names(mapping)) {
    ggchord_stop("`seq_id` mapping was removed in v0.13.1; use `accver`")
  }
  if (anyDuplicated(names(mapping))) {
    duplicated_names <- unique(names(mapping)[duplicated(names(mapping))])
    ggchord_stop(
      "Duplicated aesthetic after colour/color normalisation: ",
      paste(duplicated_names, collapse = ", ")
    )
  }
  aliases <- c(seq_color = "seq_colour", ribbon_color = "ribbon_colour", link_color = "link_colour")
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
  type <- lyr$ggchord_params$type %||% ""
  if (inherits(data, "GRanges") || inherits(data, "Seqinfo")) {
    data <- if (type %in% c("seq", "seq_label", "seq_break", "seq_center_label")) {
      as_ggchord_seq(data)
    } else if (inherits(data, "GRanges")) as_ggchord_feature(data) else {
      ggchord_stop("Seqinfo supplies sequence lengths; use GRanges for feature layers")
    }
  } else if (inherits(data, "GInteractions")) {
    data <- as_ggchord_link(data)
    if (identical(type, "link")) {
      data$qpos <- (data$qstart + data$qend) / 2
      data$spos <- (data$sstart + data$send) / 2
      data$.point_projection <- "anchor_midpoint"
    }
  }
  if (!is.data.frame(data)) {
    ggchord_stop("ggchord layer data must be a data.frame, a supported interval object, ",
      "or a function returning supported data")
  }
  data <- ggchord_normalize_accver(as.data.frame(data, stringsAsFactors = FALSE))
  mapping <- lyr$ggchord_input_mapping
  roles <- intersect(lyr$ggchord_role_aes %||% character(0), names(mapping))
  for (role in roles) {
    expr <- rlang::get_expr(mapping[[role]])
    contains_stage <- function(x) {
      if (!is.call(x)) return(FALSE)
      head <- x[[1L]]
      name <- if (is.symbol(head)) as.character(head) else if (
        is.call(head) && as.character(head[[1L]]) %in% c("::", ":::")) {
        as.character(head[[3L]])
      } else ""
      if (name %in% c("after_stat", "after_scale")) return(TRUE)
      if (identical(name, "stage")) {
        args <- as.list(x)[-1L]
        late <- args[intersect(names(args), c("after_stat", "after_scale"))]
        if (any(!vapply(late, is.null, logical(1)))) return(TRUE)
      }
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
  for (nm in setdiff(names(mapping), c(roles, names(lyr$aes_params)))) {
    expr <- rlang::get_expr(mapping[[nm]])
    text <- rlang::as_label(expr)
    if (!grepl("after_stat|after_scale", text)) {
      columns <- intersect(all.vars(expr), c("pident", "length", "anno", "type", "label"))
      ggchord_require_columns(data, columns, paste0("Mapping `", nm, "`"))
    }
  }
  transform <- lyr$ggchord_input_transform
  if (!is.null(transform)) data <- transform(data)
  if (identical(type, "gene") && !isTRUE(lyr$ggchord_params$is_feature) &&
      any(!data$strand %in% c("+", "-"))) ggchord_stop(
    "geom_gene(): directional genes need '+'/'-' strand; use geom_feature() for nondirectional intervals"
  )
  data
}

#' Resolve omitted fixed styles against local and inherited user mappings
#' @noRd
ggchord_resolve_style_defaults <- function(layer, mapping) {
  defaults <- layer$ggchord_params$style_defaults
  if (is.null(defaults)) return(layer)
  user <- layer$ggchord_input_mapping
  if (isTRUE(layer$inherit.aes)) {
    inherited <- ggchord_normalize_mapping(mapping)
    for (nm in setdiff(names(inherited), names(user))) user[[nm]] <- inherited[[nm]]
  }
  layer <- reconstruct_layer(layer, layer$data)
  # Ribbon outlines and line types use role names only when that role is mapped.
  if (identical(layer$ggchord_type, "ribbon")) {
    for (nm in c("colour", "linetype")) {
      role <- paste0("ribbon_", nm)
      if (role %in% names(user) && nm %in% names(layer$aes_params)) {
        layer$aes_params[role] <- layer$aes_params[nm]
        layer$aes_params[nm] <- NULL
        names(defaults)[names(defaults) == nm] <- role
      }
    }
    old_key <- layer$geom$draw_key
    layer$geom <- make_ribbon_geom(
      outline = "ribbon_colour" %in% names(user),
      linetype = "ribbon_linetype" %in% names(user))
    layer$geom$draw_key <- old_key
  }
  mapped_defaults <- intersect(names(defaults)[defaults], names(user))
  layer$aes_params[mapped_defaults] <- NULL
  layer
}

#' Combine computed and user visual mappings
#' @noRd
ggchord_effective_mapping <- function(lyr) {
  out <- lyr$mapping
  user <- lyr$ggchord_input_mapping
  input <- lyr$ggchord_resolved_input
  if (identical(lyr$ggchord_type, "ribbon") && !is.null(input) &&
      !"pident" %in% names(input) && !"ribbon_fill" %in% names(user)) {
    out$ribbon_fill <- NULL
  }
  if (lyr$ggchord_type %in% c("gene_poly") && !is.null(input) &&
      (!"anno" %in% names(input) || all(is.na(input$anno)))) {
    if (!"gene_fill" %in% names(user)) out$gene_fill <- NULL
    if (!"feature_fill" %in% names(user)) out$feature_fill <- NULL
  }
  if (!is.null(user)) {
    visual <- setdiff(names(user), lyr$ggchord_role_aes %||% character(0))
    for (nm in visual) out[[nm]] <- user[[nm]]
  }
  # Fixed aesthetics bypass scales, including automatically supplied mappings.
  out[intersect(names(out), names(lyr$aes_params))] <- NULL
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
    "feature_label_colour", "feature_label_orientation", "region_fill", ".component",
    "track_id", "label_text_original", "label_size_requested", "label_size_final",
    "label_size_min", "label_fit_reason", ".ggchord_resolved_size"
  )
  if ("source_row" %in% names(geometry)) {
    idx <- geometry$source_row
  } else if ("accver" %in% names(geometry) && "accver" %in% names(input)) {
    idx <- match(as.character(geometry$accver), as.character(input$accver))
    geometry <- geometry[!is.na(idx), , drop = FALSE]
    idx <- idx[!is.na(idx)]
  } else {
    return(geometry)
  }
  valid <- !is.na(idx) & idx >= 1 & idx <= nrow(input)
  # Protect computed geometry, not names that the geometry does not contain.
  # A source column such as label may still supply an ordinary visual mapping.
  cols <- setdiff(names(input), intersect(protected, names(geometry)))
  for (nm in cols) {
    source <- input[[nm]]
    if (is.list(source) && !is.data.frame(source)) {
      value <- rep(list(NULL), nrow(geometry))
      value[valid] <- source[idx[valid]]
    } else {
      value <- rep(ggchord_typed_na(source), nrow(geometry))
      value[valid] <- source[idx[valid]]
    }
    geometry[[nm]] <- value
  }
  geometry
}

#' Extract the geometry for one layer from a computed layout
#' @noRd
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
    seq_break = layout$seq_breaks %||% fallback,
    seq_region = if (nrow(layout$region_polys) > 0) layout$region_polys else fallback,
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
    line$group <- line$accver
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
  mapping <- ggplot2::aes(
    x = x, y = y, xend = xend, yend = yend, label = label,
    group = group, size = I(size), angle = angle,
    hjust = hjust, vjust = vjust, .component = I(.component)
  )
  # A line-only axis has no tick endpoints or text columns. Let the geom's
  # existing defaults supply unused aesthetics rather than evaluating absent
  # columns when breaks or labels were explicitly disabled.
  mapping <- mapping[names(mapping) %in% names(data)]
  lyr <- ggplot2::layer(
    data = data,
    mapping = mapping,
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
ggchord_feature_arc_text <- function(text, units_per_inch) {
  if (!nrow(text)) return(data.frame())
  # Every internal feature label follows the radius selected by the track
  # solver. This includes labels moved off the feature into an inner track.
  eligible <- text$feature_label_mode %in% c("inside", "adjacent") &
    !grepl("\n", text$text, fixed = TRUE) &
    !is.na(text$text) & nzchar(text$text)
  text$.draw_as_arc <- eligible
  curved <- ggchord_arc_text_layout(
    text[eligible, , drop = FALSE], units_per_inch = units_per_inch
  )
  if (nrow(curved$glyphs)) curved$glyphs$.component <- "arc_text"
  list(text = text, glyphs = curved$glyphs)
}

#' @noRd
ggchord_repel_geometry <- function(layout) {
  segment <- layout$gene_label_segments %||% data.frame()
  text <- layout$gene_labels %||% data.frame()
  arc_text <- data.frame()
  if (nrow(text) && identical(layout$gene_label_layout, "feature")) {
    curved <- ggchord_feature_arc_text(
      text, layout$text_units_per_inch %||% .35
    )
    text <- curved$text
    arc_text <- curved$glyphs
  }
  if (nrow(text) && !".draw_as_arc" %in% names(text)) {
    text$.draw_as_arc <- FALSE
  }
  if (nrow(segment)) {
    if ("occluded" %in% names(segment)) {
      segment <- segment[
        order(!segment$occluded, seq_len(nrow(segment))), , drop = FALSE
      ]
    }
    # Leader paths are stored as one or more segments per label and therefore
    # only carry the label group.  Restore the source-row identity before the
    # geometry is joined to user columns; otherwise the generic join drops all
    # segment rows while retaining the text rows.
    if (nrow(text) && "group" %in% names(segment) &&
        "group" %in% names(text)) {
      label_index <- match(segment$group, text$group)
      for (nm in intersect(
          c("source_row", "accver", "feature_track", "label_track"),
          names(text))) {
        segment[[nm]] <- text[[nm]][label_index]
      }
      if (all(c("feature_track", "label_track") %in% names(segment))) {
        segment$leader_track_start <- segment$feature_track + 1L
        segment$leader_track_end <- segment$label_track - 1L
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
    text$xend <- NA_real_
    text$yend <- NA_real_
    text$label <- text$text
    text$angle <- text$text_angle
    text$linetype <- "solid"
    text$alpha <- 1
  }
  if (nrow(arc_text)) {
    arc_text$x <- arc_text$text_x
    arc_text$y <- arc_text$text_y
    arc_text$xend <- NA_real_
    arc_text$yend <- NA_real_
    arc_text$label <- arc_text$text
    arc_text$angle <- arc_text$text_angle
    arc_text$linetype <- "solid"
    arc_text$alpha <- 1
  }
  values <- Filter(nrow, list(segment, text, arc_text))
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
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
