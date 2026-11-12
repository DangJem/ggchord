# Evidence-limited checks of one requested render. These do not certify file
# portability, font embedding, accessibility, or a user's biological choices.

ggchord_diagnostic_device <- function(path, device, width, height, dpi) {
  if (device == "png") {
    args <- list(filename = path, width = width, height = height, units = "in", res = dpi)
    if (isTRUE(capabilities("aqua"))) args$type <- "quartz"
    do.call(grDevices::png, args)
  } else if (device == "pdf") {
    grDevices::pdf(path, width = width, height = height)
  } else if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(path, width = width, height = height)
  } else grDevices::svg(path, width = width, height = height)
}

ggchord_diagnostic_guides <- function(table, canvas_mm) {
  selected <- which(grepl("^guide-box", table$layout$name))
  values <- lapply(selected, function(i) {
    grob <- table$grobs[[i]]
    if (inherits(grob, "zeroGrob")) return(NULL)
    size <- tryCatch(c(width = grid::convertWidth(grid::grobWidth(grob), "mm", valueOnly = TRUE),
      height = grid::convertHeight(grid::grobHeight(grob), "mm", valueOnly = TRUE)),
      error = function(e) c(width = NA_real_, height = NA_real_))
    data.frame(guide = table$layout$name[i], width_mm = size[["width"]],
      height_mm = size[["height"]],
      canvas_area_fraction = prod(size) / prod(canvas_mm),
      exceeds_canvas = if (all(is.finite(size))) any(size > canvas_mm) else NA,
      stringsAsFactors = FALSE)
  })
  values <- Filter(Negate(is.null), values)
  if (length(values)) do.call(rbind, values) else data.frame(guide = character(),
    width_mm = numeric(), height_mm = numeric(), canvas_area_fraction = numeric(), exceeds_canvas = logical())
}

ggchord_diagnostic_text_styles <- function(data, layer) {
  # Identity layers retain render-order sidecars that are not aesthetics.
  # In particular axis arc glyphs draw as component='text'; losing their
  # parent identity turns intentional contacts within one string into clashes.
  original <- layer$data
  if (is.data.frame(original) && nrow(original) == nrow(data) &&
      inherits(layer$stat, 'StatIdentity') && inherits(layer$position, 'PositionIdentity')) for (field in setdiff(intersect(
        c('.arc_parent_group', '.arc_glyph_index', '.axis_label', 'source_row',
          'accver', 'display_instance_id'), names(original)), names(data))) data[[field]] <- original[[field]]
  if (!'.component' %in% names(data)) return(data)
  apply_style <- function(component, style, overrides = list()) {
    at <- which(data$.component %in% component)
    if (!length(at)) return()
    if (isTRUE(style$hidden)) {
      data$label[at] <<- NA_character_
      return()
    }
    style <- c(style[setdiff(names(style), c('hidden', names(overrides)))], overrides)
    for (key in intersect(names(style), c('colour', 'alpha', 'size', 'family', 'fontface',
        'lineheight', 'angle', 'hjust', 'vjust'))) if (length(style[[key]]) == 1L) {
      if (!key %in% names(data)) data[[key]] <<- rep(NA, nrow(data))
      data[[key]][at] <<- style[[key]]
    }
  }
  if (inherits(layer$geom, 'GeomGgchordTrack')) {
    p <- layer$geom_params
    for (role in names(p$track_keys)) data[[role]] <- data[[p$track_keys[[role]]]]
    overrides <- list(alpha = 1)
    if (!is.null(p$axis_colour)) overrides$colour <- p$axis_colour
    if (!is.null(p$axis_size)) overrides$size <- p$axis_size
    apply_style('track_axis_text', p$axis_text_params, overrides)
    apply_style('track_axis_title', p$axis_title_params, overrides)
  } else if (inherits(layer$geom, 'GeomRestrictionSite')) {
    apply_style('label', layer$geom_params$text_params)
    data$.ggchord_composite_label <- rep(isTRUE(layer$geom_params$composite_labels), nrow(data))
  } else if (inherits(layer$geom, 'GeomChordAxis') || inherits(layer$geom, 'GeomChordGeneLabelRepel')) {
    style <- layer$geom_params$text_params
    if (inherits(layer$geom, 'GeomChordGeneLabelRepel') && isTRUE(style$.ggchord_resolved_size)) style$size <- NULL
    apply_style(c('text', 'arc_text'), style)
  }
  data
}

ggchord_diagnostic_text <- function(data, coordinate, panel, context, parse = FALSE, size_unit = 'mm') {
  if (!all(c("x", "y", "label", "size") %in% names(data))) return(data.frame())
  selected <- !is.na(data$label) & nzchar(as.character(data$label)) &
    is.finite(data$size) & data$size > 0
  if ('alpha' %in% names(data)) selected <- selected & (is.na(data$alpha) | data$alpha > 0)
  if ('colour' %in% names(data)) selected <- selected & vapply(data$colour, function(colour)
    !is.na(colour) && grDevices::col2rgb(colour, alpha = TRUE)[4L] > 0, logical(1))
  if (".component" %in% names(data) && any(!is.na(data$.component))) selected <- selected &
    data$.component %in% c("text", "label", "arc_text", "name", "length", "track_axis_text", "track_axis_title", "track_text")
  if (all(c(".draw_as_arc", ".component") %in% names(data))) selected <- selected &
    !(data$.draw_as_arc %in% TRUE & data$.component %in% "text")
  data <- data[selected, , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  xy <- coordinate$transform(data, panel)
  justify <- utils::getFromNamespace('compute_just', 'ggplot2')
  data$hjust <- justify(xy$hjust %||% .5, xy$x, xy$y, xy$angle %||% 0)
  data$vjust <- justify(xy$vjust %||% .5, xy$y, xy$x, xy$angle %||% 0)
  points_per_unit <- switch(size_unit, mm = 72.27/25.4, cm = 72.27/2.54,
    'in' = 72.27, pt = 1, pc = 12, 72.27/25.4)
  result <- lapply(seq_len(nrow(data)), function(i) {
    value <- function(name, fallback) if (is.null(data[[name]])) fallback else data[[name]][i]
    label <- if (isTRUE(parse) || isTRUE(value("plotmath_label", FALSE))) {
      tryCatch(parse(text = as.character(data$label[i]))[[1L]], error = function(e) as.character(data$label[i]))
    } else data$label[i]
    grob <- grid::textGrob(label, gp = grid::gpar(
      fontsize = data$size[i] * points_per_unit, fontfamily = value("family", ""),
      fontface = value("fontface", 1), lineheight = value("lineheight", 1.2)))
    w <- grid::convertWidth(grid::grobWidth(grob), "mm", valueOnly = TRUE)
    h <- grid::convertHeight(grid::grobHeight(grob), "mm", valueOnly = TRUE)
    composite <- isTRUE(value('.ggchord_composite_label', FALSE)) &&
      all(c('enzyme_label', 'coordinate_label', 'enzyme_fontface') %in% names(data))
    if (composite) {
      exact <- ggchord_site_label_metrics(value('enzyme_label', ''), value('coordinate_label', ''),
        data$size[i], value('family', ''), value('enzyme_fontface', 'plain'),
        value('lineheight', 1.2))
      w <- exact$width; h <- exact$height
    }
    angle <- value("angle", 0) * pi / 180
    hj <- if (composite) .5 else value('hjust', .5); vj <- value('vjust', .5)
    cx <- xy$x[i] * context$panel_mm[["width"]] + context$panel_bounds_mm[["xmin"]] +
      (.5 - hj) * w * cos(angle) - (.5 - vj) * h * sin(angle)
    cy <- xy$y[i] * context$panel_mm[["height"]] + context$panel_bounds_mm[["ymin"]] +
      (.5 - hj) * w * sin(angle) + (.5 - vj) * h * cos(angle)
    bw <- abs(cos(angle)) * w + abs(sin(angle)) * h
    bh <- abs(sin(angle)) * w + abs(cos(angle)) * h
    data.frame(source_row = value("source_row", NA_integer_), label = as.character(data$label[i]),
      text_object = if (value('.component', '') %in% 'arc_text' || !is.na(value('.arc_glyph_index', NA_real_)))
        paste0('arc-', value('.arc_parent_group', value('group', NA_real_)), '/', value('source_row', NA_integer_), '/',
          value('display_instance_id', value('accver', ''))) else paste0('row-', i),
      size_mm = as.numeric(data$size[i]) * points_per_unit / (72.27/25.4), xmin_mm = cx - bw / 2, xmax_mm = cx + bw / 2,
      ymin_mm = cy - bh / 2, ymax_mm = cy + bh / 2,
      centre_x_mm = cx, centre_y_mm = cy, width_mm = w, height_mm = h,
      angle_radians = angle, stringsAsFactors = FALSE)
  })
  unique(do.call(rbind, result))
}

ggchord_diagnostic_colours <- function(built, colour_vision) {
  out <- list()
  for (i in seq_along(built$data)) {
    data <- built$data[[i]]
    roles <- names(data)[grepl("(^|_)(fill|colour|color)(_|$)", names(data))]
    for (role in roles) {
      scale <- built$plot$scales$get_scales(role)
      if (is.null(scale) || !isTRUE(scale$is_discrete())) next
      colours <- unique(as.character(data[[role]]))
      colours <- colours[!is.na(colours) & colours != "transparent"]
      if (length(colours) < 2L || length(colours) > 64L) next
      rgb <- tryCatch(grDevices::col2rgb(colours, alpha = TRUE), error = function(e) NULL)
      if (is.null(rgb)) next
      colours <- colours[rgb[4L, ] == 255]
      if (length(colours) < 2L) next
      model <- list(grayscale = colours)
      if (colour_vision && requireNamespace("colorspace", quietly = TRUE)) {
        model$protan <- colorspace::protan(colours, severity = 1, linear = TRUE)
        model$deutan <- colorspace::deutan(colours, severity = 1, linear = TRUE)
        model$tritan <- colorspace::tritan(colours, severity = 1, linear = TRUE)
      }
      for (name in names(model)) {
        values <- grDevices::col2rgb(model[[name]]) / 255
        if (name == "grayscale") {
          linear <- ifelse(values <= .04045, values / 12.92, ((values + .055) / 1.055)^2.4)
          lum <- colSums(linear * c(.2126, .7152, .0722))
          distance <- abs(outer(lum, lum, `-`)); threshold <- .05
        } else {
          lab <- grDevices::convertColor(t(values), from = "sRGB", to = "Lab", scale.in = 1)
          distance <- as.matrix(stats::dist(lab)); threshold <- 10
        }
        pair <- which(upper.tri(distance) & distance < threshold, arr.ind = TRUE)
        if (!nrow(pair)) next
        out[[length(out) + 1L]] <- data.frame(layer = i, aesthetic = role, model = name,
          colour1 = colours[pair[, 1L]], colour2 = colours[pair[, 2L]],
          distance = distance[pair], threshold = threshold, stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out)) do.call(rbind, out) else data.frame(layer = integer(),
    aesthetic = character(), model = character(), colour1 = character(), colour2 = character(),
    distance = numeric(), threshold = numeric())
}

#' Inspect one target-size render without changing the plot
#'
#' Builds on a temporary device of the requested kind and measures the actual
#' single panel after gtable layout, then refines physical geometry within a
#' bounded render-local context. The original plot and its layout cache are
#' left unchanged. Standard saving continues to use [ggplot2::ggsave()].
#'
#' Results describe evidence, not a general quality certificate. Checks cover
#' registered vertex bounds, measured text, small text/lines, oversized guides,
#' and existing label/track diagnostics. Objects outside a panel with
#' `clip = "off"` are legal when they still fit the canvas. Unregistered custom
#' grobs, arrowhead extents, font embedding/fallback and cross-application
#' transparency are not certified. Device warnings remain in the report.
#'
#' Grayscale suggestions compare luminance of opaque discrete mapped colours
#' (difference below 0.05); continuous gradients are deliberately excluded.
#' With `colour_vision = TRUE`, optional colorspace public simulations for
#' protan/deutan/tritan at severity 1 in linear RGB flag Lab distances below 10.
#' These heuristic thresholds ignore spatial context, alpha compositing and
#' viewing conditions, and never automatically change the plot or fail a run.
#'
#' @param plot A ggchord plot.
#' @param width,height Positive target canvas dimensions.
#' @param units Dimensions in `"in"`, `"cm"`, `"mm"`, or `"px"`.
#' @param device Target device: `"png"`, `"pdf"`, or `"svg"`. It must be
#'   available in the current R environment. Only this device is inspected.
#' @param dpi Raster resolution; also defines pixel dimensions.
#' @param min_size,min_linewidth Advisory minimum text size and line width,
#'   both in millimetres. These do not modify the supplied style.
#' @param colour_vision Request optional simulated colour-vision hints.
#' @return A `ggchord_diagnostics` list with `measurement`, `checks`, text and
#'   guide measurements, colour suggestions, existing `label_diagnostics`,
#'   `annotation_registry`, `track_registry`, `track_diagnostics`, and captured
#'   device warnings. Field units and provenance are described in
#'   [ggchord_layout_contract]. `measurement$status` concerns panel measurement only;
#'   it never means that all quality or portability checks passed.
#' @section Label quality:
#' `label_quality` reports summary counts, text-object conflicts, registered
#' annotation polylines and measured oriented text boxes in millimetres. Its
#' `measurement` entry separates panel measurement, geometry convergence and
#' physical refinement passes. Circular calibration may try one extrapolated
#' scale, accepting it only after the usual measured convergence gate passes;
#' rejected trials retain the ordinary refinement sequence. Arc
#' glyphs belonging to one string are one text object; their intentional
#' internal contacts are not label collisions. Leader crossings exclude shared
#' endpoints; through-text checks concern box interiors, not exact glyph ink.
#' Source-label visibility follows the layer policy and may share one displayed
#' label across records. It is distinct from the number of text/glyph boxes.
#' Registered gene/feature polygon envelopes are checked against text boxes;
#' `feature_contacts` separates proven own-source placement, other-source
#' contacts and unresolved associations. Transparent polygons remain in this
#' conservative geometry check. Non-polygon strokes/markers and visibility or
#' occlusion are not checked; `feature_parts_unchecked` records skipped parts.
#' Unregistered/outer decorations remain unchecked. These measurements never
#' move, hide or change input records.
#' Degenerate annotation placeholders with coincident endpoints are not
#' counted as leaders. Circular label visibility can additionally report
#' `leader_text_status` for clear, detoured or unresolved static-text routes.
#' `leader_route_status` records bounded crossing detours or unresolved routes
#' while preserving directed endpoints and source-label positions.
#' `axis_clearance` records available automatic circular label decisions and
#' radial shifts in millimetres, positive outward and negative inward.
#' Composite site labels use the same bold name and regular coordinate
#' fragments, physical gap and text styles as the renderer.
#' @md
#' @export
#' @section Guide measurements:
#' `guide_measurements` records active roles, titles, item counts, device text
#' widths/heights, requested key sizes, direction, grid and shared side budget.
#' These complement actual guide-box measurements; requested key dimensions
#' are not claimed as final grob bounds. Font fallback remains uninspected.
#' With [adapt_ggchord()], `adaptive_layout` reports effective theme styles,
#' selected default fixed geom line widths and scale factors. In adaptive mode,
#' `guide_layout` reports the bounded joint
#' search, selected guide grids/positions, explicit locks and measured guide-box
#' bounds in mm. `feasible` concerns the guide canvas and reserved panel budget;
#' it does not certify every label or custom grob. `insufficient_space` retains
#' all categories and reports the failure instead of silently shrinking text.
#'
#' @examples
#' library(dplyr)
#' data(seq_data_example)
#' seqs <- seq_data_example |>
#'   select(accver, length)
#' p <- ggchord(seq_data = seqs) +
#'   geom_seq()
#' diagnose_ggchord(p, width = 8, height = 6, units = "in", device = "png")
diagnose_ggchord <- function(plot, width = 7, height = 5,
    units = c("in", "cm", "mm", "px"), device = c("png", "pdf", "svg"),
    dpi = 144, min_size = 2, min_linewidth = .2, colour_vision = FALSE) {
  caller <- "diagnose_ggchord()"
  if (!inherits(plot, "ggchord")) ggchord_stop(caller, ": plot must be a ggchord plot")
  units <- match.arg(units); device <- match.arg(device)
  for (name in c("width", "height", "dpi")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0)
      ggchord_stop(caller, ": ", name, " must be one positive finite number")
  }
  for (name in c("min_size", "min_linewidth")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0)
      ggchord_stop(caller, ": minimum sizes must be non-negative finite scalars")
  }
  if (!is.logical(colour_vision) || length(colour_vision) != 1L || is.na(colour_vision))
    ggchord_stop(caller, ": colour_vision must be TRUE or FALSE")
  inches <- function(value) switch(units, "in" = value, cm = value / 2.54,
    mm = value / 25.4, px = value / dpi)
  target <- c(width = inches(width), height = inches(height))
  path <- tempfile("ggchord-diagnostic-", fileext = paste0(".", device))
  previous <- grDevices::dev.cur(); owned <- NULL
  on.exit({
    if (!is.null(owned) && owned %in% grDevices::dev.list()) grDevices::dev.off(owned)
    if (previous %in% grDevices::dev.list()) grDevices::dev.set(previous)
    unlink(path)
  }, add = TRUE)
  warnings <- character()
  table <- withCallingHandlers({
    ggchord_diagnostic_device(path, device, target[["width"]], target[["height"]], dpi)
    if (grDevices::dev.cur() == previous) ggchord_stop(caller, ": target device did not open")
    owned <- grDevices::dev.cur()
    built <- ggplot2::ggplot_build(ggchord_render_plot_copy(plot))
    ggchord_render_table(built, force = TRUE)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  context <- attr(table, "ggchord_render_context")
  built <- attr(table, "ggchord_render_built")
  layout <- attr(table, "ggchord_render_layout")
  withCallingHandlers({
    grid::grid.newpage()
    ggchord_with_render_context(context, grid::grid.draw(table))
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  checks <- list(); text <- list()
  add <- function(code, status, detail, layer_id = NA_character_, source_row = NA_integer_) {
    checks[[length(checks) + 1L]] <<- data.frame(code = code, status = status,
      layer_id = layer_id, source_row = as.integer(source_row), detail = detail,
      stringsAsFactors = FALSE)
  }
  measured <- identical(context$status, "measured")
  add("panel_measurement", if (measured) "measured" else "unavailable", context$basis)
  if (measured && !isTRUE(context$converged)) add("physical_layout", "warning",
    "Bounded panel refinement did not converge; physical geometry retains estimated accuracy")
  canvas <- context$canvas_mm
  if (measured) for (i in seq_along(built$data)) {
    data <- built$data[[i]]
    if (!is.data.frame(data) || !nrow(data)) next
    layer <- built$plot$layers[[i]]
    original <- layer$data
    if (!"source_row" %in% names(data) && is.data.frame(original) && nrow(original) == nrow(data)) {
      source_row <- original$source_row
      input <- layer$ggchord_resolved_input
      if (is.null(source_row) && "accver" %in% names(original) &&
          is.data.frame(input) && "accver" %in% names(input))
        source_row <- match(original$accver, input$accver)
      data$source_row <- source_row %||% rep(NA_integer_, nrow(data))
    }
    id <- layer$ggchord_layer_id %||% paste0("ggplot-layer-", i)
    panel <- built$layout$panel_params[[1L]]
    coordinate <- built$plot$coordinates
    if (all(c("x", "y") %in% names(data))) {
      xy <- coordinate$transform(data, panel)
      x <- context$panel_bounds_mm[["xmin"]] + xy$x * context$panel_mm[["width"]]
      y <- context$panel_bounds_mm[["ymin"]] + xy$y * context$panel_mm[["height"]]
      outside <- is.finite(x) & is.finite(y) & (x < 0 | x > canvas[["width"]] |
        y < 0 | y > canvas[["height"]])
      if (any(outside)) add("canvas_geometry_bounds", "warning",
        "Registered geometry vertices extend beyond the canvas", id)
    }
    styled <- ggchord_diagnostic_text_styles(data, layer)
    txt <- ggchord_diagnostic_text(styled, coordinate, panel, context,
      parse = layer$geom_params$parse %||% FALSE, size_unit = layer$geom_params$size.unit %||% 'mm')
    if (nrow(txt)) {
      txt$layer_id <- id; text[[length(text) + 1L]] <- txt
      for (row in which(txt$size_mm < min_size)) add("small_text", "suggestion",
        paste0("Text size is ", signif(txt$size_mm[row], 3), " mm; advisory minimum is ", min_size), id, txt$source_row[row])
      outside <- txt$xmin_mm < 0 | txt$xmax_mm > canvas[["width"]] |
        txt$ymin_mm < 0 | txt$ymax_mm > canvas[["height"]]
      for (row in which(outside)) add("canvas_text_bounds", "warning",
        "Measured text box extends beyond the canvas", id, txt$source_row[row])
    }
    for (column in intersect(c("linewidth", "link_linewidth", "track_linewidth"), names(data))) {
      values <- data[[column]]
      small <- which(is.finite(values) & values > 0 & values < min_linewidth)
      if (length(small)) add("small_linewidth", "suggestion",
        paste0("At least one ", column, " is below ", min_linewidth, " mm"), id)
    }
  }
  guides <- ggchord_diagnostic_guides(table, canvas)
  if (identical(layout$length_mode, "equal")) add("normalized_sequence_lengths", "information",
    "Opt-in equal-length display; source lengths/positions are retained and no shared absolute scale is implied")
  layer_sources <- ggchord_annotation_source_maps(layout)
  annotation_duplicates <- ggchord_annotation_duplicate_sources(layer_sources)
  if (nrow(annotation_duplicates)) add("duplicate_annotation_display", "suggestion",
    "Constructor annotations appear in multiple geometry layers; use disjoint layer data if duplication is unintended")
  if (any(guides$exceeds_canvas %in% TRUE)) add("guide_size", "warning",
    "At least one guide is wider or taller than the target canvas")
  colours <- ggchord_diagnostic_colours(built, colour_vision)
  if (nrow(colours)) add("colour_similarity", "suggestion",
    "Some opaque discrete colours are close under the stated heuristic; inspect colour_suggestions")
  add("font_embedding_and_fallback", "not_checked", "Selected-device text metrics do not establish embedding or identify substituted fonts")
  add("cross_application_transparency", "not_checked", "Transparency in third-party applications was not inspected")
  add("unregistered_grobs_and_arrowheads", "not_checked", "Canvas checks cover registered vertices and text boxes, not arbitrary grobs or arrowhead extents")
  add("colour_vision", if (colour_vision && requireNamespace("colorspace", quietly = TRUE)) "simulated" else "not_checked",
    if (!colour_vision) "Optional simulations were not requested" else
      if (!requireNamespace("colorspace", quietly = TRUE)) "Optional colorspace is unavailable" else
        "colorspace severity=1, linear RGB; Lab-distance hints do not guarantee readability")
  text <- if (length(text)) do.call(rbind, text) else data.frame()
  quality <- ggchord_label_quality(layout, text, built, context)
  if (nrow(quality$conflicts)) add("label_geometry_contacts", "warning",
    "Measured panel text/annotation leaders have contacts; inspect label_quality conflicts and scope")
  structure(list(measurement = context,
    adaptive_layout = context$adaptive_layout, guide_layout = context$guide_layout,
    checks = unique(do.call(rbind, checks)),
    text = text, guides = guides, label_quality = quality,
    guide_measurements=ggchord_measure_guides(built$plot,guides),
    colour_suggestions = colours, label_diagnostics = layout$collinear_label_diagnostics %||% data.frame(),
    label_visibility = ggchord_label_visibility(layout),
    annotation_registry = layout$circular_annotation_registry %||% data.frame(),
    layer_sources = layer_sources, annotation_duplicates = annotation_duplicates,
    track_registry = layout$track_registry %||% data.frame(),
    track_diagnostics = layout$track_diagnostics %||% data.frame(),
    warnings = unique(warnings)), class = "ggchord_diagnostics")
}

#' @export
print.ggchord_diagnostics <- function(x, ...) {
  cat("ggchord target-size diagnostics\n")
  cat("Panel measurement: ", x$measurement$status, "\n", sep = "")
  cat("Warnings: ", sum(x$checks$status == "warning"),
    "; suggestions: ", sum(x$checks$status == "suggestion"),
    "; unchecked areas: ", sum(x$checks$status == "not_checked"), "\n", sep = "")
  invisible(x)
}
