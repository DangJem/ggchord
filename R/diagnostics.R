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

ggchord_diagnostic_text <- function(data, coordinate, panel, context, parse = FALSE) {
  if (!all(c("x", "y", "label", "size") %in% names(data))) return(data.frame())
  selected <- !is.na(data$label) & nzchar(as.character(data$label)) &
    is.finite(data$size) & data$size > 0
  if (".component" %in% names(data) && any(!is.na(data$.component))) selected <- selected &
    data$.component %in% c("text", "label", "arc_text", "name", "length", "track_axis_text")
  if (all(c(".draw_as_arc", ".component") %in% names(data))) selected <- selected &
    !(data$.draw_as_arc %in% TRUE & data$.component %in% "text")
  data <- data[selected, , drop = FALSE]
  if (!nrow(data)) return(data.frame())
  xy <- coordinate$transform(data, panel)
  result <- lapply(seq_len(nrow(data)), function(i) {
    value <- function(name, fallback) if (is.null(data[[name]])) fallback else data[[name]][i]
    label <- if (isTRUE(parse) || isTRUE(value("plotmath_label", FALSE))) {
      tryCatch(parse(text = as.character(data$label[i]))[[1L]], error = function(e) as.character(data$label[i]))
    } else data$label[i]
    grob <- grid::textGrob(label, gp = grid::gpar(
      fontsize = data$size[i] * (72.27 / 25.4), fontfamily = value("family", ""),
      fontface = value("fontface", 1), lineheight = value("lineheight", 1.2)))
    w <- grid::convertWidth(grid::grobWidth(grob), "mm", valueOnly = TRUE)
    h <- grid::convertHeight(grid::grobHeight(grob), "mm", valueOnly = TRUE)
    angle <- value("angle", 0) * pi / 180
    justify <- function(z, at) {
      if (is.numeric(z) && is.finite(z)) return(z)
      known <- c(left = 0, bottom = 0, centre = .5, center = .5, right = 1, top = 1)
      if (as.character(z) %in% names(known)) return(unname(known[as.character(z)]))
      if (identical(as.character(z), "inward")) return(if (at < .5) 0 else if (at > .5) 1 else .5)
      if (identical(as.character(z), "outward")) return(if (at < .5) 1 else if (at > .5) 0 else .5)
      .5
    }
    hj <- justify(value("hjust", .5), xy$x[i]); vj <- justify(value("vjust", .5), xy$y[i])
    cx <- xy$x[i] * context$panel_mm[["width"]] + context$panel_bounds_mm[["xmin"]] +
      (.5 - hj) * w * cos(angle) - (.5 - vj) * h * sin(angle)
    cy <- xy$y[i] * context$panel_mm[["height"]] + context$panel_bounds_mm[["ymin"]] +
      (.5 - hj) * w * sin(angle) + (.5 - vj) * h * cos(angle)
    bw <- abs(cos(angle)) * w + abs(sin(angle)) * h
    bh <- abs(sin(angle)) * w + abs(cos(angle)) * h
    data.frame(source_row = value("source_row", NA_integer_), label = as.character(data$label[i]),
      size_mm = as.numeric(data$size[i]), xmin_mm = cx - bw / 2, xmax_mm = cx + bw / 2,
      ymin_mm = cy - bh / 2, ymax_mm = cy + bh / 2, stringsAsFactors = FALSE)
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
#'   device warnings. `measurement$status` concerns panel measurement only;
#'   it never means that all quality or portability checks passed.
#' @md
#' @export
#' @examples
#' p <- ggchord(data.frame(accver = "A", length = 1000)) +
#'   geom_seq() + coord_collinear()
#' report <- diagnose_ggchord(p, width = 7, height = 4, device = "pdf")
#' report$measurement$status
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
    txt <- ggchord_diagnostic_text(data, coordinate, panel, context, parse = layer$geom_params$parse %||% FALSE)
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
  structure(list(measurement = context,
    checks = unique(do.call(rbind, checks)),
    text = if (length(text)) do.call(rbind, text) else data.frame(), guides = guides,
    colour_suggestions = colours, label_diagnostics = layout$collinear_label_diagnostics %||% data.frame(),
    annotation_registry = layout$circular_annotation_registry %||% data.frame(),
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
