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

#' Diagnose plot references, layout and target-device fit
#' @inheritParams export_ggchord_layout
#' @param action Return a report, warn on issues, or stop on error-level issues.
#' @return A structured issue table and measurement conditions.
#' @export
#' @examples
#' seqs <- data.frame(
#'   accver = c("A", "B"), length = c(1000, 800),
#'   sequence_unit = "bp", topology = "linear"
#' )
#'
#' genes <- data.frame(
#'   accver = c("A", "A", "B"), feature_id = c("g1", "g2", "g3"),
#'   start = c(100, 400, 150), end = c(240, 620, 380),
#'   strand = c("+", "-", "+"), label = c("g1", "g2", "g3"),
#'   type = c("CDS", "tRNA", "CDS")
#' )
#'
#' p_basic <- ggchord(seqs) + geom_seq() +
#'   geom_gene(data = genes) + coord_collinear()
#'
#' diagnose_ggchord(p_basic, width = 8, height = 4, device = "pdf")
diagnose_ggchord <- function(plot, width = NULL, height = NULL, units = "in",
    device = NULL, dpi = 300, action = "report") {
  caller <- "diagnose_ggchord()"
  ggchord_output_plot(plot, caller)
  action <- match.arg(action, c("report", "warn", "error"))
  request <- ggchord_output_request(width, height, units, device, dpi, caller)
  warnings <- character()
  failure <- NULL
  built <- tryCatch(withCallingHandlers({
    if (is.null(request$size)) {
      ggchord_diagnose_sources(plot)
    } else {
      ggchord_with_output_device(request, function() {
        result <- ggchord_output_build_geometry(ggchord_query_clone(plot))
        ggchord_diagnose_built(result, request$size * 25.4)
      })
    }
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    failure <<- conditionMessage(e)
    NULL
  })
  measurement <- ggchord_output_measurement(request, !is.null(request$size), failure,
    actual = attr(built, "measurement_device"))
  severe <- if (!is.null(failure)) ggchord_output_issue("error", "build_failed", failure)
  warning_rows <- lapply(unique(warnings), function(message)
    ggchord_output_issue("warning", "plot_warning", message))
  problems <- ggchord_rbind_fill(c(list(built, severe), warning_rows)) %||%
    ggchord_output_issue("information", "", "")[FALSE, ]
  errors <- problems$message[problems$severity == "error"]
  warned <- problems$message[problems$severity %in% c("warning", "error")]
  if (action == "error" && length(errors)) ggchord_stop(caller, ": ", paste(errors, collapse = "\n"))
  if (action == "warn" && length(warned)) warning(paste(warned, collapse = "\n"), call. = FALSE)
  structure(list(issues = problems, measurement = measurement), class = "ggchord_diagnostics")
}

ggchord_diagnose_sources <- function(plot) {
  issues <- list()
  check <- function(expr, path, layer = NA_character_) {
    tryCatch(expr, error = function(e) {
      issues[[length(issues) + 1L]] <<- ggchord_output_issue(
        "error", "source_failed", conditionMessage(e), path, layer, "source")
      NULL
    })
  }
  walk <- function(plot, path) {
    if (inherits(plot, "ggchord_fusion")) {
      for (layer in plot$ribbons) check(ggchord_fusion_input(layer), path, layer$id)
      for (id in names(plot$panels))
        walk(plot$panels[[id]]$plot, c(path, ggchord_fusion_path_ids(id)))
      return(invisible(NULL))
    }
    if (!inherits(plot, "ggchord")) {
      ids <- ggchord_output_layer_ids(plot)
      for (i in seq_along(plot$layers))
        check(plot$layers[[i]]$layer_data(plot$data), path, ids[i])
      return(invisible(NULL))
    }
    check({
      plot <- ggchord_render_inputs(plot)
      plot <- ggchord_register_layers(plot)
      ggchord_scale_scopes(plot)
      if (isTRUE(plot$coordinates$ggchord_dotplot)) {
        ggchord_dotplot_geometry(plot)
      } else {
        plot <- ggchord_prepare_geometry_inputs(plot)
        for (layer in plot$layers) {
          if (is.null(layer$ggchord_type)) {
            check(layer$layer_data(plot$data), path, layer$ggchord_layer_id %||% NA_character_)
            next
          }
          check({
            input <- ggchord_resolve_layer_input(layer,
              ggchord_layer_fallback(layer, plot$ggchord$data))
            if (identical(layer$ggchord_params$type, "gene")) {
              ggchord_require_columns(input, c("accver", "start", "end"), "Feature layer")
              ggchord_check_reference_ids(input, "accver",
                plot$ggchord$data$seq_data$accver, "Feature layer")
            }
          }, path, layer$ggchord_layer_id)
        }
      }
    }, path)
    invisible(NULL)
  }
  walk(plot, character())
  ggchord_rbind_fill(issues) %||% ggchord_output_issue("information", "", "")[FALSE, ]
}

ggchord_diagnose_built <- function(result, canvas_mm) {
  issues <- list()
  add <- function(severity, code, message, path = character(), layer = NA_character_,
      component = "panel", suggestion = NA_character_) {
    issues[[length(issues) + 1L]] <<- ggchord_output_issue(severity, code, message, path,
      layer, component, suggestion)
  }
  for (i in seq_along(result$leaves)) {
    key <- names(result$leaves)[i]
    path <- ggchord_fusion_path_ids(key)
    leaf <- result$leaves[[i]]
    context <- leaf$context
    if (!identical(context$status, "measured")) {
      add("warning", "panel_measurement_unavailable", context$basis %||% "Unable to measure panel", path)
      next
    }
    if (identical(context$converged, FALSE)) add("warning", "geometry_not_converged",
      "The physical geometry did not converge", path, suggestion = "Increase available panel size")
    built <- leaf$built
    for (j in seq_along(built$data)) {
      data <- built$data[[j]]
      layer <- built$plot$layers[[j]]
      id <- layer$ggchord_layer_id %||% paste0("ggplot-layer-", j)
      if (!is.data.frame(data) || !nrow(data)) next
      styled <- ggchord_diagnostic_text_styles(data, layer)
      groups <- split(seq_len(nrow(styled)), as.integer(styled$PANEL))
      for (rows in groups) {
        panel <- as.integer(styled$PANEL[rows[1L]])
        text <- ggchord_diagnostic_text(styled[rows, , drop = FALSE], built$layout$coord,
          built$layout$panel_params[[panel]], context$facet_contexts[[panel]] %||% context,
          parse = layer$geom_params$parse %||% FALSE, size_unit = layer$geom_params$size.unit %||% "mm")
        if (nrow(text) && any(text$xmin_mm < 0 | text$ymin_mm < 0 |
            text$xmax_mm > canvas_mm[1L] | text$ymax_mm > canvas_mm[2L]))
          add("warning", "canvas_text_bounds", "Text extends beyond the target canvas", path, id,
            "text", "Increase output dimensions or revise label placement")
      }
    }
    guides <- leaf$guides
    if (any(guides$exceeds_canvas %in% TRUE)) add("warning", "guide_size",
      "A guide exceeds the target canvas", path, component = "guide",
      suggestion = "Change the guide layout or increase the output dimensions")
  }
  ggchord_rbind_fill(issues) %||% ggchord_output_issue("information", "", "")[FALSE, ]
}

#' @export
print.ggchord_diagnostics <- function(x, ...) {
  cat("ggchord diagnostics\n")
  cat("Measurement: ", x$measurement$status, "\n", sep = "")
  print(x$issues, row.names = FALSE)
  invisible(x)
}
