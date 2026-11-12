# A render-local measurement context. It is restored after nested builds and
# never stored in global R options or reused by another plot/device.
.ggchord_render_state <- new.env(parent = emptyenv())
.ggchord_render_state$context <- NULL

ggchord_render_context <- function() .ggchord_render_state$context

ggchord_render_cache_key <- function() {
  context <- ggchord_render_context()
  list(status = context$status %||% "estimated", units_per_mm = context$units_per_mm,
    canvas_mm = ggchord_render_dimensions(), device = grDevices::dev.cur())
}

ggchord_with_render_context <- function(context, expr) {
  previous <- .ggchord_render_state$context
  .ggchord_render_state$context <- context
  on.exit(.ggchord_render_state$context <- previous, add = TRUE)
  eval(substitute(expr), envir = parent.frame())
}

ggchord_render_plot_copy <- function(plot) {
  plot <- ggchord_register_layers(plot)
  plot$ggchord$ref <- new.env(parent = emptyenv())
  plot$ggchord$layout <- NULL
  plot$scales <- plot$scales$clone()
  plot
}

ggchord_render_dimensions <- function() {
  available <- ggchord_render_context()$available_canvas_mm
  if (length(available) == 2L && all(is.finite(available))) return(available)
  if (grDevices::dev.cur() == 1L) return(c(width = NA_real_, height = NA_real_))
  value <- tryCatch(grDevices::dev.size("in") * 25.4,
    error = function(e) c(NA_real_, NA_real_))
  stats::setNames(value, c("width", "height"))
}

ggchord_estimated_render_context <- function(layout = NULL) {
  unit <- layout$text_units_per_inch %||% NA_real_
  list(status = if (is.finite(unit)) "estimated" else "unavailable",
    panel_mm = c(width = NA_real_, height = NA_real_),
    canvas_mm = ggchord_render_dimensions(),
    units_per_inch = unit, axis_units_per_inch = c(x = unit, y = unit),
    units_per_mm = c(x = unit, y = unit) / 25.4,
    geometry_metric_status = "estimated", converged = NA,
    device = if (grDevices::dev.cur() == 1L) "none" else names(grDevices::dev.cur()),
    font_conditions = "Font fallback and embedding have not been inspected",
    basis = "Layout-only estimate; actual panel has not been measured")
}

ggchord_panel_context <- function(table, built) {
  unavailable <- ggchord_estimated_render_context(built$plot$ggchord$ref$layout)
  if (grDevices::dev.cur() == 1L) return(unavailable)
  panel <- which(grepl("^panel($|[-.])", table$layout$name))
  if (length(panel) != 1L) {
    unavailable$basis <- "Actual-panel measurement currently supports one panel"
    return(unavailable)
  }
  entry <- table$layout[panel, , drop = FALSE]
  range <- built$layout$panel_params[[1L]]
  xr <- range$x.range %||% range$x$continuous_range
  yr <- range$y.range %||% range$y$continuous_range
  if (length(xr) != 2L || length(yr) != 2L || any(!is.finite(c(xr, yr)))) return(unavailable)
  result <- tryCatch({
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      length(table$heights), length(table$widths), widths = table$widths,
      heights = table$heights, respect = table$respect)))
    on.exit(grid::popViewport(), add = TRUE)
    grid::pushViewport(grid::viewport(layout.pos.row = seq.int(entry$t, entry$b),
      layout.pos.col = seq.int(entry$l, entry$r)))
    on.exit(grid::popViewport(), add = TRUE)
    size <- c(width = grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE),
      height = grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE))
    lower <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
    upper <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
    list(size = size, bounds = c(xmin = lower$x, xmax = upper$x,
      ymin = lower$y, ymax = upper$y) * 25.4)
  }, error = function(e) NULL)
  if (is.null(result) || any(!is.finite(result$size)) || any(result$size <= 0)) return(unavailable)
  units <- c(x = abs(diff(xr)) / result$size[["width"]],
    y = abs(diff(yr)) / result$size[["height"]])
  list(status = "measured", panel_mm = result$size,
    canvas_mm = ggchord_render_dimensions(), panel_bounds_mm = result$bounds,
    available_canvas_mm = ggchord_render_context()$available_canvas_mm,
    x_range = xr, y_range = yr, units_per_mm = units,
    axis_units_per_inch = units * 25.4, units_per_inch = max(units) * 25.4,
    geometry_metric_status = "estimated", converged = FALSE,
    device = names(grDevices::dev.cur()),
    font_conditions = "Metrics use the selected device; font fallback and embedding remain unchecked",
    basis = "Resolved gtable grid.layout panel, including null units and aspect constraints")
}

ggchord_physical_requested <- function(plot) {
  any(vapply(plot$layers, function(layer) {
    params <- layer$ggchord_params
    if (is.null(params)) return(FALSE)
    params$type %in% c("seq", "track") || isTRUE(params$physical_layout) ||
      identical(params$gene_label_fit, "shrink") || identical(params$feature_label_fit, "shrink")
  }, logical(1)))
}

ggchord_base_gtable <- function(built) {
  class(built) <- setdiff(class(built), "ggchord_built")
  ggplot2::ggplot_gtable(built)
}

ggchord_render_table <- function(built, force = FALSE) {
  source <- attr(built, "ggchord_source_plot")
  if (is.null(source)) source <- built$plot
  guided <- function(value) {
    if (identical(source$ggchord$layout_policy$mode,"adaptive"))
      ggchord_guide_joint_layout(value,source) else list(built=value,table=ggchord_base_gtable(value))
  }
  result <- guided(built)
  built <- result$built; table <- result$table %||% ggchord_base_gtable(built)
  context <- ggchord_panel_context(table, built)
  physical <- isTRUE(force) || !is.null(source$ggchord$layout_policy) || ggchord_physical_requested(source)
  passes <- 0L
  evaluations <- 0L
  history <- list()
  acceleration <- list(attempted = FALSE, accepted = FALSE)
  if (physical && identical(context$status, "measured")) {
    passes_allowed <- 10L
    for (pass in seq_len(passes_allowed)) {
      # Try one scalar extrapolation only for an observed circular contraction.
      if (isTRUE(source$coordinates$ggchord_circular) && !acceleration$attempted &&
          length(history) >= 2L && length(context$units_per_mm) == 2L &&
          all(is.finite(context$units_per_mm)) && all(context$units_per_mm > 0) &&
          abs(diff(context$units_per_mm))/max(context$units_per_mm) < 1e-6) {
        h0 <- history[[length(history)-1L]]; h1 <- history[[length(history)]]
        delta <- h1[['input']]-h0[['input']]
        slope <- if (abs(delta) > 1e-12) (h1[['output']]-h0[['output']])/delta else NA_real_
        predicted <- if (is.finite(slope) && slope >= 0 && slope < .95)
          (h1[['output']]-slope*h1[['input']])/(1-slope) else NA_real_
        current <- max(context$units_per_mm)
        if (is.finite(predicted) && predicted > 0 && predicted >= current*.5 && predicted <= current*3) {
          acceleration$attempted <- TRUE
          proposed <- context
          proposed$units_per_mm <- stats::setNames(rep(predicted, 2L), c('x', 'y'))
          proposed$axis_units_per_inch <- proposed$units_per_mm*25.4
          proposed$units_per_inch <- predicted*25.4
          proposed$geometry_metric_status <- 'estimated'
          proposed$basis <- 'Extrapolated scale candidate; actual panel must be rebuilt and measured'
          evaluations <- evaluations+1L
          trial <- tryCatch({
            trial_plot <- ggchord_render_plot_copy(source)
            trial_built <- ggchord_with_render_context(proposed, ggplot2::ggplot_build(trial_plot))
            trial_result <- ggchord_with_render_context(proposed, guided(trial_built))
            trial_built <- trial_result$built
            trial_table <- trial_result$table %||% ggchord_base_gtable(trial_built)
            trial_context <- ggchord_panel_context(trial_table, trial_built)
            error <- if (identical(trial_context$status, 'measured'))
              max(abs(trial_context$units_per_mm/proposed$units_per_mm-1)) else Inf
            list(built = trial_built, result = trial_result, table = trial_table,
              context = trial_context, relative = error)
          }, error = function(e) list(relative = Inf, message = conditionMessage(e)))
          accepted <- is.finite(trial$relative) && trial$relative <= .005 &&
            !identical(trial$result$report$status, 'insufficient_space')
          acceleration$accepted <- accepted
          acceleration$proposed_units_per_mm <- predicted
          acceleration$measured_relative <- trial$relative
          if (!is.null(trial$message)) acceleration$failure_message <- trial$message
          if (accepted) {
            built <- trial$built; result <- trial$result; table <- trial$table
            context <- trial$context
            context$converged <- TRUE; context$geometry_metric_status <- 'measured'
            passes <- evaluations
            break
          }
        }
      }
      input <- context
      next_plot <- ggchord_render_plot_copy(source)
      built <- ggchord_with_render_context(input, ggplot2::ggplot_build(next_plot))
      result <- ggchord_with_render_context(input, guided(built))
      built <- result$built; table <- result$table %||% ggchord_base_gtable(built)
      context <- ggchord_panel_context(table, built)
      evaluations <- evaluations+1L
      passes <- evaluations
      if (!identical(context$status, "measured")) break
      relative <- max(abs(context$units_per_mm / input$units_per_mm - 1))
      context$converged <- is.finite(relative) && relative <= .005
      context$geometry_metric_status <- if (isTRUE(context$converged)) "measured" else "estimated"
      history <- c(history, list(c(input = max(input$units_per_mm), output = max(context$units_per_mm))))
      if (isTRUE(context$converged)) break
    }
  }
  context$refinement_passes <- passes
  context$refinement_acceleration <- acceleration
  if (physical && identical(context$status, 'measured') && !isTRUE(context$converged)) warning(
    'ggchord(): physical geometry did not converge within the bounded refinement budget; geometry measurements remain estimates. Inspect diagnose_ggchord() and consider a larger canvas or revised label/guide settings.',
    call. = FALSE)
  context$adaptive_layout <- built$plot$ggchord$adaptive_report
  context$guide_layout <- result$report
  if (identical(result$report$status,"insufficient_space")) warning(
    "adapt_ggchord(): guides cannot fit the canvas and panel budget under the supplied constraints; increase width/height or revise explicit guide layout. See guide_layout diagnostics.",call.=FALSE)
  layout <- built$plot$ggchord$ref$layout
  if (!is.null(layout)) {
    layout$render_context <- context
    layout$adaptive_layout <- context$adaptive_layout
    layout$guide_layout <- context$guide_layout
    if (is.data.frame(layout$track_registry) && nrow(layout$track_registry))
      layout$track_registry$measurement_status <- context$geometry_metric_status
  }
  attr(table, "ggchord_render_context") <- context
  attr(table, "ggchord_render_layout") <- layout
  attr(table, "ggchord_render_built") <- built
  table
}


#' @export
#' @importFrom ggplot2 ggplot_gtable
ggplot_gtable.ggchord_built <- function(data) ggchord_render_table(data)
