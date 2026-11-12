# Isolated measured-refinement/target-metrics candidate on canonical 8x8.
helper <- readLines('tools/long-term-audit/check-label-quality.R', warn = FALSE)
end <- grep('^# Explicit geometry controls', helper)[1L]-1L
eval(parse(text = helper[seq_len(end)]))
lines <- readLines('tools/documentation-examples.R', warn = FALSE)
marks <- grep('^# ---- [a-z]+ ----$', lines)
chunks <- setNames(lapply(seq_len(length(marks)-1L), function(i)
  lines[seq.int(marks[i]+1L, marks[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$', '\\1', lines[marks[-length(marks)]]))
for (name in c('docstyle', 'showcaseplasmiddata', 'showcaseplasmid')) eval(parse(text = chunks[[name]]))
namespace <- asNamespace('ggchord')
prototype_refinement_log <- new.env(parent = emptyenv())
prototype_refinement_log$events <- list()
candidate_render <- function(built, force = FALSE) {
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
    passes_allowed <- if (isTRUE(source$coordinates$ggchord_circular)) 10L else 3L
    for (pass in seq_len(passes_allowed)) {
      # Try one scalar extrapolation only for an observed circular contraction.
      if (isTRUE(source$coordinates$ggchord_circular) && !acceleration$attempted &&
          length(history) >= 2L && abs(diff(context$units_per_mm))/max(context$units_per_mm) < 1e-6) {
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
          prototype_refinement_log$events <- c(prototype_refinement_log$events,
            list(list(input = proposed$units_per_mm, slope = slope, relative = trial$relative,
              accepted = accepted, message = trial$message)))
          cat('Verified physical candidate:', predicted, 'residual', trial$relative,
            'accepted', accepted, '\n'); flush.console()
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

private <- new.env(parent = namespace)
private$prototype_refinement_log <- prototype_refinement_log
environment(candidate_render) <- private
unlockBinding('ggchord_render_table', namespace)
assign('ggchord_render_table', candidate_render, envir = namespace)
lockBinding('ggchord_render_table', namespace)
text_override <- '--target' %in% commandArgs(trailingOnly = TRUE)
if (text_override) {
  text_candidate <- get('ggchord_text_boxes', envir = namespace)
  rewrite <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (identical(expr, quote(ggchord_measurement_device(text = texts[valid]))))
      return(quote(if (grDevices::dev.cur() > 1L) function() invisible(NULL) else
        ggchord_measurement_device(text = texts[valid])))
    as.call(lapply(as.list(expr), rewrite))
  }
  body(text_candidate) <- rewrite(body(text_candidate))
  unlockBinding('ggchord_text_boxes', namespace)
  assign('ggchord_text_boxes', text_candidate, envir = namespace)
  lockBinding('ggchord_text_boxes', namespace)
}

prefix <- 'v0.27.0-flagship-refinement-prototype-'
result <- measure(p_show_plasmid, 'square', 8, 8)
sources <- list.files('R', '[.]R$', full.names = TRUE)
jsonlite::write_json(list(scope = 'In-process candidate on canonical 8x8; not native production certification',
  production_overrides = list(physical_refinement = TRUE, text_metrics = text_override),
  result = result, events = prototype_refinement_log$events,
  baseline_runtime_sha256 = as.list(setNames(vapply(sources,
    function(file) digest::digest(file = file, algo = 'sha256'), character(1)), sources))),
  'tools/long-term-audit/v0.27.0-flagship-refinement-prototype-results.json',
  auto_unbox = TRUE, pretty = TRUE, na = 'null')
cat('Canonical 8x8 measured-refinement candidate recorded\n')
