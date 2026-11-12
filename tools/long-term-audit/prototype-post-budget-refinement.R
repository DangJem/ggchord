# One bounded axis prediction only after ordinary non-circular refinement fails.
helper <- readLines('tools/long-term-audit/check-label-quality.R',warn=FALSE)
cut <- grep('^# Explicit geometry controls',helper)[1L]-1L
helper <- helper[seq_len(cut)]
stopifnot(sum(grepl('width, height, 110',helper,fixed=TRUE))==1L)
helper <- sub('width, height, 110','width, height, 144',helper,fixed=TRUE)
eval(parse(text=helper))
audit_session <- audit_source_snapshot(c(audit_session$additional,
  'tools/long-term-audit/prototype-post-budget-refinement.R',
  'tools/long-term-audit/linear-lesson-source.R'))
source('tools/long-term-audit/linear-lesson-source.R')
lesson <- load_linear_lesson()
original_features <- lesson$show_features
namespace <- asNamespace('ggchord')
candidate <- get('ggchord_render_table',namespace)
replacements <- 0L
rewrite <- function(expr) {
  if(!is.call(expr)) return(expr)
  if(identical(expr,quote(history <- c(history,list(c(input=max(input$units_per_mm),output=max(context$units_per_mm))))))) {
    replacements <<- replacements+1L
    return(quote(history <- c(history,list(c(input=max(input$units_per_mm),output=max(context$units_per_mm),
      input_x=input$units_per_mm[[1L]],input_y=input$units_per_mm[[2L]],
      output_x=context$units_per_mm[[1L]],output_y=context$units_per_mm[[2L]])))))
  }
  if(identical(expr,quote(context$refinement_passes <- passes))) {
    replacements <<- replacements+1L
    return(quote({
      if(physical && !isTRUE(source$coordinates$ggchord_circular) &&
          identical(context$status,'measured') && !isTRUE(context$converged) &&
          length(history)>=2L && length(context$units_per_mm)==2L &&
          all(is.finite(context$units_per_mm)) && all(context$units_per_mm>0)) {
        h0 <- history[[length(history)-1L]];h1 <- history[[length(history)]]
        a <- unname(h0[c('input_x','input_y')]);b <- unname(h1[c('input_x','input_y')])
        out0 <- unname(h0[c('output_x','output_y')]);out1 <- unname(h1[c('output_x','output_y')])
        delta <- b-a
        slope <- (out1-out0)/delta
        predicted <- (out1-slope*b)/(1-slope)
        stable <- abs(delta)<=1e-12*pmax(1,abs(a),abs(b)) &
          abs(out1-b)/pmax(abs(b),1e-12)<=.005
        predicted[stable] <- context$units_per_mm[stable]
        allowed <- stable | (is.finite(slope) & slope>=0 & slope<.95)
        if(all(allowed) && all(is.finite(predicted)) && all(predicted>0) &&
            all(predicted>=context$units_per_mm*.5) && all(predicted<=context$units_per_mm*3)) {
          proposed <- context
          proposed$units_per_mm <- stats::setNames(predicted,c('x','y'))
          proposed$axis_units_per_inch <- proposed$units_per_mm*25.4
          proposed$units_per_inch <- max(predicted)*25.4
          proposed$geometry_metric_status <- 'estimated'
          proposed$basis <- 'Post-budget axis candidate; actual panel must be rebuilt and measured'
          trial <- tryCatch({
            trial_plot <- ggchord_render_plot_copy(source)
            trial_built <- ggchord_with_render_context(proposed,ggplot2::ggplot_build(trial_plot))
            trial_result <- ggchord_with_render_context(proposed,guided(trial_built))
            trial_built <- trial_result$built
            trial_table <- trial_result$table %||% ggchord_base_gtable(trial_built)
            trial_context <- ggchord_panel_context(trial_table,trial_built)
            relative <- if(identical(trial_context$status,'measured'))
              max(abs(trial_context$units_per_mm/proposed$units_per_mm-1)) else Inf
            list(built=trial_built,result=trial_result,table=trial_table,context=trial_context,relative=relative)
          },error=function(e) list(relative=Inf,message=conditionMessage(e)))
          accepted <- is.finite(trial$relative) && trial$relative<=.005 &&
            !identical(trial$result$report$status,'insufficient_space')
          passes <- passes+1L
          acceleration <- list(attempted=TRUE,accepted=accepted,phase='after_ordinary_budget',
            proposed_units_per_mm=max(predicted),proposed_axis_units_per_mm=proposed$units_per_mm,
            measured_relative=trial$relative)
          if(!is.null(trial$message)) acceleration$failure_message <- trial$message
          if(accepted) {
            built <- trial$built;result <- trial$result;table <- trial$table;context <- trial$context
            context$converged <- TRUE;context$geometry_metric_status <- 'measured'
          }
        }
      }
      context$refinement_passes <- passes
    }))
  }
  as.call(lapply(as.list(expr),rewrite))
}
body(candidate) <- rewrite(body(candidate));stopifnot(replacements==2L)
unlockBinding('ggchord_render_table',namespace)
assign('ggchord_render_table',candidate,namespace)
lockBinding('ggchord_render_table',namespace)
prefix <- 'v0.27.0-linear-budget-prototype-'
warnings <- character()
result <- withCallingHandlers(measure(lesson$p_linear_labels,'linearlabels',10,6),warning=function(w) {
  warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
})
stopifnot(identical(original_features,lesson$show_features))
audit_assert_source_snapshot(audit_session)
jsonlite::write_json(list(scope='In-process post-budget non-circular axis prediction on the unchanged 10x6 inch/144 dpi linear-label lesson; not native certification',
  production_override=TRUE,measurement_helper_dpi=144,
  candidate_body_sha256=digest::digest(deparse(body(candidate)),algo='sha256'),
  result=result,warnings=warnings,source_features_preserved=TRUE,
  baseline_runtime_sha256=as.list(audit_session$engine)),
  'tools/long-term-audit/v0.27.0-linear-budget-prototype-results.json',pretty=TRUE,auto_unbox=TRUE,na='null')
cat('Post-budget axis candidate observations recorded\n')
