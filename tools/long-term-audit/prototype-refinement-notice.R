# Candidate failure explanation over real completed measurement contexts.
# Does not render, change source labels or modify the installed namespace.
source('tools/audit-source-snapshot.R')
snapshot <- audit_source_snapshot('tools/long-term-audit/prototype-refinement-notice.R')
notice <- function(context, physical) {
  context$refinement_reason <- if (!isTRUE(physical)) 'not_requested' else
    if (!identical(context$status, 'measured')) 'panel_not_measured' else
      if (isTRUE(context$converged)) 'tolerance_met' else 'bounded_refinement_exhausted'
  if (identical(context$refinement_reason, 'bounded_refinement_exhausted')) warning(
    'ggchord(): physical geometry did not converge within the bounded refinement budget; geometry measurements remain estimates. Inspect diagnose_ggchord() and consider a larger canvas or revised label/guide settings.',
    call. = FALSE)
  context
}
folder <- 'tools/long-term-audit'
inputs <- c('v0.27.0-flagship-device-published-matched-square-results.json',
  'v0.27.0-flagship-device-compact-results.json',
  'v0.27.0-label-transition-device-58-results.json',
  'v0.27.0-label-transition-device-10x10-58-results.json')
observations <- list()
for (input in inputs) {
  report <- jsonlite::fromJSON(file.path(folder, input), simplifyVector = FALSE)
  cases <- if (!is.null(report$profiles)) report$profiles else report$cases
  for (id in names(cases)) {
    record <- cases[[id]]; messages <- character()
    context <- list(status = record$measurement, converged = record$geometry_converged,
      geometry_metric_status = record$geometry_metric_status)
    outcome <- withCallingHandlers(notice(context, TRUE), warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart('muffleWarning')
    })
    stopifnot(length(messages) == if (isTRUE(context$converged)) 0L else 1L)
    observations[[paste(input,id,sep='/')]] <- list(input_sha256 = digest::digest(
      file = file.path(folder,input),algo='sha256'), outcome = outcome, warnings = messages)
  }
}
audit_assert_source_snapshot(snapshot)
jsonlite::write_json(list(scope = 'Candidate notice over six recorded real contexts; native render integration remains unverified',
  production_override = FALSE, production_implementation = FALSE,
  observations = observations, baseline_runtime_sha256 = as.list(snapshot$engine),
  candidate_body_sha256 = digest::digest(deparse(body(notice)),algo='sha256')),
  file.path(folder,'v0.27.0-refinement-notice-prototype-results.json'),
  pretty=TRUE,auto_unbox=TRUE,na='null')
cat('Candidate explanation observations recorded; native integration not verified\n')
