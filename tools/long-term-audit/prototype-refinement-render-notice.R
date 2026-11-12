# Actual PNG render using only an in-process final notice insertion.
script <- readLines('tools/long-term-audit/check-flagship-label-quality.R', warn = FALSE)
cut <- grep('^results <- list', script)[1L]
eval(parse(text = script[seq_len(cut)]))
namespace <- asNamespace('ggchord')
candidate <- get('ggchord_render_table', namespace)
replacements <- 0L
rewrite <- function(expr) {
  if (!is.call(expr)) return(expr)
  if (identical(expr, quote(context$refinement_acceleration <- acceleration))) {
    replacements <<- replacements+1L
    return(quote({
      context$refinement_acceleration <- acceleration
      if (physical && identical(context$status, 'measured') && !isTRUE(context$converged)) warning(
        'ggchord(): physical geometry did not converge within the bounded refinement budget; geometry measurements remain estimates. Inspect diagnose_ggchord() and consider a larger canvas or revised label/guide settings.',
        call. = FALSE)
    }))
  }
  as.call(lapply(as.list(expr), rewrite))
}
body(candidate) <- rewrite(body(candidate)); stopifnot(replacements == 1L)
unlockBinding('ggchord_render_table', namespace)
assign('ggchord_render_table', candidate, namespace)
lockBinding('ggchord_render_table', namespace)
prefix <- 'v0.27.0-refinement-notice-render-prototype-'
notices <- character()
withCallingHandlers(eval(parse(text = script[-seq_len(cut)])), warning = function(w) {
  notices <<- c(notices, conditionMessage(w)); invokeRestart('muffleWarning')
})
file <- paste0('tools/long-term-audit/',prefix,paste(names(profiles),collapse='-'),'-results.json')
report <- jsonlite::fromJSON(file, simplifyVector=FALSE)
stopifnot(length(grep('physical geometry did not converge',notices,fixed=TRUE)) == 1L,
  identical(report$profiles$compact$geometry_converged,FALSE),
  report$profiles$compact$n_visible_source_labels == 164L)
report$baseline_runtime_sha256 <- report$runtime_sha256; report$runtime_sha256 <- NULL
report$production_override <- 'Notice-only insertion, actual compact PNG render; no native implementation certification'
report$warnings <- notices
report$candidate_body_sha256 <- digest::digest(deparse(body(candidate)),algo='sha256')
# One real positive build verifies that ordinary convergence remains quiet.
data(plasmid_example_pUC19)
positive <- ggchord(plasmid_example_pUC19[c('accver','length')]) + geom_seq() + coord_circular()
positive_warnings <- character()
positive_path <- tempfile(fileext='.png')
ggchord_diagnostic_device(positive_path,'png',6,6,110)
table <- withCallingHandlers(ggplot2::ggplotGrob(positive), warning=function(w) {
  positive_warnings <<- c(positive_warnings,conditionMessage(w)); invokeRestart('muffleWarning')
})
positive_context <- attr(table,'ggchord_render_context')
grDevices::dev.off(); unlink(positive_path)
stopifnot(isTRUE(positive_context$converged),
  !any(grepl('physical geometry did not converge',positive_warnings,fixed=TRUE)))
report$positive <- list(geometry_converged=positive_context$converged,warnings=positive_warnings)
audit_assert_source_snapshot(audit_session)
jsonlite::write_json(report,file,pretty=TRUE,auto_unbox=TRUE,na='null')
cat('Actual render notice candidate checked; production sources unchanged\n')
