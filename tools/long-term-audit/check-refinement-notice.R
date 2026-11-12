# Actual native PNG failure notice plus one converged public plot.
script <- readLines('tools/long-term-audit/check-flagship-label-quality.R', warn = FALSE)
cut <- grep('^results <- list', script)[1L]
eval(parse(text = script[seq_len(cut)]))
audit_assert_source_snapshot(audit_session)
audit_session <- audit_source_snapshot(c(audit_session$additional, 'tools/long-term-audit/check-refinement-notice.R'))
prefix <- 'v0.27.0-refinement-notice-native-'
notices <- character()
withCallingHandlers(eval(parse(text = script[-seq_len(cut)])), warning = function(w) {
  notices <<- c(notices, conditionMessage(w)); invokeRestart('muffleWarning')
})
file <- paste0('tools/long-term-audit/',prefix,paste(names(profiles),collapse='-'),'-results.json')
report <- jsonlite::fromJSON(file, simplifyVector=FALSE)
stopifnot(length(grep('physical geometry did not converge',notices,fixed=TRUE)) == 1L,
  identical(report$profiles$compact$geometry_converged,FALSE),
  report$profiles$compact$n_visible_source_labels == 164L)
report$native_notice_verified <- TRUE
report$production_override <- FALSE
report$warnings <- notices
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
cat('Native drawing notice verified on a real failure and a converged plot\n')
