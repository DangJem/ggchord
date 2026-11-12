# Same original three progressive plots; canvas size is the only change.
helper <- readLines('tools/long-term-audit/check-label-quality.R',warn=FALSE)
cut <- grep('^# Explicit geometry controls',helper)[1L]-1L
helper <- helper[seq_len(cut)]
stopifnot(sum(grepl('width, height, 110',helper,fixed=TRUE))==1L)
helper <- sub('width, height, 110','width, height, 144',helper,fixed=TRUE)
eval(parse(text=helper))
audit_session <- audit_source_snapshot(c(audit_session$additional,
  'tools/long-term-audit/check-linear-canvas.R','tools/long-term-audit/linear-lesson-source.R'))
source('tools/long-term-audit/linear-lesson-source.R')
lesson <- load_linear_lesson()
original <- list(features=lesson$show_features,links=lesson$show_links,sequences=lesson$show_sequences)
prefix <- 'v0.27.0-linear-canvas-'
objects <- c(linearseq='p_linear_seq',linearlinks='p_linear_links',linearlabels='p_linear_labels')
observations <- list()
for(id in names(objects)) {
  cat('Common linear canvas:',id,'12 x 8\n');flush.console()
  warnings <- character()
  result <- withCallingHandlers(measure(lesson[[objects[id]]],paste0(id,'-12x8'),12,8),warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')
  })
  result$warnings <- warnings
  observations[[id]] <- result
}
stopifnot(identical(original,list(features=lesson$show_features,links=lesson$show_links,sequences=lesson$show_sequences)))
audit_assert_source_snapshot(audit_session)
jsonlite::write_json(list(scope='Original progressive linear plots at a common 12x8 inch/144 dpi canvas; no layer/data/font/guide/ratio changes. The original 10x6 failure remains separate.',
  observations=observations,source_tables_preserved=TRUE,per_case_manual_tuning=0L,
  runtime_sha256=as.list(audit_session$engine)),
  'tools/long-term-audit/v0.27.0-linear-canvas-12x8-results.json',pretty=TRUE,auto_unbox=TRUE,na='null')
cat('Common linear canvas observations recorded\n')
