# Locate a real new drawing notice without replacing repository PNGs/manifests.
source('tools/audit-source-snapshot.R')
refinement_observer_snapshot <- audit_source_snapshot(c(
  'tools/long-term-audit/observe-documentation-refinement.R',
  'tools/documentation-plots.R','tools/documentation-examples.R',
  'tools/render-teaching-comparisons.R','tools/documentation-fingerprints.py',
  'website/parameter-comparisons.json'))
script <- readLines('tools/documentation-plots.R',warn=FALSE)
cut <- grep('^for\\(i in selected\\)',script)[1L]-1L
stopifnot(is.finite(cut))
preparation_warnings <- list()
withCallingHandlers(eval(parse(text=script[seq_len(cut)])),warning=function(w) {
  preparation_warnings[[length(preparation_warnings)+1L]] <<- list(
    snippet=if(exists('name',envir=globalenv(),inherits=FALSE)) get('name',envir=globalenv()) else NA_character_,
    warning=conditionMessage(w))
  invokeRestart('muffleWarning')
})
observations <- list();found <- FALSE
for(i in selected) {
  id <- assets$id[i];messages <- character()
  cat('Refinement attribution:',id,'\n');flush.console()
  plot <- get(assets$object[i],env)
  preview <- withCallingHandlers(view_ggchord(plot,width=assets$width[i],height=assets$height[i],
    units='in',dpi=assets$dpi[i],device='png',bg=assets$background[i],viewer='none'),warning=function(w) {
      messages <<- c(messages,conditionMessage(w));invokeRestart('muffleWarning')
    })
  source_file <- paste0('man/figures/',id,'.png')
  observations[[id]] <- list(width=assets$width[i],height=assets$height[i],dpi=assets$dpi[i],
    warnings=messages,original_png_sha256=digest::digest(file=source_file,algo='sha256'),
    observed_png_sha256=digest::digest(file=preview,algo='sha256'))
  found <- any(grepl('physical geometry did not converge',messages,fixed=TRUE))
  if(found) {cat('Native refinement warning located:',id,'\n');break}
}
audit_assert_source_snapshot(refinement_observer_snapshot)
jsonlite::write_json(list(scope='Explicit warning attribution; unchanged canonical preparations and PNG dimensions. Repository PNGs/manifests are not written. Stops at the first observed physical-refinement warning; not a layout-quality pass.',
  original_input_preparation_warnings=preparation_warnings,
  observations=observations,native_warning_found=found,
  full_figure_set_observed=length(observations)==length(selected),
  runtime_sha256=as.list(refinement_observer_snapshot$engine)),
  'tools/long-term-audit/v0.27.0-documentation-refinement-attribution.json',pretty=TRUE,auto_unbox=TRUE,na='null')
cat('Refinement attribution observations recorded\n')
