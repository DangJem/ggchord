# Targeted diagnosis of declared documentation comparisons; no source changes.
source('tools/audit-source-snapshot.R')
snapshot <- audit_source_snapshot(c('tools/long-term-audit/inspect-comparison-quality.R',
  'tools/documentation-examples.R','website/parameter-comparisons.json'))
pkgload::load_all(quiet=TRUE)
data(plasmid_example_pUC19)
registry <- jsonlite::fromJSON('website/parameter-comparisons.json',simplifyVector=FALSE)
requested <- commandArgs(trailingOnly=TRUE)
if (!length(requested)) requested <- 'compareexternaloff'
stopifnot(all(requested %in% vapply(registry,`[[`,character(1),'id')))
lines <- readLines('tools/documentation-examples.R',warn=FALSE)
markers <- grep('^# ---- [a-z]+ ----$',lines)
chunks <- setNames(lapply(seq_len(length(markers)-1L),function(i)
  lines[seq.int(markers[i]+1L,markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$','\\1',lines[markers[-length(markers)]]))
observations <- list()
for(id in requested) {
  row <- Filter(function(x) x$id==id,registry)[[1L]]
  environment <- new.env(parent=globalenv())
  for(key in row$snippets) eval(parse(text=chunks[[key]]),environment)
  plot <- environment[[row$object]]
  original <- environment$finish_labels
  cat('Comparison diagnosis:',id,'\n');flush.console()
  result <- diagnose_ggchord(plot,width=row$width,height=row$height,dpi=144,device='png')
  stopifnot(identical(original,environment$finish_labels))
  observations[[id]] <- list(settings=row,measurement=result$measurement,
    checks=result$checks,warnings=result$warnings,label_quality=result$label_quality,
    visibility=result$label_visibility,
    source_labels_preserved=TRUE,documented_png_sha256=digest::digest(
      file=paste0('man/figures/',id,'.png'),algo='sha256'))
}
audit_assert_source_snapshot(snapshot)
jsonlite::write_json(list(scope='Targeted original comparison settings; report observations, not assumed quality passes',
  observations=observations,runtime_sha256=as.list(snapshot$engine)),
  paste0('tools/long-term-audit/v0.27.0-comparison-quality-',paste(requested,collapse='-'),'.json'),
  pretty=TRUE,auto_unbox=TRUE,na='null')
cat('Targeted comparison observations recorded\n')
