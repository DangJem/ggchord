# Refresh the single registered PNG of each potentially affected Chord/Radial
# figure. Other coordinates keep their original files and provenance.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional = c('tools/documentation-examples.R',
  'tools/site-audit/refresh-ribbon-figures.R','tools/documentation-fingerprints.py',
  list.files('data','[.]rda$',full.names=TRUE), 'R/sysdata.rda'))
pkgload::load_all(quiet = TRUE)
lines <- readLines('tools/documentation-examples.R', warn = FALSE)
markers <- grep('^# ---- [a-z]+ ----$',lines)
chunks <- setNames(lapply(seq_len(length(markers)-1L),function(i)
  lines[seq.int(markers[i]+1L,markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$','\\1',lines[markers[-length(markers)]]))
env <- new.env(parent = globalenv())
for (key in setdiff(names(chunks),c('save','preview','frontsave','frontpreview','fusionsave','imports','firstsave','firstsavebasic'))) {
  cat('Prepare:',key,'\n')
  invisible(eval(parse(text=chunks[[key]]),env))
}
env$p_labelauto <- env$label_plots$auto
env$p_labelradial <- env$label_plots$radial
env$p_labelarc <- env$label_plots$arc
assets <- read.csv('website/assets/plots.csv',stringsAsFactors=FALSE)
selected <- which(assets$coordinate %in% c('chord','radial_synteny'))
stopifnot(length(selected)>0L)
fingerprints <- jsonlite::fromJSON(system2('python3','tools/documentation-fingerprints.py',stdout=TRUE))
result <- list(status='in progress',runtime_sha256=as.list(audit$engine),
  inputs_sha256=as.list(audit$inputs),figures=list(),
  scope='Single normal-size PNG per registered Chord/Radial figure; no release device matrix. Other coordinates retained under the scoped style change.')
path <- 'tools/site-audit/ribbon-figure-refresh.json'
for (i in selected) {
  id <- assets$id[i];cat('Refresh:',id,'\n');flush.console()
  before <- unname(tools::md5sum(assets$file[i]))
  plot <- get(assets$object[i],env)
  preview <- view_ggchord(plot,width=assets$width[i],height=assets$height[i],
    units='in',dpi=assets$dpi[i],device='png',bg=assets$background[i],viewer='none')
  stopifnot(file.copy(preview,assets$file[i],overwrite=TRUE))
  audit_assert_source_snapshot(audit)
  assets$md5[i] <- unname(tools::md5sum(assets$file[i]))
  assets$code_md5[i] <- fingerprints$code
  assets$dependency_md5[i] <- fingerprints$assets[[id]]
  assets$generator[i] <- 'tools/site-audit/refresh-ribbon-figures.R'
  result$figures[[id]] <- list(file=assets$file[i],before_md5=before,md5=assets$md5[i],
    sha256=digest::digest(file=assets$file[i],algo='sha256',serialize=FALSE),
    width=assets$width[i],height=assets$height[i],dpi=assets$dpi[i],
    code_md5=assets$code_md5[i],dependency_md5=assets$dependency_md5[i],status='actual PNG rendered')
  write.csv(assets,'website/assets/plots.csv',row.names=FALSE)
  jsonlite::write_json(result,path,auto_unbox=TRUE,pretty=TRUE)
}
result$status <- 'complete'
jsonlite::write_json(result,path,auto_unbox=TRUE,pretty=TRUE)
cat('Refreshed',length(selected),'actual single-PNG figures\n')
