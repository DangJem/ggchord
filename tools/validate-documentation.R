# Scope: current documentation imports and one explicit-size PNG preview/diagnostic.
# Does not run testthat, build/check the package, knit vignettes or test device matrices.
pkgload::load_all(quiet = TRUE)
text <- readLines('tools/documentation-examples.R', warn = FALSE)
markers <- grep('^# ---- [a-z]+ ----$', text)
chunks <- setNames(lapply(seq_len(length(markers)-1L),function(i)
  text[seq.int(markers[i]+1L, markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$', '\\1', text[markers[-length(markers)]]))
env <- new.env(parent = globalenv())
for (name in c('setup','sequencepalette','docstyle','collinear','imports')) eval(parse(text=chunks[[name]]),env)
stopifnot(nrow(env$seq_files)>0L, nrow(env$link_files)>0L, nrow(env$feature_files)>0L,
  all(env$link_files$qaccver %in% env$seq_files$accver),
  all(env$link_files$saccver %in% env$seq_files$accver))
image <- view_ggchord(env$p_collinear,width=8,height=5,units='in',viewer='none')
stopifnot(file.exists(image),file.info(image)$size>0)
report <- diagnose_ggchord(env$p_collinear,width=8,height=5,units='in',device='png')
stopifnot(inherits(report,'ggchord_diagnostics'),is.list(report$measurement))
record <- list(scope='Documentation only: three original-file readers and one target-size PNG preview/diagnostic',
  source='tools/documentation-examples.R', code_md5=unname(tools::md5sum('tools/documentation-examples.R')),
  R=R.version.string, ggplot2=as.character(utils::packageVersion('ggplot2')),
  input_rows=list(sequences=nrow(env$seq_files),links=nrow(env$link_files),features=nrow(env$feature_files)),
  preview=list(width=8,height=5,units='in',device='png',result='PASS'),
  panel_status=report$measurement$status,
  limitations='No complete vignette rendering, package check or release-device matrix; measured panel is not visual certification')
writeLines(jsonlite::toJSON(record,auto_unbox=TRUE,pretty=TRUE),'website/verification.json')
cat('PASS: documented FASTA/BLAST/GFF3 imports and one 8 x 5 in PNG preview/diagnostic.\n')
