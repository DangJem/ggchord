# One source-backed before/after gallery pair; no package release checks.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional = c('tools/documentation-examples.R',
  'tools/long-term-audit/check-identity-default.R'))
pkgload::load_all(quiet = TRUE)
text <- readLines('tools/documentation-examples.R')
markers <- grep('^# ---- [a-z]+ ----$', text)
chunks <- setNames(lapply(seq_len(length(markers)-1L), function(i)
  text[seq.int(markers[i]+1L,markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$','\\1',text[markers[-length(markers)]]))
env <- new.env(parent = globalenv())
for (name in c('setup','sequencepalette','docstyle','showcasedata','learnseq','learnlinks'))
  eval(parse(text = chunks[[name]]), env)
current <- env$p_learn_links + coord_collinear()
# Extract the exact earlier palette/guide rather than inventing a baseline.
old <- system2('git', c('show', 'bffa858:tools/documentation-examples.R'), stdout=TRUE)
at <- grep('^# ---- docstyle ----$',old); end <- grep('^# ---- chord ----$',old)
oldenv <- new.env(parent=env)
eval(parse(text=old[seq.int(at+1L,end-1L)]),oldenv)
stopifnot(exists('identity_palette',oldenv,inherits=FALSE))
baseline <- current + scale_ribbon_fill_gradientn(
  colours=oldenv$identity_palette,values=oldenv$identity_values,
  limits=c(0,100),breaks=oldenv$identity_breaks,name='Identity (%)',
  guide=oldenv$identity_colourbar)
observed <- list()
for (stage in c('before','after')) {
  plot <- if(stage=='before') baseline else current
  path <- paste0('tools/long-term-audit/v0.27.0-identity-',stage,'.png')
  ggchord_diagnostic_device(path,'png',12,8,144)
  table <- ggplot2::ggplotGrob(plot)
  grid::grid.newpage(); grid::grid.draw(table)
  built <- attr(table,'ggchord_render_built')
  scale <- built$plot$scales$get_scales('ribbon_fill')
  stopifnot(!is.null(scale))
  samples <- c(1,25,49,51,79,81,89,91,94,96,99)
  colours <- scale$map(samples)
  observed[[stage]] <- list(file=path,breaks=scale$get_breaks(),limits=scale$get_limits(),
    samples=samples,colours=colours,converged=attr(table,'ggchord_render_context')$converged)
  if(stage=='after') {
    stopifnot(identical(scale$get_limits(),c(0,100)),
      identical(as.numeric(scale$get_breaks()),c(0,50,80,90,95,100)),
      length(unique(colours[1:3]))==1L, length(unique(colours))==5L)
  } else stopifnot(length(unique(colours[1:3]))>1L)
  grDevices::dev.off()
}
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(cases=observed,runtime_sha256=as.list(audit$engine),
  input_sha256=as.list(audit$inputs),scope='Same real HSPs, 12x8in/144dpi; public binned mapping and explicit continuous override'),
  'tools/long-term-audit/v0.27.0-identity-results.json',pretty=TRUE,auto_unbox=TRUE)
cat('Identity default and explicit continuous override checked\n')
