# Explicitly revisit the three measured noncircular profiles that hit the old cap.
script <- readLines('tools/long-term-audit/check-adaptive-quality.R', warn = FALSE)
cut <- grep('^profiles <- list', script)[1L]
eval(parse(text = script[seq_len(cut)]))
namespace <- asNamespace('ggchord')
candidate <- get('ggchord_render_table', envir = namespace)
changed <- 0L
rewrite <- function(expr) {
  if (!is.call(expr)) return(expr)
  if (identical(expr[[1L]], as.name('<-')) && identical(expr[[2L]], as.name('passes_allowed'))) {
    expr[[3L]] <- 10L; changed <<- changed+1L
    return(expr)
  }
  as.call(lapply(as.list(expr), rewrite))
}
body(candidate) <- rewrite(body(candidate))
stopifnot(changed == 1L)
unlockBinding('ggchord_render_table', namespace)
assign('ggchord_render_table', candidate, envir = namespace)
lockBinding('ggchord_render_table', namespace)
suffix <- '-noncircular-refinement-prototype'
cases <- list(collinear = c(6, 9), radial = c(9, 6), sankey = c(6, 9))
observations <- list()
for (id in names(cases)) {
  cat('Noncircular refinement candidate:', id, '\n'); flush.console()
  p <- plots[[id]] + adapt_ggchord('adaptive')
  dimensions <- cases[[id]]
  observations[[id]] <- measure(p, id, dimensions[1L], dimensions[2L], TRUE)
}
sources <- list.files('R', '[.]R$', full.names = TRUE)
jsonlite::write_json(list(scope = 'Three in-process candidate profiles only; source, explicit styles and convergence measured, no native certification',
  production_override = 'Use the existing circular ordinary refinement budget for the noncircular candidate',
  candidate_body_sha256 = digest::digest(deparse(body(candidate)), algo = 'sha256'),
  baseline_runtime_sha256 = as.list(setNames(vapply(sources,
    function(file) digest::digest(file = file, algo = 'sha256'), character(1)), sources)),
  observations = observations), 'tools/long-term-audit/v0.27.0-noncircular-refinement-prototype-results.json',
  pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Noncircular measured refinement candidates recorded\n')
