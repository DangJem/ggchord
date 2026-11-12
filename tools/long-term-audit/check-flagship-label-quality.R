# Same published plot, annotations and settings; only physical size changes.
# Reuse the real-device measurement helper without executing its catalogue loop.
helper <- readLines('tools/long-term-audit/check-label-quality.R', warn = FALSE)
end <- grep('^# Explicit geometry controls', helper)[1L]-1L
stopifnot(is.finite(end))
eval(parse(text = helper[seq_len(end)]))
joint <- '--joint' %in% commandArgs(trailingOnly = TRUE)
prefix <- if (joint) 'v0.27.0-flagship-joint-' else 'v0.27.0-flagship-quality-'
routes <- '--routes' %in% commandArgs(trailingOnly = TRUE)
if (routes) prefix <- 'v0.27.0-flagship-routes-'
if ('--precise' %in% commandArgs(trailingOnly = TRUE)) {
  routes <- TRUE
  prefix <- 'v0.27.0-flagship-precise-'
}
if ('--target' %in% commandArgs(trailingOnly = TRUE)) {
  routes <- TRUE
  prefix <- 'v0.27.0-flagship-target-'
}
if ('--perimeter' %in% commandArgs(trailingOnly = TRUE)) {
  routes <- TRUE
  prefix <- 'v0.27.0-flagship-perimeter-'
}
if ('--styled' %in% commandArgs(trailingOnly = TRUE)) {
  routes <- TRUE
  prefix <- 'v0.27.0-flagship-styled-'
}
if ('--final' %in% commandArgs(trailingOnly = TRUE)) {
  routes <- TRUE
  prefix <- 'v0.27.0-flagship-final-'
}
if ('--ordered' %in% commandArgs(trailingOnly = TRUE)) {
  routes <- TRUE
  prefix <- 'v0.27.0-flagship-ordered-'
}
lines <- readLines('tools/documentation-examples.R', warn = FALSE)
markers <- grep('^# ---- [a-z]+ ----$', lines)
chunks <- stats::setNames(lapply(seq_len(length(markers)-1L), function(i)
  lines[seq.int(markers[i]+1L, markers[i+1L]-1L)]),
  sub('^# ---- ([a-z]+) ----$', '\\1', lines[markers[-length(markers)]]))
for (name in c('docstyle', 'showcaseplasmiddata', 'showcaseplasmid'))
  eval(parse(text = chunks[[name]]))
source <- p_show_plasmid$ggchord$data
layer_inputs <- lapply(p_show_plasmid$layers, function(layer) layer$ggchord_input_data)
profiles <- list(published = c(16.5, 8), matched = c(10, 8), square = c(8, 8), compact = c(6, 6))
selected <- sub('^--profiles=', '', grep('^--profiles=', commandArgs(trailingOnly = TRUE), value = TRUE))
if (length(selected)) {
  stopifnot(length(selected) == 1L)
  selected <- strsplit(selected, ',', fixed = TRUE)[[1L]]
  stopifnot(length(selected) > 0L, all(selected %in% names(profiles)))
  profiles <- profiles[selected]
}
results <- list()
for (name in names(profiles)) {
  cat('Flagship labels:', name, '\n'); flush.console()
  size <- profiles[[name]]
  results[[name]] <- measure(p_show_plasmid, name, size[1L], size[2L])
  stopifnot(identical(source, p_show_plasmid$ggchord$data),
    identical(layer_inputs, lapply(p_show_plasmid$layers, function(layer) layer$ggchord_input_data)))
}
sources <- list.files('R', '[.]R$', full.names = TRUE)
jsonlite::write_json(list(scope = 'Published pETDuet-1 plot from canonical documentation chunks; same annotation tables and fixed layer settings across the selected physical PNG sizes',
  source_annotations = list(features = nrow(show_plasmid_features), primers = nrow(show_plasmid_primers),
    selected_restriction_rows = nrow(show_plasmid_sites)), profiles = results,
  source_and_layer_inputs_preserved = TRUE, per_profile_manual_tuning = 0L,
  runtime_sha256 = as.list(stats::setNames(vapply(sources,
    function(path) digest::digest(file = path, algo = 'sha256'), character(1)), sources)),
  caveats = c('Contact counts and visibility reasons are observations, not assumed passes.',
    'Exact ink, non-polygon strokes, arbitrary canvases and publication device/platform matrix remain unverified.')),
  if (routes) paste0('tools/long-term-audit/', prefix, paste(names(profiles), collapse = '-'), '-results.json') else
  if (joint) 'tools/long-term-audit/v0.27.0-flagship-joint-results.json' else
    'tools/long-term-audit/v0.27.0-flagship-quality-results.json', pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Flagship controlled size observations recorded\n')
