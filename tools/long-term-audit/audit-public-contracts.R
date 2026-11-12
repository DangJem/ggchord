# Current public surface only. Never rewrite historical frozen candidates or
# infer semantic aliases just because two R objects compare identical.
pkgload::load_all(quiet = TRUE)
listed <- jsonlite::read_json('tools/teaching-audit/api-signatures.json')
exports <- getNamespaceExports('ggchord')
stopifnot(setequal(exports, vapply(listed, `[[`, '', 'name')))
canonical <- function(text) if (!nzchar(text)) '<required-or-dots>' else
  paste(deparse(parse(text = text)[[1L]], width.cutoff = 500), collapse = ' ')
rows <- list(); differences <- character()
for (entry in listed) {
  object <- getExportedValue('ggchord', entry$name)
  arguments <- if (is.function(object)) formals(object) else list()
  actual <- vapply(seq_along(arguments), function(i)
    paste(deparse(arguments[[i]], width.cutoff = 500), collapse = ' '), character(1))
  expected <- vapply(entry$parameters, canonical, character(1))
  local <- identical(entry$owner, 'ggchord')
  matched <- identical(names(arguments), names(expected)) &&
    identical(unname(vapply(actual, canonical, character(1))), unname(expected))
  if (local && !matched) differences <- c(differences, entry$name)
  upstream <- if (entry$owner == 'object') 'ggplot2' else entry$owner
  upstream_equal <- if (!local) identical(object, getExportedValue(upstream, entry$name)) else NA
  if (!local) stopifnot(upstream_equal)
  rows[[length(rows)+1L]] <- data.frame(name = entry$name, owner = entry$owner,
    kind = if (is.function(object)) 'function' else 'upstream object',
    parameter_count = length(arguments),
    signature = paste(paste0(names(arguments), '=', actual), collapse = '; '),
    local_documented_signature_matches = if (local) matched else NA,
    exact_upstream_binding = upstream_equal)
}
stopifnot(length(differences) == 0L, length(formals(theme_ggchord)) == 109L,
  identical(names(formals(ggchord))[1:3], c('seq_data', 'ribbon_data', 'gene_data')))
# Inspect declaration syntax, not arbitrary prose, defaults or third-party calls.
brands <- 'snapgene|geneious|circos|clinker|gggenomes|genoplotr|gggenes|circlize|ggview|pfam|interpro|mcscanx|mummer|odgi|ggplot2'
branded <- character()
scan <- function(expr, path) {
  if (!is.call(expr) && !is.expression(expr) && !is.pairlist(expr)) return()
  if (is.call(expr) && identical(expr[[1L]], as.name('function'))) {
    found <- names(expr[[2L]])
    found <- found[grepl(brands, found, ignore.case = TRUE)]
    if (length(found)) branded <<- c(branded, paste(path, found, sep = ':'))
  }
  if (is.call(expr) && as.character(expr[[1L]])[1L] %in% c('<-', '=' ) && length(expr) == 3L &&
      is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name('function'))) {
    name <- paste(deparse(expr[[2L]]), collapse = '')
    if (grepl(brands, name, ignore.case = TRUE)) branded <<- c(branded, paste(path, name, sep = ':'))
  }
  # Pairlists contain missing formal arguments; recurse only through call bodies.
  if (!is.pairlist(expr)) for (i in seq_along(expr)) {
    if (is.call(expr) && i == 1L) next
    if (is.call(expr[[i]]) || is.expression(expr[[i]])) scan(expr[[i]], path)
  }
}
for (path in list.files('R', '[.]R$', full.names = TRUE)) scan(parse(path), path)
stopifnot(!length(branded))
folder <- 'tools/long-term-audit'
write.csv(do.call(rbind, rows), file.path(folder, 'v0.27.0-public-contracts.csv'), row.names = FALSE)
jsonlite::write_json(list(status = 'passed signature/source-surface review only',
  base_commit = system2('git', c('rev-parse', 'HEAD'), stdout = TRUE),
  working_tree_note = 'Development review can include uncommitted implementation; source hashes, not HEAD alone, identify loaded inputs.',
  source_sha256 = as.list(stats::setNames(vapply(list.files('R', '[.]R$', full.names = TRUE),
    function(path) digest::digest(file = path, algo = 'sha256'), character(1)),
    list.files('R', '[.]R$', full.names = TRUE))),
  exports = length(exports), local_exports = sum(vapply(listed, function(x) x$owner == 'ggchord', logical(1))),
  theme_formals = length(formals(theme_ggchord)), normalized_signature_drift = differences,
  upstream_bindings_unchanged = TRUE, old_constructor_third_position_retained = TRUE,
  third_party_branded_declarations = branded,
  caveats = c('Introspection is not functional validation or a 1.0 API freeze.',
    'Whitespace from deparse is normalized through parsed expressions, not treated as a semantic API change.',
    'No semantic alias is inferred from object equality; after_stat and after_scale retain distinct staged meanings.',
    'Historical frozen snapshots and release checks are not rewritten or executed.')),
  file.path(folder, 'v0.27.0-public-contract-results.json'), pretty = TRUE, auto_unbox = TRUE)
cat(length(exports), 'public bindings match current inventory; 109 theme formals retained\n')
