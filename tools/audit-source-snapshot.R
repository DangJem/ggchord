# Freeze the code actually loaded by an explicit audit, before rendering starts.
# This is development tooling and is not part of the installed package.
audit_source_snapshot <- function(additional = character()) {
  engine <- list.files('R', '[.]R$', full.names = TRUE)
  if (!length(engine)) stop('No R sources found in the project root', call. = FALSE)
  paths <- sort(unique(c(engine, 'DESCRIPTION', 'NAMESPACE', additional,
    'tools/audit-source-snapshot.R')))
  if (any(!file.exists(paths))) stop('An audit source input is missing', call. = FALSE)
  hashes <- vapply(paths, function(path) digest::digest(file = path, algo = 'sha256'), character(1))
  list(engine = hashes[engine], inputs = hashes, additional = additional)
}

audit_assert_source_snapshot <- function(snapshot) {
  current <- audit_source_snapshot(snapshot$additional)
  if (!identical(snapshot, current)) stop(
    'Audit source inputs changed during execution. Retain partial observations, but do not mark the run complete.',
    call. = FALSE)
  invisible(snapshot)
}
