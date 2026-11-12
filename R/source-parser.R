# Shared optional parser process boundary; source paths are data, never code.
ggchord_parse_external_source <- function(script_text, configuration, python, caller) {
  ggchord_conversion_require('jsonlite', caller)
  if (!is.character(python) || length(python) != 1L || is.na(python) || !nzchar(python))
    ggchord_stop(caller, ': python must name one existing interpreter')
  executable <- if (file.exists(python)) file.path(normalizePath(dirname(python)), basename(python)) else unname(Sys.which(python))
  if (!nzchar(executable)) ggchord_stop(caller, ': selected Python interpreter was not found')
  temporary <- tempfile('ggchord-source-'); dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  script <- file.path(temporary, 'parse.py'); config <- file.path(temporary, 'input.json')
  output <- file.path(temporary, 'output.json'); log <- file.path(temporary, 'parser.log')
  writeLines(script_text, script, useBytes = TRUE)
  jsonlite::write_json(configuration, config, auto_unbox = TRUE)
  status <- suppressWarnings(system2(executable, c(shQuote(script), shQuote(config), shQuote(output)), stdout = log, stderr = log))
  messages <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
  if (!identical(status, 0L) || !file.exists(output)) ggchord_stop(caller,
    ': optional source parser failed: ', paste(messages, collapse = '\n'))
  list(data = jsonlite::read_json(output, simplifyVector = FALSE), messages = messages)
}

