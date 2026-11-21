# Confirm that every published invocation has the actual callable's formals.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1L]])
pkgload::load_all(root, export_all = FALSE, quiet = TRUE)
manual <- jsonlite::read_json(file.path(root, "website/reference-manual/structure.json"))
rows <- lapply(names(manual$topics), function(topic) {
  entry <- manual$topics[[topic]]
  package <- if (entry$native) "ggplot2" else "ggchord"
  actual <- formals(getExportedValue(package, topic))
  defaults <- entry$defaults
  parameters <- names(defaults)
  terms <- mapply(function(key, value) if (nzchar(value)) paste(key, "=", value) else key,
                  parameters, defaults, USE.NAMES = FALSE)
  documented <- parse(text = paste0("function(", paste(terms, collapse = ", "), ") NULL"))[[1L]][[2L]]
  if (!identical(actual, documented)) stop("Reference signature differs from ", package, "::", topic)
  if (!length(parameters)) return(NULL)
  data.frame(api = topic, parameter = parameters, default = unlist(defaults, use.names = FALSE),
             topic = topic, anchor = ifelse(parameters == "...", "para-ellipsis", paste0("parameter-", parameters)))
})
jsonlite::write_json(do.call(rbind, rows), args[[2L]], dataframe = "rows", auto_unbox = TRUE)
