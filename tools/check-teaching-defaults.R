# Compare documented usages with current public function formals, without tests/checks.
pkgload::load_all(quiet = TRUE)
ns <- asNamespace("ggchord")
exports <- getNamespaceExports(ns)
results <- list()
for (path in list.files("man", "\\.Rd$", full.names = TRUE)) {
  rd <- tools::parse_Rd(path)
  usage <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\usage"), rd)
  if (!length(usage)) next
  code <- paste(unlist(usage[[1L]]), collapse = "")
  expressions <- tryCatch(parse(text = code), error = function(e) NULL)
  if (is.null(expressions)) next
  for (expression in expressions) {
    if (!is.call(expression) || !is.symbol(expression[[1L]])) next
    name <- as.character(expression[[1L]])
    if (!name %in% exports) next
    object <- getExportedValue("ggchord", name)
    if (!is.function(object) || !identical(environmentName(environment(object)), "ggchord")) next
    values <- as.list(expression)[-1L]
    argument_names <- names(values)
    if (is.null(argument_names)) argument_names <- rep("", length(values))
    for (i in which(!nzchar(argument_names))) {
      argument_names[i] <- as.character(values[[i]])
      values[i] <- alist(required = )
    }
    names(values) <- argument_names
    documented <- as.pairlist(values)
    results[[name]] <- list(rd = path, matches = identical(formals(object), documented))
  }
}
writeLines(jsonlite::toJSON(results, auto_unbox = TRUE, pretty = TRUE),
  "tools/teaching-audit/defaults-check.json")
failed <- names(Filter(function(x) !isTRUE(x$matches), results))
if (length(failed)) stop("Usage/formal mismatch: ", paste(failed, collapse = ", "))
cat(length(results), "documented exported function usages match current formals.\n")
