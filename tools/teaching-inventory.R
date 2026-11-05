# Inspect current exports and signatures without installing or checking the package.
pkgload::load_all(quiet = TRUE)
ns <- asNamespace("ggchord")
exports <- sort(getNamespaceExports(ns))
objects <- lapply(exports, function(name) getExportedValue("ggchord", name))
rows <- lapply(seq_along(exports), function(index) {
  name <- exports[index]
  object <- getExportedValue("ggchord", name)
  owner <- if (is.function(object)) environmentName(environment(object)) else "object"
  defaults <- if (is.function(object)) lapply(formals(object), function(x)
    paste(deparse(x), collapse = " ")) else list()
  same <- which(vapply(objects, function(other) identical(object, other), logical(1)))
  aliases <- setdiff(exports[same], name)
  list(name = name, owner = owner, alias = paste(aliases, collapse = ";"), parameters = defaults)
})
cat(jsonlite::toJSON(rows, auto_unbox = TRUE, pretty = TRUE))
