# Read R syntax with source references so the site can copy only the definitions
# a lesson uses, preserving the original formatting and assignment order.
chunks <- jsonlite::fromJSON(file("stdin"), simplifyVector = FALSE)
inspect <- function(code) {
  exprs <- parse(text = code, keep.source = TRUE)
  refs <- attr(exprs, "srcref")
  lapply(seq_along(exprs), function(i) {
    expr <- exprs[[i]]
    thunk <- function() NULL
    body(thunk) <- expr
    used <- codetools::findGlobals(thunk, merge = FALSE)
    defined <- codetools::findLocalsList(list(expr))
    # findGlobals treats every assignment as local, including x <- f(x)
    # and x$field <- value. Both forms require the earlier definition of x.
    if (is.call(expr) && identical(expr[[1]], as.name("<-"))) {
      lhs <- expr[[2]]
      rhs <- expr[[3]]
      old <- intersect(defined, all.names(rhs, functions = FALSE))
      if (is.call(lhs)) old <- union(old, defined)
      used$variables <- union(used$variables, old)
    }
    value <- expr
    assigned <- is.call(value) && identical(value[[1]], as.name("<-"))
    if (assigned) value <- value[[3]]
    while (is.call(value) && identical(value[[1]], as.name("+"))) value <- value[[2]]
    output <- if (is.symbol(value)) as.character(value) else if (is.call(value)) {
      name <- value[[1]]
      if (is.call(name) && identical(name[[1]], as.name("::"))) as.character(name[[3]])
      else if (is.symbol(name)) as.character(name) else ""
    } else ""
    list(output = output, assigned = assigned, start = refs[[i]][1], end = refs[[i]][3], defines = as.list(defined),
         needs = as.list(unique(c(used$variables, used$functions))),
         calls = as.list(used$functions))
  })
}
result <- lapply(chunks, inspect)
packages <- c("ggplot2", "dplyr", "purrr", "tidyr")
exports <- setNames(lapply(packages, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.list(getNamespaceExports(pkg)) else list()
}), packages)
cat(jsonlite::toJSON(list(chunks = result, packages = exports), auto_unbox = TRUE))
