#!/usr/bin/env Rscript
# Static API inventory for ggchord. No package installation, execution,
# network access, repository mutation, or publish step is performed.
# Usage: Rscript audit_ggchord_api.R /path/to/ggchord /path/to/output
# Review the generated tables manually; zero lexical references DOES NOT prove
# a function is unused (ggproto, S3, eval(quote()), NSE, strings, and callback
# dispatch require additional investigation).

args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1L) normalizePath(args[[1]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
out <- if (length(args) >= 2L) args[[2]] else file.path(repo, "design-roadmap", "audit-generated")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
nsfile <- file.path(repo, "NAMESPACE")
if (!file.exists(nsfile)) stop("NAMESPACE is missing; run this from ggchord checkout")
ns <- readLines(nsfile, warn = FALSE)
strip <- function(pattern) sub(pattern, "\\1", grep(pattern, ns, value = TRUE))
exports <- strip("^export\\(([^)]*)\\)$")
method_records <- grep("^S3method\\(", ns, value = TRUE)
import_records <- grep("^importFrom\\(", ns, value = TRUE)
imports <- do.call(rbind, lapply(import_records, function(line) {
  m <- regmatches(line, regexec("^importFrom\\(([^,]*),([^)]*)\\)", line))[[1]]
  if (length(m) != 3L) return(data.frame(package = "", name = ""))
  data.frame(package = m[[2]], name = m[[3]], stringsAsFactors = FALSE)
}))
if (is.null(imports)) imports <- data.frame(package = character(), name = character())
source_files <- sort(list.files(file.path(repo, "R"), pattern = "\\.R$", full.names = TRUE))
code <- list(); formal_records <- list(); defs <- list()
format_expr <- function(x) {
  s <- paste(deparse(x, width.cutoff = 400L), collapse = " ")
  if (nchar(s) > 1000L) paste0(substr(s, 1L, 997L), "...") else s
}
walk <- function(expr, file, scope = "<top-level>", depth = 0L) {
  if (!is.call(expr)) return(invisible(NULL))
  h <- as.character(expr[[1]])[1]
  if (h %in% c("<-", "=", "<<-") && length(expr) == 3L &&
      is.symbol(expr[[2]]) && is.call(expr[[3]]) &&
      identical(as.character(expr[[3]][[1]])[1], "function")) {
    name <- as.character(expr[[2]])
    fullname <- if (scope == "<top-level>") name else paste0(scope, "::", name)
    fn <- expr[[3]]
    sr <- attr(expr, "srcref")
    line <- if (is.null(sr)) NA_integer_ else as.integer(sr[[1]])
    arglist <- as.list(fn[[2]])
    defs[[length(defs) + 1L]] <<- data.frame(
      name = name, qualified_name = fullname,
      source_file = file, source_line = line,
      scope = scope, is_top_level = scope == "<top-level>",
      is_exported = name %in% exports && scope == "<top-level>",
      parameter_count = length(arglist),
      calls_lexical = paste(sort(unique(all.names(fn[[3]], functions = TRUE))), collapse = ";"),
      stringsAsFactors = FALSE)
    if (length(arglist)) for (j in seq_along(arglist)) {
      default <- tryCatch({
        raw_default <- paste(deparse(arglist[[j]], width.cutoff = 400L), collapse = " ")
        if (!nzchar(trimws(raw_default))) "<required>" else raw_default
      }, error = function(e) "<required>")
      formal_records[[length(formal_records) + 1L]] <<- data.frame(
        function_name = name, qualified_function = fullname,
        source_file = file, source_line = line,
        is_exported = name %in% exports && scope == "<top-level>",
        parameter_index = j, parameter = names(arglist)[[j]],
        default_expr = default, stringsAsFactors = FALSE)
    }
    walk(fn[[3]], file, scope = fullname, depth = depth + 1L)
    return(invisible(NULL))
  }
  if (h == "function") {
    # Anonymous functions are relevant but usually not stable names.
    walk(expr[[3]], file, scope = paste0(scope, "::<anonymous>"), depth = depth + 1L)
    return(invisible(NULL))
  }
  if (length(expr) >= 2L) for (i in seq.int(2L, length(expr))) walk(expr[[i]], file, scope = scope, depth = depth + 1L)
  invisible(NULL)
}
for (file in source_files) {
  rel <- substring(file, nchar(repo) + 2L)
  parsed <- tryCatch(parse(file = file, keep.source = TRUE), error = function(e) {
    warning("Cannot parse ", rel, ": ", conditionMessage(e)); NULL
  })
  if (is.null(parsed)) next
  for (expr in parsed) walk(expr, rel)
}
functions <- if (length(defs)) do.call(rbind, defs) else data.frame()
arguments <- if (length(formal_records)) do.call(rbind, formal_records) else data.frame()
if (nrow(functions)) functions <- functions[order(functions$source_file, functions$source_line, functions$qualified_name), ]
if (nrow(arguments)) arguments <- arguments[order(arguments$source_file, arguments$source_line, arguments$qualified_function, arguments$parameter_index), ]
reexport <- exports %in% imports$name & exports %in% imports$name[imports$package %in% c("ggplot2", "grid")]
export_table <- data.frame(symbol = exports, is_third_party_reexport = reexport,
  namespace_source = "NAMESPACE", stringsAsFactors = FALSE)
if (nrow(functions)) {
  export_table$definition_candidates <- vapply(exports, function(name) paste(
    unique(functions$source_file[functions$name == name & functions$is_top_level]), collapse = ";"), "")
} else export_table$definition_candidates <- ""
export_table$review_decision <- "UNREVIEWED"
export_table$review_evidence <- ""
export_table$future_public_contract <- "UNDECIDED"
s3 <- data.frame(registration = method_records, stringsAsFactors = FALSE)
write.csv(export_table, file.path(out, "api_exports.csv"), row.names = FALSE, na = "")
write.csv(imports, file.path(out, "api_imports.csv"), row.names = FALSE, na = "")
write.csv(s3, file.path(out, "api_s3.csv"), row.names = FALSE, na = "")
write.csv(functions, file.path(out, "api_functions.csv"), row.names = FALSE, na = "")
write.csv(arguments, file.path(out, "api_arguments.csv"), row.names = FALSE, na = "")
# Heuristic only: lexical names of other functions may be found in quoted code,
# and callbacks accessed by string will be missed. DO NOT use for auto-delete.
if (nrow(functions)) {
  roots <- functions[functions$is_top_level, c("name", "source_file", "calls_lexical")]
  lexical <- unique(unlist(strsplit(roots$calls_lexical, ";", fixed = TRUE), use.names = FALSE))
  usage <- unique(roots[, c("name", "source_file")])
  usage$lexically_mentioned_in_function_bodies <- usage$name %in% lexical
  usage$is_exported <- usage$name %in% exports
  usage$human_review <- "REQUIRED"
  write.csv(usage, file.path(out, "api_lexical_reference_HEURISTIC.csv"), row.names = FALSE, na = "")
}
cat(sprintf("Repository: %s\nR source files: %d\nExports: %d\nThird-party reexports: %d\nS3 registrations: %d\nNamed function definitions (including nested): %d\nNamed formal arguments: %d\nOutput: %s\n",
  repo, length(source_files), length(exports), sum(reexport), length(method_records),
  nrow(functions), nrow(arguments), normalizePath(out)))
cat("IMPORTANT: Review aliases, ggproto objects, dynamic wrappers, ellipsis, quote/eval and callbacks manually.\n")
