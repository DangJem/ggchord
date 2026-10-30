#!/usr/bin/env Rscript
# Static source inventory only. Never loads/executes the package or deletes API.
# Usage: Rscript tools/audit-ggchord-api.R [repo] [output] [audit-start-revision]
# Derived from the supplied v0.15 audit handoff, extended for anonymous/ggproto
# methods, aliases, dynamic construction, documentation and coverage accounting.
options(expressions = 50000L)
argv <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(argv)) argv[[1L]] else ".", mustWork = TRUE)
out <- if (length(argv) >= 2L) argv[[2L]] else file.path(repo, "design-roadmap/audit-generated")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
git <- function(...) {
  x <- suppressWarnings(system2("git", c("-C", shQuote(repo), ...), stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(x, "status"))) stop("Cannot read Git metadata: ", paste(x, collapse = " "))
  x
}
head_revision <- git("rev-parse", "HEAD")[[1L]]
audit_start <- if (length(argv) >= 3L) argv[[3L]] else head_revision
if (!identical(git("rev-parse", shQuote(audit_start))[[1L]], audit_start))
  stop("audit-start-revision must be a full resolvable commit SHA")
relpath <- function(path) substring(path, nchar(repo) + 2L)
source_input_paths <- function() c(
  file.path(repo, "DESCRIPTION"), file.path(repo, "NAMESPACE"),
  sort(list.files(file.path(repo, "R"), "\\.R$", full.names = TRUE)),
  sort(list.files(file.path(repo, "man"), "\\.Rd$", full.names = TRUE)),
  sort(list.files(file.path(repo, "design-roadmap/audit-reviewed"), "\\.csv$", full.names = TRUE)),
  sort(list.files(file.path(repo, "design-roadmap/audit-reviewed-details"), "\\.csv$", full.names = TRUE))
)
# Freeze the input set before the first read. A concurrent edit must abort the
# run before output is written, rather than attach new hashes to an older AST.
source_files <- source_input_paths()
source_hashes_before <- tools::md5sum(source_files)
if (anyNA(source_hashes_before)) stop("Cannot fingerprint all source inputs; rerun after resolving missing inputs")
source_r_files <- source_files[startsWith(source_files, paste0(file.path(repo, "R"), "/"))]
source_rd_files <- source_files[startsWith(source_files, paste0(file.path(repo, "man"), "/"))]
review_files <- source_files[startsWith(source_files, paste0(file.path(repo, "design-roadmap/audit-reviewed"), "/"))]
script_path <- file.path(repo, "tools/audit-ggchord-api.R")
generator_hash_before <- unname(tools::md5sum(script_path))
package_version_before <- read.dcf(file.path(repo, "DESCRIPTION"))[1L, "Version"]
verify_snapshot <- function() {
  after_paths <- source_input_paths()
  if (!identical(source_files, after_paths) ||
      !identical(source_hashes_before, tools::md5sum(after_paths)) ||
      !identical(head_revision, git("rev-parse", "HEAD")[[1L]]) ||
      !identical(generator_hash_before, unname(tools::md5sum(script_path))))
    stop("Audit input files, contents, HEAD or generator changed during parsing; no output written. Rerun with a stable source snapshot.")
  invisible(TRUE)
}
fmt <- function(x, max_chars = Inf) {
  s <- paste(deparse(x, width.cutoff = 400L), collapse = " ")
  if (nchar(s) > max_chars) paste0(substr(s, 1L, max_chars), " [TRUNCATED]") else s
}
head_name <- function(x) if (is.symbol(x)) as.character(x) else fmt(x, 160L)
call_name <- function(x) {
  if (!is.call(x)) return("")
  h <- x[[1L]]
  if (is.symbol(h)) return(as.character(h))
  if (is.call(h) && head_name(h[[1L]]) %in% c("::", ":::"))
    return(paste0(as.character(h[[2L]]), "::", as.character(h[[3L]])))
  fmt(h, 160L)
}
leaf_name <- function(x) sub("^.*::", "", x)
row <- function(...) data.frame(..., stringsAsFactors = FALSE, check.names = FALSE)
records <- new.env(parent = emptyenv())
for (name in c("functions", "arguments", "objects", "calls", "references", "dynamic", "aesthetics"))
  records[[name]] <- list()
add <- function(bucket, ...) records[[bucket]][[length(records[[bucket]]) + 1L]] <- row(...)
collapse_records <- function(name, columns = character()) {
  x <- records[[name]]
  if (length(x)) return(do.call(rbind, x))
  as.data.frame(setNames(rep(list(character()), length(columns)), columns), check.names = FALSE)
}
write_table <- function(x, name) {
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "", fileEncoding = "UTF-8")
}
ns <- readLines(file.path(repo, "NAMESPACE"), warn = FALSE)
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
imports <- do.call(rbind, lapply(grep("^importFrom\\(", ns, value = TRUE), function(x) {
  z <- strsplit(sub("^importFrom\\((.*)\\)$", "\\1", x), ",", fixed = TRUE)[[1L]]
  row(package = z[[1L]], symbol = z[[2L]])
}))
s3 <- do.call(rbind, lapply(grep("^S3method\\(", ns, value = TRUE), function(x) {
  z <- strsplit(sub("^S3method\\((.*)\\)$", "\\1", x), ",", fixed = TRUE)[[1L]]
  generic <- gsub('^"|"$', "", z[[1L]])
  row(generic = generic, class = z[[2L]], method = paste0(generic, ".", z[[2L]]), registration = x)
}))
reexport_names <- intersect(exports, imports$symbol)
category <- function(name) {
  if (name %in% reexport_names) return("reexport")
  if (name == "ggchord") return("constructor")
  for (prefix in c("geom", "stat", "position", "coord", "scale", "guide", "theme", "read", "as"))
    if (startsWith(name, paste0(prefix, "_"))) return(prefix)
  if (name %in% s3$method) return("S3")
  "utility_or_internal"
}
# Read saved Rd as documentation evidence. Do not evaluate examples or macros.
rd_docs <- list()
rd_text <- function(x) gsub("[[:space:]]+", " ", trimws(paste(unlist(x, use.names = FALSE), collapse = "")))
for (file in source_rd_files) {
  rd <- tools::parse_Rd(file)
  tags <- vapply(rd, function(x) { y <- attr(x, "Rd_tag"); if (is.null(y)) "" else y }, "")
  aliases <- vapply(rd[tags == "\\alias"], rd_text, "")
  parameters <- list()
  for (argblock in rd[tags == "\\arguments"]) for (item in argblock)
    if (identical(attr(item, "Rd_tag"), "\\item")) {
      for (nm in trimws(strsplit(rd_text(item[[1L]]), ",", fixed = TRUE)[[1L]]))
        parameters[[nm]] <- rd_text(item[[2L]])
    }
  value <- if (any(tags == "\\value")) rd_text(rd[tags == "\\value"][[1L]]) else "NOT_DOCUMENTED"
  for (alias in aliases) rd_docs[[alias]] <- list(file = relpath(file), parameters = parameters, value = value)
}
doc_for <- function(name) {
  x <- rd_docs[[name]]
  if (is.null(x)) list(file = "NOT_DOCUMENTED", parameters = list(), value = "NOT_DOCUMENTED") else x
}
# Stores syntax only, not evaluated functions/objects.
fn_syntax <- list(); top_bindings <- list(); parse_counts <- list()
operators <- c("{", "(", "[", "[[", "$", "@", "<-", "=", "<<-", "->", "->>",
  "if", "for", "while", "repeat", "break", "next", "function", "::", ":::",
  "+", "-", "*", "/", "^", "%%", "%/%", "==", "!=", "<", ">", "<=", ">=", "!", "&", "&&", "|", "||", "~", ":")
dynamic_calls <- c("quote", "substitute", "bquote", "expression", "eval", "evalq", "parse",
  "assign", "get", "mget", "do.call", "match.fun", "getS3method", "registerS3method",
  "register_theme_elements", "ggproto", "body", "formals", "environment")
literal_chars <- function(x) {
  if (is.character(x)) return(x)
  if (is.call(x) && leaf_name(call_name(x)) == "c") {
    parts <- as.list(x)[-1L]
    if (all(vapply(parts, is.character, logical(1L)))) return(unlist(parts, use.names = FALSE))
  }
  character()
}
# Every FUNCTION token must be visited exactly once. Token coordinates refer to
# syntax start, unlike the original helper's mostly absent call srcref attributes.
for (file in source_r_files) {
  file_rel <- relpath(file)
  parsed <- parse(file, keep.source = TRUE) # fail rather than silently skipping a file
  pd <- getParseData(parsed)
  loc <- pd[pd$token == "FUNCTION", c("line1", "col1", "line2", "col2")]
  loc <- loc[order(loc$line1, loc$col1), , drop = FALSE]
  fn_index <- 0L
  syntax_calls <- 0L
  top_sr <- attr(parsed, "srcref")
  walk <- function(x, scope = "<top-level>", owner = "", line = 1L,
                   quoted = FALSE, label = "", assigned = "", top = FALSE) {
    if (is.symbol(x)) {
      add("references", function_id = scope, symbol = as.character(x), source_file = file_rel,
          quoted = quoted)
      return(invisible(NULL))
    }
    if (is.pairlist(x) || is.expression(x)) {
      if (length(x)) for (i in seq_along(x)) {
        text <- paste(deparse(x[[i]]), collapse = " ")
        if (nzchar(text)) walk(x[[i]], scope, owner, line, quoted)
      }
      return(invisible(NULL))
    }
    if (!is.call(x)) return(invisible(NULL))
    h <- call_name(x)
    leaf <- leaf_name(h)
    syntax_calls <<- syntax_calls + 1L
    if (h == "function") {
      fn_index <<- fn_index + 1L
      stopifnot(fn_index <= nrow(loc))
      start <- loc[fn_index, ]
      fnline <- start$line1
      fname <- if (nzchar(label)) label else paste0("<anonymous@", fnline, ":", start$col1, ">")
      id <- paste0(file_rel, ":", fnline, ":", start$col1)
      qualified <- if (scope == "<top-level>") fname else paste0(scope, "::", fname)
      # A function assigned to a top-level symbol is a direct named definition.
      is_top <- top && nzchar(assigned) && !grepl("[$@\\[]", assigned)
      kind <- if (is_top) "named_top_level" else if (nzchar(assigned)) "assigned_nested_or_method" else if (nzchar(label)) "named_callback_or_method" else "anonymous_callback"
      args <- as.list(x[[2L]])
      body <- x[[3L]]
      tail_expr <- if (is.call(body) && call_name(body) == "{") body[[length(body)]] else body
      add("functions", function_id = id, name = fname, qualified_name = qualified,
          source_file = file_rel, source_line = fnline, source_column = start$col1,
          scope = scope, owner = owner, definition_kind = kind, quoted_context = quoted,
          is_top_level = is_top, is_exported = is_top && assigned %in% exports,
          category = category(if (is_top) assigned else owner), parameter_count = length(args),
          tail_expression = fmt(tail_expr, 240L), evidence = "AST_FUNCTION_LITERAL",
          semantic_review = "REVIEW_REQUIRED")
      fn_syntax[[id]] <<- list(args = args, body = body, name = fname, source_file = file_rel,
                              source_line = fnline, is_top = is_top)
      # Missing defaults are language missing-argument objects: indexing the
      # pairlist through a one-element list keeps them safe for deparse().
      if (length(args)) for (j in seq_along(args)) {
        nm <- names(args)[[j]]
        default <- paste(deparse(args[[j]], width.cutoff = 400L), collapse = " ")
        required <- !nzchar(default) && nm != "..."
        add("arguments", function_id = id, `function` = fname, qualified_function = qualified,
            source_file = file_rel, source_line = fnline, is_exported = is_top && assigned %in% exports,
            category = category(if (is_top) assigned else owner), parameter_index = j,
            parameter = nm, default_expr = if (nm == "...") "<ellipsis>" else if (required) "<required>" else default,
            is_required_syntax = required, signature_evidence = "AST_FUNCTION_LITERAL")
      }
      # Traverse defaults before body to match source token order; anonymous
      # functions in default expressions must also be covered.
      if (length(x[[2L]])) for (j in seq_along(x[[2L]])) {
        # Missing defaults cannot be bound to a local variable.
        text <- paste(deparse(x[[2L]][[j]]), collapse = " ")
        if (nzchar(text)) walk(x[[2L]][[j]], id, owner, fnline, quoted)
      }
      walk(body, id, if (is_top) assigned else owner, fnline, quoted)
      return(invisible(NULL))
    }
    if (h %in% c("<-", "=", "<<-") && length(x) == 3L) {
      lhs <- if (is.symbol(x[[2L]])) as.character(x[[2L]]) else fmt(x[[2L]], 200L)
      rhs <- x[[3L]]
      rhs_call <- call_name(rhs)
      binding_kind <- if (is.call(rhs) && rhs_call == "function") "function_literal" else if (is.symbol(rhs)) "symbol_alias_or_binding" else if (leaf_name(rhs_call) == "ggproto") "ggproto_construction" else if (is.call(rhs)) "call_result" else "literal_or_expression"
      if (scope == "<top-level>") {
        add("objects", name = lhs, source_file = file_rel, source_line = line,
            binding_kind = binding_kind, rhs = fmt(rhs, 240L),
            rhs_symbol = if (is.symbol(rhs)) as.character(rhs) else "", factory = rhs_call,
            is_exported = lhs %in% exports, quoted_context = quoted,
            binding_context = if (quoted) "QUOTED_SYNTAX" else if (top) "TOP_LEVEL_EXPRESSION" else "ENCLOSED_AST_CANDIDATE",
            location_precision = "ENCLOSING_TOP_LEVEL_EXPRESSION")
        if (top && is.symbol(x[[2L]])) top_bindings[[lhs]] <<- list(rhs = rhs, file = file_rel, line = line)
      }
      if (grepl("(^|[^a-zA-Z])(body|formals|environment)\\(", lhs))
        add("dynamic", source_file = file_rel, source_line = line, function_id = scope,
            operation = "FUNCTION_MUTATION", expression = paste0(lhs, " <- ", fmt(rhs, 180L)),
            quoted = quoted, evidence = "AST_SYNTAX_CANDIDATE")
      if (grepl("[$]", lhs) && grepl("default_aes|geom_defaults", lhs)) {
        # names(geom$default_aes)[...] is a collection transformation, not an
        # aesthetic called 'default_aes)' or a deparsed subscript fragment.
        lhs_expr <- x[[2L]]
        direct_field <- is.call(lhs_expr) && call_name(lhs_expr) == "$" &&
          is.symbol(lhs_expr[[3L]]) && !as.character(lhs_expr[[3L]]) %in% c("default_aes", "geom_defaults")
        aes_owner <- if (nzchar(owner)) owner else if (is.call(lhs_expr) && call_name(lhs_expr) == "$")
          fmt(lhs_expr[[2L]], 160L) else "<top-level-transform>"
        add("aesthetics", source_file = file_rel, source_line = line, owner = aes_owner,
            declaration = if (direct_field) "default_aes_assignment" else "default_aes_dynamic_transform",
            aesthetic = if (direct_field) as.character(lhs_expr[[3L]]) else "<dynamic:default_aes_transform>",
            default_or_expression = paste0(lhs, " <- ", fmt(rhs, 180L)), review_status = "REVIEW_REQUIRED")
      }
      walk(x[[2L]], scope, owner, line, quoted)
      walk(rhs, scope, if (scope == "<top-level>") lhs else owner, line, quoted,
           label = lhs, assigned = lhs, top = top)
      return(invisible(NULL))
    }
    if (!h %in% operators)
      add("calls", caller = scope, target = h, source_file = file_rel,
          source_line = line, relation = if (quoted) "QUOTED_CALL_SYNTAX" else "CALL_SYNTAX_CANDIDATE")
    if (leaf %in% dynamic_calls)
      add("dynamic", source_file = file_rel, source_line = line, function_id = scope,
          operation = leaf, expression = fmt(x, 220L), quoted = quoted,
          evidence = "AST_SYNTAX_CANDIDATE")
    if (leaf == "ggchord_capture_layer_input" && length(x) >= 5L) {
      role_names <- literal_chars(x[[5L]])
      if (!length(role_names)) role_names <- "<dynamic:input_roles>"
      for (ae in role_names) add("aesthetics", source_file = file_rel,
        source_line = line, owner = owner, declaration = "input_role",
        aesthetic = ae, default_or_expression = fmt(x[[5L]], 220L), review_status = "REVIEW_REQUIRED")
    }
    names_x <- names(x)
    children <- if (length(x) >= 2L) seq.int(2L, length(x)) else integer()
    for (j in children) {
      argname <- if (length(names_x) >= j) names_x[[j]] else ""
      if (is.na(argname)) argname <- ""
      if (argname %in% c("default_aes", "required_aes", "optional_aes", "non_missing_aes", "aesthetics", "renames")) {
        val <- x[[j]]
        aes_names <- if (argname == "default_aes" && is.call(val) && leaf_name(call_name(val)) == "aes") names(as.list(val)[-1L]) else literal_chars(val)
        if (!length(aes_names)) aes_names <- "<dynamic:REVIEW_REQUIRED>"
        for (ae in aes_names) add("aesthetics", source_file = file_rel, source_line = line,
          owner = owner, declaration = argname, aesthetic = ae,
          default_or_expression = fmt(val, 220L), review_status = "REVIEW_REQUIRED")
      }
    }
    # Call heads can themselves contain functions (e.g. (function(...) ...)()).
    if (is.call(x[[1L]]) || (is.symbol(x[[1L]]) && !h %in% operators))
      walk(x[[1L]], scope, owner, line, quoted)
    next_quoted <- quoted || leaf %in% c("quote", "substitute", "bquote", "expression")
    for (j in children) {
      # Arguments like f(x = ) are missing syntax, not NULL expressions.
      txt <- paste(deparse(x[[j]]), collapse = " ")
      if (!nzchar(txt)) next
      argname <- if (length(names_x) >= j) names_x[[j]] else ""
      if (is.na(argname)) argname <- ""
      walk(x[[j]], scope, owner, line, next_quoted, label = argname)
    }
    invisible(NULL)
  }
  for (i in seq_along(parsed)) {
    start_line <- if (length(top_sr) >= i) as.integer(top_sr[[i]][1L]) else 1L
    walk(parsed[[i]], line = start_line, top = TRUE)
  }
  if (fn_index != nrow(loc)) stop("Function-token coverage mismatch in ", file_rel, ": ", fn_index, " / ", nrow(loc))
  parse_counts[[length(parse_counts) + 1L]] <- row(source_file = file_rel,
    top_level_expressions = length(parsed), function_tokens = nrow(loc), visited_functions = fn_index,
    call_nodes = syntax_calls, parse_status = "OK")
}
functions <- collapse_records("functions")
arguments <- collapse_records("arguments")
objects <- collapse_records("objects")
calls <- unique(collapse_records("calls"))
references <- collapse_records("references")
dynamic <- unique(collapse_records("dynamic"))
aesthetics <- unique(collapse_records("aesthetics"))
coverage <- do.call(rbind, parse_counts)
stopifnot(nrow(functions) == sum(coverage$function_tokens))
s3$function_id <- functions$function_id[match(s3$method, functions$name)]
s3$definition_source <- paste0(functions$source_file[match(s3$method, functions$name)], ":",
                               functions$source_line[match(s3$method, functions$name)])
if (anyNA(s3$function_id)) stop("Unlocated S3 method: ", paste(s3$method[is.na(s3$function_id)], collapse = ", "))
# Resolve direct definitions and aliases to their source signature. This is
# syntax resolution, not proof of dynamic runtime dispatch or stable contracts.
resolve <- function(name, seen = character()) {
  if (name %in% seen) return(list(id = "", mode = "CYCLIC_OR_UNRESOLVED", origin = name))
  k <- which(functions$is_top_level & functions$name == name)
  if (length(k)) return(list(id = functions$function_id[[k[[1L]]]], mode = "DIRECT_FUNCTION", origin = name))
  b <- top_bindings[[name]]
  if (is.null(b)) return(list(id = "", mode = "NO_LOCAL_SIGNATURE", origin = name))
  if (is.symbol(b$rhs)) {
    z <- resolve(as.character(b$rhs), c(seen, name))
    z$mode <- paste0("SYMBOL_ALIAS/", z$mode)
    return(z)
  }
  # Explicit source-backed factory case: clones theme_ggchord, changes body
  # and environment, returns out; no formals(out) assignment exists.
  if (is.call(b$rhs) && identical(call_name(b$rhs), "ggchord_theme_variant")) {
    factory <- top_bindings[["ggchord_theme_variant"]]
    text <- if (is.null(factory)) "" else fmt(factory$rhs)
    if (grepl("out <- theme_ggchord", text, fixed = TRUE) &&
        grepl("body(out) <-", text, fixed = TRUE) && !grepl("formals(out)", text, fixed = TRUE)) {
      z <- resolve("theme_ggchord", c(seen, name))
      z$mode <- "SOURCE_TRACED_FACTORY_CLONE_REQUIRES_RUNTIME_REVIEW"
      return(z)
    }
  }
  list(id = "", mode = "DYNAMIC_SIGNATURE_REVIEW_REQUIRED", origin = name)
}
export_rows <- list()
for (name in exports) {
  imported <- name %in% reexport_names
  z <- resolve(name)
  b <- top_bindings[[name]]
  source <- if (imported) paste0(imports$package[match(name, imports$symbol)], "::", name) else if (!is.null(b)) paste0(b$file, ":", b$line) else "REVIEW_REQUIRED"
  export_rows[[length(export_rows) + 1L]] <- row(symbol = name, category = category(name),
    is_third_party_reexport = imported, definition_source = source,
    signature_function_id = z$id, signature_origin = z$origin,
    signature_evidence = if (imported) "NAMESPACE_IMPORT_FROM" else z$mode,
    documentation_source = doc_for(name)$file, decision = "REVIEW_REQUIRED",
    confidence = "STATIC_BASELINE_ONLY", source_revision = head_revision)
  if (!imported && nzchar(z$id) && !identical(name, z$origin)) {
    template <- arguments[arguments$function_id == z$id, , drop = FALSE]
    template$function_id <- paste0("export:", name)
    template[["function"]] <- name
    template$qualified_function <- name
    template$is_exported <- TRUE
    template$category <- category(name)
    template$signature_evidence <- z$mode
    if (!is.null(b)) { template$source_file <- b$file; template$source_line <- b$line }
    arguments <- rbind(arguments, template)
  }
}
export_table <- do.call(rbind, export_rows)
# Function argument semantics are deliberately not guessed from parameter names.
# Documentation text and syntax facts are separate fields from reviewed contract.
semantic_cols <- c("type_contract", "allowed_values", "unit", "null_semantics", "na_semantics",
  "applies_to_coords", "applies_to_geoms", "validation_location", "forwarded_to", "duplicated_with",
  "deprecation_status", "proposed_name", "proposed_default", "decision", "rationale", "test_case", "migration_example")
for (col in semantic_cols) arguments[[col]] <- "REVIEW_REQUIRED"
arguments$behavior_source <- paste0(arguments$source_file, ":", arguments$source_line)
arguments$documentation_source <- "NOT_DOCUMENTED"
arguments$documented_parameter <- "NOT_DOCUMENTED"
arguments$actual_usage <- "REVIEW_REQUIRED; source occurrence is not runtime use"
arguments$missing_semantics <- ifelse(arguments$parameter == "...",
  "ELLIPSIS_COLLECTOR; zero arguments allowed; forwarding contract REVIEW_REQUIRED",
  ifelse(arguments$is_required_syntax, "NO_DEFAULT_IN_SIGNATURE; runtime requiredness REVIEW_REQUIRED",
         "HAS_DEFAULT; semantic meaning REVIEW_REQUIRED"))
arguments$auto_none_semantics <- "REVIEW_REQUIRED"
arguments$confidence <- "STATIC_SIGNATURE_ONLY"
arguments$source_revision <- head_revision
reference_counts <- table(paste(references$function_id, references$symbol, sep = "\t"))
for (i in seq_len(nrow(arguments))) {
  doc <- doc_for(arguments[["function"]][[i]])
  text <- doc$parameters[[arguments$parameter[[i]]]]
  arguments$documentation_source[[i]] <- doc$file
  if (!is.null(text)) arguments$documented_parameter[[i]] <- text
  fid <- arguments$function_id[[i]]
  if (startsWith(fid, "export:")) fid <- resolve(arguments[["function"]][[i]])$id
  n <- unname(reference_counts[match(paste(fid, arguments$parameter[[i]], sep = "\t"), names(reference_counts))])
  if (is.na(n)) n <- 0L
  arguments$actual_usage[[i]] <- paste0("AST body/default symbol occurrences=", n, "; runtime use REVIEW_REQUIRED")
}
# Apply curated evidence to unique public function+parameter rows, or private
# function_id+parameter rows. AST defaults, order and current source provenance
# are immutable; historical review provenance is recorded in separate columns.
arguments$static_confidence <- arguments$confidence
arguments$static_actual_usage <- arguments$actual_usage
arguments$review_status <- "REVIEW_REQUIRED"
arguments$review_evidence <- "NOT_REVIEWED"
arguments$review_input_path <- "NONE"
arguments$review_source_revision <- "NONE"
review_allowed <- c("function", "function_id", "parameter", semantic_cols, "actual_usage", "confidence",
                    "missing_semantics", "auto_none_semantics",
                    "review_status", "review_evidence", "source_revision")
review_required <- c("function", "parameter", "review_status", "review_evidence", "source_revision")
review_inputs <- lapply(review_files, function(file) {
  d <- read.csv(file, check.names = FALSE, na.strings = character(), colClasses = "character")
  if (anyDuplicated(names(d))) stop("Duplicate review column in ", relpath(file))
  unknown <- setdiff(names(d), review_allowed)
  missing <- setdiff(review_required, names(d))
  if (length(unknown) || length(missing))
    stop("Invalid review schema in ", relpath(file), "; unknown columns: ", paste(unknown, collapse = ", "),
         "; missing columns: ", paste(missing, collapse = ", "))
  if (nrow(d) && any(!nzchar(as.matrix(d[, review_required, drop = FALSE]))))
    stop("Empty required review provenance or key in ", relpath(file))
  if (nrow(d) && any(!grepl("^[0-9a-f]{40}$", d$source_revision)))
    stop("Review source_revision must be a full commit SHA in ", relpath(file))
  list(path = relpath(file), data = d)
})
review_keys <- unlist(lapply(review_inputs, function(x) {
  d <- x$data
  key <- if ("function_id" %in% names(d)) paste0("id:", d$function_id) else paste0("export:", d[["function"]])
  paste(key, d$parameter, sep = "\t")
}), use.names = FALSE)
if (anyDuplicated(review_keys))
  stop("Duplicate function+parameter review keys: ", paste(unique(review_keys[duplicated(review_keys)]), collapse = ", "))
for (input in review_inputs) {
  d <- input$data
  fields <- intersect(names(d), c(semantic_cols, "actual_usage", "confidence", "review_status", "review_evidence",
                                  "missing_semantics", "auto_none_semantics"))
  for (i in seq_len(nrow(d))) {
    target <- if ("function_id" %in% names(d)) {
      which(!arguments$is_exported & arguments$function_id == d$function_id[[i]] &
            arguments[["function"]] == d[["function"]][[i]] & arguments$parameter == d$parameter[[i]])
    } else {
      which(arguments$is_exported & arguments[["function"]] == d[["function"]][[i]] &
            arguments$parameter == d$parameter[[i]])
    }
    if (length(target) != 1L)
      stop("Review key must match exactly one parameter in its declared scope: ", d[["function"]][[i]], " / ",
           d$parameter[[i]], " in ", input$path, "; matches=", length(target))
    for (field in fields) arguments[[field]][[target]] <- d[[field]][[i]]
    arguments$review_input_path[[target]] <- input$path
    arguments$review_source_revision[[target]] <- d$source_revision[[i]]
  }
}
reviewed <- arguments$review_status %in% c("SOURCE_REVIEWED", "RUNTIME_REVIEWED")
reviewed_argument_rows <- sum(reviewed)
reviewed_export_names <- length(unique(arguments[["function"]][reviewed & arguments$is_exported]))
# A local identifier in a call position is only a dependency candidate: lexical
# shadowing, quote/eval, ggproto inheritance and dynamic registration remain open.
local_names <- unique(c(functions$name[functions$is_top_level], objects$name, exports))
internal <- calls[leaf_name(calls$target) %in% local_names, , drop = FALSE]
internal$target_definition_candidates <- vapply(leaf_name(internal$target), function(nm) {
  paste(functions$function_id[functions$is_top_level & functions$name == nm], collapse = ";")
}, "")
internal$resolution <- "LEXICAL_CANDIDATE_NOT_RUNTIME_CALL_GRAPH"
return_rows <- lapply(exports, function(nm) {
  z <- resolve(nm); doc <- doc_for(nm)
  tail <- if (nzchar(z$id)) functions$tail_expression[match(z$id, functions$function_id)] else "EXTERNAL_OR_DYNAMIC"
  row(`function` = nm, source = export_table$definition_source[match(nm, exports)],
    documented_return = doc$value, documentation_source = doc$file,
    tail_expression = tail, stable_fields = "REVIEW_REQUIRED", source_tracking = "REVIEW_REQUIRED",
    mutation_and_cache_ownership = "REVIEW_REQUIRED", confidence = "DOCUMENTATION_AND_SYNTAX_ONLY")
})
compat <- do.call(rbind, lapply(exports[grepl("^(geom|stat)_", exports)], function(nm)
  row(layer = nm, coordinate = c("coord_chord", "coord_circular", "coord_collinear"),
      support_status = "REVIEW_REQUIRED", evidence = export_table$definition_source[match(nm, exports)],
      documentation_source = doc_for(nm)$file, runtime_validation = "NOT_RUN")))
reexports <- imports[imports$symbol %in% reexport_names, , drop = FALSE]
reexports$retention_decision <- "REVIEW_REQUIRED"
reexports$replacement_call <- paste0(reexports$package, "::", reexports$symbol)
reexports$runtime_signature_scope <- "EXTERNAL_PACKAGE_NOT_EVALUATED"
# Reviewed entry contracts use their own key space, not parameter overlays.
# Keep current syntax/documentation and historical review provenance separate.
read_detail <- function(name, keys, required) {
  path <- file.path(repo, "design-roadmap/audit-reviewed-details", paste0(name, ".csv"))
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path, check.names = FALSE, na.strings = character(), colClasses = "character")
  needed <- unique(c(keys, required, "review_status", "source_revision"))
  if (anyDuplicated(names(d)) || length(setdiff(needed, names(d)))) stop("Invalid detail schema: ", path)
  if (nrow(d) && (any(!nzchar(as.matrix(d[, needed, drop = FALSE]))) ||
                 any(!grepl("^[0-9a-f]{40}$", d$source_revision)))) stop("Empty detail or invalid provenance: ", path)
  key <- do.call(paste, c(d[keys], sep = "\t"))
  if (anyDuplicated(key)) stop("Duplicate detail key: ", path)
  attr(d, "input_path") <- relpath(path)
  d
}
import_reviews <- read_detail("import-contracts", c("package", "symbol"),
  c("import_role", "decision", "rationale", "evidence", "unknowns_next_phase"))
if (!is.null(import_reviews)) {
  key <- function(x) paste(x$package, x$symbol, sep = "\t")
  if (!setequal(key(imports), key(import_reviews))) stop("Import review must cover exactly all declarations")
  d <- import_reviews[match(key(imports), key(import_reviews)), , drop = FALSE]
  for (field in c("import_role", "decision", "rationale", "evidence", "unknowns_next_phase", "review_status"))
    imports[[field]] <- d[[field]]
  imports$review_source_revision <- d$source_revision
}
entry_reviews <- read_detail("export-decisions", "function", c("decision", "rationale", "return_contract",
  "source_tracking", "mutation_and_cache_ownership", "coordinate_scope", "evidence", "unknowns_next_phase"))
returns <- do.call(rbind, return_rows)
if (!is.null(entry_reviews)) {
  if (!setequal(entry_reviews[["function"]], exports)) stop("Entry decisions must cover exactly all exports")
  d <- entry_reviews[match(exports, entry_reviews[["function"]]), , drop = FALSE]
  for (field in c("decision", "rationale", "coordinate_scope", "unknowns_next_phase")) export_table[[field]] <- d[[field]]
  export_table$confidence <- d$review_status
  for (x in c("export_table", "returns")) {
    value <- get(x)
    value$review_status <- d$review_status
    value$review_evidence <- d$evidence
    value$review_source_revision <- d$source_revision
    value$review_input_path <- attr(entry_reviews, "input_path")
    assign(x, value)
  }
  returns$stable_fields <- d$return_contract
  returns$source_tracking <- d$source_tracking
  returns$mutation_and_cache_ownership <- d$mutation_and_cache_ownership
  returns$unknowns_next_phase <- d$unknowns_next_phase
  returns$confidence <- d$review_status
  ix <- match(reexports$symbol, d[["function"]])
  reexports$retention_decision <- d$decision[ix]
  reexports$rationale <- d$rationale[ix]
  reexports$review_source_revision <- d$source_revision[ix]
}
coordinate_reviews <- read_detail("coordinate-compatibility", c("layer", "coordinate"),
  c("support_status", "evidence", "runtime_validation"))
if (!is.null(coordinate_reviews)) {
  key <- function(x) paste(x$layer, x$coordinate, sep = "\t")
  if (!setequal(key(compat), key(coordinate_reviews))) stop("Coordinate review must cover exactly the layer/coordinate matrix")
  d <- coordinate_reviews[match(key(compat), key(coordinate_reviews)), , drop = FALSE]
  for (field in c("support_status", "evidence", "runtime_validation", "review_status")) compat[[field]] <- d[[field]]
  compat$review_source_revision <- d$source_revision
  compat$review_input_path <- attr(coordinate_reviews, "input_path")
}
function_reviews <- read_detail("internal-functions", "function_id", c("role", "return_contract",
  "state_ownership", "dynamic_dependency", "decision", "rationale", "unknowns_next_phase", "evidence"))
if (!is.null(function_reviews)) {
  ix <- match(function_reviews$function_id, functions$function_id)
  if (anyNA(ix)) stop("Unknown function_id in internal function review")
  if (any(functions$is_exported[ix])) stop("Internal function review must not overwrite public entry evidence")
  fields <- c("role", "return_contract", "state_ownership", "dynamic_dependency", "decision",
    "rationale", "unknowns_next_phase")
  for (field in fields) {
    functions[[field]] <- "REVIEW_REQUIRED"
    functions[[field]][ix] <- function_reviews[[field]]
  }
  functions$semantic_review[ix] <- function_reviews$review_status
  functions$review_source_revision <- "NONE"
  functions$review_source_revision[ix] <- function_reviews$source_revision
  functions$review_evidence <- "NOT_REVIEWED"
  functions$review_evidence[ix] <- function_reviews$evidence
  functions$review_input_path <- "NONE"
  functions$review_input_path[ix] <- attr(function_reviews, "input_path")
  hit <- match(s3$function_id, function_reviews$function_id)
  s3$review_status <- ifelse(is.na(hit), "REVIEW_REQUIRED", function_reviews$review_status[hit])
}
if (!is.null(entry_reviews)) {
  # Only direct exported function definitions have a matching syntax row.
  # Alias/factory entry views retain their separate entry provenance.
  ix <- which(functions$is_exported)
  hit <- match(functions$name[ix], entry_reviews[["function"]])
  if (anyNA(hit)) stop("Exported function syntax lacks an entry review")
  for (field in c("return_contract", "decision", "rationale", "unknowns_next_phase")) {
    if (is.null(functions[[field]])) functions[[field]] <- "REVIEW_REQUIRED"
    functions[[field]][ix] <- entry_reviews[[field]][hit]
  }
  for (field in c("review_source_revision", "review_evidence", "review_input_path"))
    if (is.null(functions[[field]])) functions[[field]] <- "NONE"
  functions$semantic_review[ix] <- entry_reviews$review_status[hit]
  functions$role[ix] <- paste0("公开入口 ", functions$name[ix], "；职责见返回合同和逐参数 actual_usage")
  functions$state_ownership[ix] <- entry_reviews$mutation_and_cache_ownership[hit]
  functions$dynamic_dependency[ix] <- paste0("公开调用边界；别名/工厂、NSE与延迟求值见逐参数合同；", entry_reviews$coordinate_scope[hit])
  functions$review_source_revision[ix] <- entry_reviews$source_revision[hit]
  functions$review_evidence[ix] <- entry_reviews$evidence[hit]
  functions$review_input_path[ix] <- attr(entry_reviews, "input_path")
}
stage_reviews <- read_detail("layout-stages", "stage", c("source_file", "reads", "writes",
  "upstream", "state_owner", "risk", "evidence"))
if (!is.null(stage_reviews)) {
  for (field in c("review_status", "stage_contract", "review_source_revision")) {
    objects[[field]] <- if (field == "review_status") "REVIEW_REQUIRED" else "NONE"
    dynamic[[field]] <- if (field == "review_status") "REVIEW_REQUIRED" else "NONE"
  }
  for (i in seq_len(nrow(stage_reviews))) {
    d <- stage_reviews[i, , drop = FALSE]
    hit <- which(objects$name == d$stage & objects$source_file == d$source_file &
                 objects$binding_context == "TOP_LEVEL_EXPRESSION")
    if (length(hit) != 1L) stop("Layout stage must match a unique top-level binding: ", d$stage)
    objects$review_status[hit] <- d$review_status
    objects$stage_contract[hit] <- paste0("audit-reviewed-details/layout-stages.csv#", d$stage)
    objects$review_source_revision[hit] <- d$source_revision
    hit <- which(dynamic$source_file == d$source_file & dynamic$operation == "quote")
    if (length(hit) != 1L) stop("Layout stage must match a unique quote: ", d$stage)
    dynamic$review_status[hit] <- d$review_status
    dynamic$stage_contract[hit] <- paste0("audit-reviewed-details/layout-stages.csv#", d$stage)
    dynamic$review_source_revision[hit] <- d$source_revision
    # eval calls are keyed by the exact called stage, not by nearby line.
    hit <- which(dynamic$source_file == "R/layout.R" & dynamic$operation == "eval" &
                 startsWith(dynamic$expression, paste0("eval(", d$stage, ",")))
    if (length(hit) != 1L) stop("Layout stage must match a unique eval: ", d$stage)
    dynamic$review_status[hit] <- d$review_status
    dynamic$stage_contract[hit] <- paste0("audit-reviewed-details/layout-stages.csv#", d$stage)
    dynamic$review_source_revision[hit] <- d$source_revision
  }
}
aes_keys <- c("source_file", "owner", "declaration", "aesthetic", "default_or_expression")
aes_fields <- c("semantic_role", "type_contract", "unit", "evaluation_stage", "null_na_semantics",
  "scale_or_geometry_owner", "decision", "rationale", "unknowns_next_phase")
aes_reviews <- read_detail("aesthetic-contracts", aes_keys, c(aes_fields, "evidence"))
if (!is.null(aes_reviews)) {
  key <- function(x) do.call(paste, c(x[aes_keys], sep = "\t"))
  if (!setequal(key(aesthetics), key(aes_reviews))) stop("Aesthetic review must match exactly current syntax records")
  d <- aes_reviews[match(key(aesthetics), key(aes_reviews)), , drop = FALSE]
  for (field in c(aes_fields, "review_status")) aesthetics[[field]] <- d[[field]]
  aesthetics$review_evidence <- d$evidence
  aesthetics$review_source_revision <- d$source_revision
  aesthetics$review_input_path <- attr(aes_reviews, "input_path")
}
# Attach reviewed ownership to binding/call views without converting syntax
# candidates into a runtime call graph or a deletion decision.
if (is.null(objects$review_status)) objects$review_status <- "REVIEW_REQUIRED"
if (is.null(objects$review_source_revision)) objects$review_source_revision <- "NONE"
objects$binding_semantics <- "AST_BINDING_ONLY; semantic ownership not reviewed yet"
objects$review_evidence <- "NOT_REVIEWED"
objects$unknowns_next_phase <- "本版仍须核查所属模块与间接使用；不能依据零引用删除"
if (!is.null(stage_reviews)) for (i in seq_len(nrow(stage_reviews))) {
  d <- stage_reviews[i, , drop = FALSE]
  hit <- which(objects$source_file == d$source_file & objects$binding_context == "QUOTED_SYNTAX")
  objects$review_status[hit] <- "DELAYED_STAGE_SCOPE_REVIEWED"
  objects$binding_semantics[hit] <- paste0("延迟表达式中的赋值语法，eval后属于当次共享frame；阶段：", d$stage,
    "；不是独立namespace对象或已冻结返回字段")
  objects$review_evidence[hit] <- d$evidence
  objects$review_source_revision[hit] <- d$source_revision
  objects$unknowns_next_phase[hit] <- paste0("v0.16.3拆阶段时核查具体临时值的生产/消费与失效；", d$risk)
  hit <- which(objects$source_file == d$source_file & objects$name == d$stage &
               objects$binding_context == "TOP_LEVEL_EXPRESSION")
  objects$binding_semantics[hit] <- paste0("quote表达式绑定；由compute_chord_layout eval到当次共享frame；阶段：", d$stage)
  objects$review_evidence[hit] <- d$evidence
  objects$unknowns_next_phase[hit] <- paste0("v0.16.3阶段拆分前保持当前读写/顺序可比；", d$risk)
}
for (i in which(functions$is_top_level & functions$semantic_review %in% c("SOURCE_REVIEWED", "RUNTIME_REVIEWED"))) {
  f <- functions[i, , drop = FALSE]
  hit <- which(objects$name == f$name & objects$source_file == f$source_file &
               objects$binding_context == "TOP_LEVEL_EXPRESSION" & objects$binding_kind == "function_literal")
  if (length(hit) != 1L) stop("Reviewed top-level function must have one binding: ", f$function_id)
  objects$review_status[hit] <- "FUNCTION_CONTRACT_REVIEWED"
  objects$binding_semantics[hit] <- paste0("函数绑定；契约见api_functions.csv#", f$function_id)
  objects$review_evidence[hit] <- f$review_evidence
  objects$review_source_revision[hit] <- f$review_source_revision
  objects$unknowns_next_phase[hit] <- f$unknowns_next_phase
}
if (!is.null(entry_reviews)) for (i in seq_len(nrow(entry_reviews))) {
  d <- entry_reviews[i, , drop = FALSE]
  hit <- which(objects$name == d[["function"]] & objects$binding_context == "TOP_LEVEL_EXPRESSION" &
               objects$binding_kind != "function_literal")
  if (!length(hit)) next
  if (length(hit) != 1L) stop("Export alias/factory must have one top-level binding: ", d[["function"]])
  objects$review_status[hit] <- "ENTRY_ALIAS_OR_FACTORY_REVIEWED"
  objects$binding_semantics[hit] <- d$return_contract
  objects$review_evidence[hit] <- d$evidence
  objects$review_source_revision[hit] <- d$source_revision
  objects$unknowns_next_phase[hit] <- d$unknowns_next_phase
}
object_reviews <- read_detail("object-contracts", c("source_file", "name", "binding_context"),
  c("binding_semantics", "decision", "unknowns_next_phase", "evidence"))
if (!is.null(object_reviews)) {
  binding_key <- function(x) paste(x$source_file, x$name, x$binding_context, sep = "\t")
  ix <- match(binding_key(object_reviews), binding_key(objects))
  if (anyNA(ix) || anyDuplicated(binding_key(objects)[ix])) stop("Unknown or ambiguous reviewed object binding")
  if (any(objects$review_status[ix] != "REVIEW_REQUIRED")) stop("Object detail must not overwrite another review scope")
  for (field in c("binding_semantics", "unknowns_next_phase", "review_status"))
    objects[[field]][ix] <- object_reviews[[field]]
  objects$review_evidence[ix] <- object_reviews$evidence
  objects$review_source_revision[ix] <- object_reviews$source_revision
}
for (name in c("internal", "dynamic")) {
  value <- get(name)
  key <- if (name == "internal") value$caller else value$function_id
  ix <- match(key, functions$function_id)
  value$caller_review_status <- ifelse(is.na(ix), "NO_FUNCTION_OWNER", functions$semantic_review[ix])
  value$caller_contract <- ifelse(is.na(ix), "TOP_LEVEL_OR_DYNAMIC_CONTEXT",
    paste0("api_functions.csv#", functions$function_id[ix]))
  value$runtime_scope <- "SYNTAX_ONLY; caller review does not certify dynamic exhaustiveness or unusedness"
  assign(name, value)
}
# All call candidates keep syntax-only scope. Attach the already reviewed
# enclosing contract, never infer execution, purity or unusedness from it.
ix <- which(dynamic$review_status == "REVIEW_REQUIRED" &
  dynamic$caller_review_status %in% c("SOURCE_REVIEWED", "RUNTIME_REVIEWED"))
hit <- match(dynamic$function_id[ix], functions$function_id)
dynamic$review_status[ix] <- "CALLER_SOURCE_SCOPE_REVIEWED"
dynamic$stage_contract[ix] <- dynamic$caller_contract[ix]
dynamic$review_source_revision[ix] <- functions$review_source_revision[hit]
if (!is.null(stage_reviews)) for (i in seq_len(nrow(stage_reviews))) {
  d <- stage_reviews[i, , drop = FALSE]
  ix <- which(dynamic$review_status == "REVIEW_REQUIRED" & dynamic$quoted &
    dynamic$source_file == d$source_file & dynamic$function_id == "<top-level>")
  dynamic$review_status[ix] <- "DELAYED_STAGE_SCOPE_REVIEWED"
  dynamic$stage_contract[ix] <- paste0("audit-reviewed-details/layout-stages.csv#", d$stage)
  dynamic$review_source_revision[ix] <- d$source_revision
}
if (!is.null(object_reviews)) for (i in seq_len(nrow(object_reviews))) {
  d <- object_reviews[i, , drop = FALSE]
  at <- which(objects$source_file == d$source_file & objects$name == d$name &
    objects$binding_context == d$binding_context)
  ix <- which(dynamic$review_status == "REVIEW_REQUIRED" & dynamic$function_id == "<top-level>" &
    dynamic$source_file == d$source_file & dynamic$source_line == objects$source_line[at])
  dynamic$review_status[ix] <- "OBJECT_CONSTRUCTION_SCOPE_REVIEWED"
  dynamic$stage_contract[ix] <- paste0("audit-reviewed-details/object-contracts.csv#", d$name)
  dynamic$review_source_revision[ix] <- d$source_revision
}
source_review_complete <- reviewed_argument_rows == nrow(arguments) &&
  all(functions$semantic_review %in% c("SOURCE_REVIEWED", "RUNTIME_REVIEWED")) &&
  all(objects$review_status != "REVIEW_REQUIRED") &&
  all(dynamic$review_status != "REVIEW_REQUIRED") &&
  !is.null(entry_reviews) && !is.null(coordinate_reviews) && !is.null(aes_reviews) && !is.null(import_reviews)
ref_summary <- aggregate(row(function_id = references$function_id),
                         list(function_id = references$function_id, symbol = references$symbol,
                              quoted = references$quoted), length)
names(ref_summary)[[4L]] <- "occurrences"
verify_snapshot()
hashes <- row(path = vapply(source_files, relpath, ""), md5 = unname(source_hashes_before))
tmp <- tempfile(); writeLines(paste(hashes$path, hashes$md5, sep = "\t"), tmp, useBytes = TRUE)
source_fingerprint <- unname(tools::md5sum(tmp)); unlink(tmp)
write_table(export_table, "api_exports"); write_table(imports, "api_imports"); write_table(s3, "api_s3")
write_table(functions, "api_functions"); write_table(arguments, "api_arguments")
write_table(objects, "api_objects"); write_table(aesthetics, "api_aesthetics")
write_table(returns, "api_return_contracts")
write_table(internal, "api_internal_dependencies"); write_table(reexports, "api_reexports")
write_table(compat, "api_compatibility_matrix"); write_table(dynamic, "api_dynamic_registration")
write_table(ref_summary, "api_symbol_references_HEURISTIC")
write_table(coverage, "source_coverage"); write_table(hashes, "source_hashes")
source_dirty <- git("status", "--short", "--", "R", "NAMESPACE", "DESCRIPTION", "man",
  shQuote("design-roadmap/audit-reviewed/*.csv"), shQuote("design-roadmap/audit-reviewed-details/*.csv"))
manifest <- row(key = c("audit_start_revision", "git_head_at_generation", "source_fingerprint_md5",
  "source_worktree_status", "source_snapshot_status", "generator_md5", "package_version", "R_version",
  "source_file_count", "export_names", "own_export_names", "reexports", "S3_registrations",
  "function_syntax_tokens", "function_rows", "anonymous_or_method_rows", "argument_rows_including_export_signature_views",
  "unresolved_own_export_signatures", "review_input_files", "human_reviewed_argument_rows", "human_reviewed_export_names",
  "entry_decision_rows", "coordinate_review_rows", "internal_function_review_rows", "layout_stage_review_rows",
  "aesthetic_syntax_rows", "aesthetic_review_rows",
  "bindings_with_reviewed_contract_or_scope", "S3_with_reviewed_function_contract",
  "semantic_review_status", "package_execution", "rebuild_command"),
 value = c(audit_start, head_revision, source_fingerprint,
   if (length(source_dirty)) paste(source_dirty, collapse = "; ") else "CLEAN_SOURCE_INPUTS",
   "VERIFIED_UNCHANGED_BEFORE_WRITE", generator_hash_before, package_version_before, R.version.string,
   nrow(coverage), length(exports), length(setdiff(exports, reexport_names)), length(reexport_names), nrow(s3),
   sum(coverage$function_tokens), nrow(functions), sum(!functions$is_top_level), nrow(arguments),
   sum(!export_table$is_third_party_reexport & !nzchar(export_table$signature_function_id)),
   length(review_files), reviewed_argument_rows, reviewed_export_names,
   if (is.null(entry_reviews)) 0L else nrow(entry_reviews),
   if (is.null(coordinate_reviews)) 0L else nrow(coordinate_reviews),
   if (is.null(function_reviews)) 0L else nrow(function_reviews),
   if (is.null(stage_reviews)) 0L else nrow(stage_reviews),
   nrow(aesthetics), if (is.null(aes_reviews)) 0L else nrow(aes_reviews),
   sum(objects$review_status != "REVIEW_REQUIRED"),
   if (is.null(s3$review_status)) 0L else sum(s3$review_status != "REVIEW_REQUIRED"),
   if (source_review_complete) "SOURCE_REVIEW_COVERED; runtime unknowns and deferred decisions remain explicit; NOT_RELEASE_ACCEPTANCE" else
     "OPEN; syntax coverage and partial semantic review are not whole-audit completion", "NOT_LOADED_OR_EXECUTED",
   paste("Rscript tools/audit-ggchord-api.R . design-roadmap/audit-generated", audit_start)))
write_table(manifest, "manifest")
readme <- c("# ggchord API 静态审计基线", "",
  "此目录由 `tools/audit-ggchord-api.R` 生成。完整覆盖状态见 manifest.csv；源码复核、局部运行验证与发布验收分别记录。", "",
  paste0("起始提交：`", audit_start, "`；生成时 HEAD：`", head_revision, "`。"),
  paste0("输入指纹（MD5，仅作可重复性标识）：`", source_fingerprint, "`。"),
  "源码、已保存 Rd 与人工复核输入的逐文件散列见 source_hashes.csv；工作树状态单独记录，不能把带改动源码冒充原始提交。",
  "解析前固定输入文件列表、散列与 HEAD，写表前再次核对；并发改动会使生成失败并要求重跑，不写入不一致的新产物。", "",
  "## 重建", "", "```sh", paste("Rscript tools/audit-ggchord-api.R . design-roadmap/audit-generated", audit_start), "```", "",
  "读取 CSV 时使用 `read.csv(path, na.strings = character(), check.names = FALSE)`，以免把代码默认表达式 NA 误读为 CSV 缺失值。", "",
  "## 覆盖", "", paste0("- ", nrow(coverage), " 个 R 文件全部成功解析；", sum(coverage$function_tokens), " 个 FUNCTION token 与访问数逐文件相等。"),
  paste0("- ", length(exports), " 个导出名称，", length(reexport_names), " 个重导出，", nrow(s3), " 个 S3 注册；导出不等于独立函数。"),
  paste0("- ", nrow(functions), " 个函数语法定义（含匿名、嵌套及 ggproto 方法）；", nrow(arguments), " 行参数含导出别名/预设的签名视图。"),
  "- api_objects.csv 补充顶层绑定、别名、ggproto 构造和工厂调用；api_dynamic_registration.csv 标记 quote/eval、注册、函数体/环境修改等。",
  "- 颜色拼写别名通过源绑定追踪；theme 预设的签名来自工厂复制源码，运行时形参另由108入口签名脚本核对。",
  "- api_arguments.csv 保存全部计划字段、默认表达式及 Rd 参数说明；未知语义显式 REVIEW_REQUIRED，文档说明不等于行为证据。",
  "- 默认表达式占位：<required> 表示无语法默认值，真实是否必填仍须核查 missing() 等分支；<ellipsis> 表示可为空的 ...，不标为必填。",
  "- api_aesthetics.csv 记录静态美学声明/映射或待查表达式；api_return_contracts.csv 保留文档返回说明及语法尾表达式，不推定稳定字段。",
  "- api_compatibility_matrix.csv 合并人工三坐标源码合同；runtime_validation 单独说明有限例子或 NOT_RUN，没有自动执行绘图。", "",
  "## 人工语义复核", "",
  paste0("- 本次纳入 ", length(review_files), " 份 audit-reviewed/*.csv，", reviewed_argument_rows,
         " 行参数、", reviewed_export_names, " 个公开入口标记为已做来源或运行复核；其余语义仍待复核。"),
  "- 参数表按 function+parameter 匹配公开参数；有 function_id 时必须精确匹配唯一内部函数及其参数。重复键、未知字段、签名字段或未定位键使生成失败。",
  "- 默认表达式、形参顺序、签名证据、当前 source_revision 不接受覆盖；原始 static_confidence/static_actual_usage 与语义复核分列。",
  "- review_source_revision、review_input_path、review_status、review_evidence 保留人工复核来源；复核提交与当前提交不同须核查相关源文件差异，不能自动把旧证据升级为当前行为。",
  "- audit-reviewed-details/ 下的 export-decisions 与 coordinate-compatibility 按不同键合并返回、去留和兼容记录；要求完整且无多余键，保留独立 review SHA。",
  "- internal-functions 按私有 function_id 合并角色、返回、状态与动态依赖证据；未覆盖的函数仍保留 REVIEW_REQUIRED，不将整文件解析冒充语义复核。",
  "- layout-stages 记录七个共享帧阶段的读写、顺序、状态和风险；精确挂接对应的 top-level quote 绑定及 eval 调用，其他动态记录仍待查。",
  "- aesthetic-contracts 按文件/owner/声明种类/美学或占位/表达式精确对应当前声明；input_role、draw aesthetic与动态default集合分开，源码复核不等于所有运行集合已穷举。",
  "- 对象表区分已复核函数绑定、公开别名/工厂和延迟阶段内的局部赋值语法；依赖/动态表只挂接caller合同，保持SYNTAX_ONLY，不将其变成删除或运行可达性结论。",
  "- object-contracts逐项记录非函数原型/常量/环境及local临时绑定；动态记录分别挂接函数、延迟阶段或原型构造作用域。状态表示源码拥有者已审，不表示动态结果已穷举。",
  "- 没有相应人工覆盖文件时照常生成，语义保持 REVIEW_REQUIRED；两个人工目录均纳入输入快照校验。", "",
  "## 局限与后续工作", "",
  "- 全部文件被解析不等于全部动态对象被穷举。运行生成的函数、动态注册、字符串回调、NSE、ggproto 继承及语言对象注入必须人工/定向运行复核。",
  "- 函数行号为 FUNCTION token 起点；顶层绑定和调用证据行号为所属表达式/函数的起点，location_precision 或证据类型说明精度。",
  "- 普通调用语法、被引用调用语法和符号出现次数均不是真实运行调用图。零引用既不证明无用，也不构成删除授权。",
  "- 声明型 default/required aesthetic 可提取；经运行修改/继承得来的最终集合仍待复核，重复声明不机械合并。",
  "- 不加载包，不解析或执行原始业务数据，不运行测试、包检查、vignette 或设备矩阵。语义、验证路径、NULL/NA、单位、迁移及去留决定仍需逐项证据。",
  "- 部分公开参数通过外部 ggplot2/grid 的 ... 透传，其完整外部签名不由本仓库静态函数列表重建。")
writeLines(readme, file.path(out, "README.md"), useBytes = TRUE)
cat(paste(manifest$key, manifest$value, sep = ": "), sep = "\n")
