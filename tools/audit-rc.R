# v0.16.6 current public contract coverage. Read-only introspection of the API;
# writes only rc-*-v0.16.6.csv, never the historical candidate snapshots.
# Run from the checkout root: Rscript tools/audit-rc.R
pkgload::load_all(quiet = TRUE)
out <- "design-roadmap/audit-generated"
read <- function(name) read.csv(file.path(out, paste0("candidate-", name,
  "-v0.16.4.csv")), na.strings = character(), check.names = FALSE)
write <- function(x, name) write.csv(x, file.path(out,
  paste0("rc-", name, "-v0.16.6.csv")), row.names = FALSE, na = "")
revision <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
ns <- asNamespace("ggchord")
exports <- sort(getNamespaceExports("ggchord"))
prior <- read("exports")
stopifnot(setequal(exports, prior$name))
sources <- list()
for (file in list.files("R", "[.]R$", full.names = TRUE)) {
  expr <- parse(file, keep.source = TRUE)
  refs <- attr(expr, "srcref")
  for (i in seq_along(expr)) {
    x <- expr[[i]]
    if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]]))
      sources[[as.character(x[[2L]])]] <- paste0(file, ":", refs[[i]][1L])
  }
}
rd_text <- function(x) gsub("[[:space:]]+", " ", trimws(paste(
  unlist(x, use.names = FALSE), collapse = "")))
docs <- list()
for (file in list.files("man", "[.]Rd$", full.names = TRUE)) {
  rd <- tools::parse_Rd(file)
  tags <- vapply(rd, function(x) attr(x, "Rd_tag") %||% "", character(1))
  args <- list()
  for (block in rd[tags == "\\arguments"]) for (item in block)
    if (identical(attr(item, "Rd_tag"), "\\item")) for (name in trimws(
      strsplit(rd_text(item[[1L]]), ",", fixed = TRUE)[[1L]]))
        args[[name]] <- rd_text(item[[2L]])
  for (name in vapply(rd[tags == "\\alias"], rd_text, character(1)))
    docs[[name]] <- list(file = file, args = args,
      value = if (any(tags == "\\value")) rd_text(rd[tags == "\\value"][[1L]])
        else "Upstream object; see its original package documentation")
}
default <- function(args, i) {
  if (names(args)[i] == "...") return("<ellipsis>")
  x <- paste(deparse(args[[i]], width.cutoff = 500), collapse = " ")
  if (nzchar(x)) x else "<required>"
}
prior$source <- vapply(prior$name, function(n) sources[[n]] %||%
  paste0(prior$ownership[match(n, prior$name)], " original namespace"), character(1))
prior$review <- "KEEP; current signatures/defaults and documented contract reviewed"
for (i in seq_len(nrow(prior))) {
  n <- prior$name[i]; value <- getExportedValue("ggchord", n)
  if (prior$ownership[i] != "ggchord") {
    stopifnot(identical(value, getExportedValue(prior$ownership[i], n)))
    prior$review[i] <- "KEEP_UPSTREAM_IDENTITY; upstream owns its evolving signature"
  } else {
    args <- formals(value)
    sig <- paste(paste0(names(args), "=", vapply(seq_along(args),
      function(j) default(args, j), character(1))), collapse = "; ")
    stopifnot(identical(sig, prior$signature[i]))
  }
  if (nzchar(prior$alias_of[i])) stopifnot(identical(value,
    getExportedValue("ggchord", prior$alias_of[i])))
}
write(prior[c("name", "family", "ownership", "alias_of", "source", "documentation",
  "review")], "exports")
parameters <- read("parameters")
for (i in seq_len(nrow(parameters))) {
  n <- parameters$function_name[i]; p <- parameters$parameter[i]
  args <- formals(getExportedValue("ggchord", n))
  stopifnot(p %in% names(args), identical(default(args, match(p, names(args))),
    parameters$default_expression[i]), nzchar(docs[[n]]$args[[p]]))
}
parameters$rd_changed <- vapply(seq_len(nrow(parameters)), function(i)
  !identical(parameters$current_documentation[i],
    docs[[parameters$function_name[i]]]$args[[parameters$parameter[i]]]), logical(1))
parameters$current_rd_text <- vapply(seq_len(nrow(parameters)), function(i)
  docs[[parameters$function_name[i]]]$args[[parameters$parameter[i]]], character(1))
parameters$source <- vapply(parameters$function_name, function(n) sources[[n]], character(1))
parameters$review <- "STATIC_CONTRACT_REVIEW; semantic baseline + current Rd, not individual runtime PASS"
parameters$baseline <- "candidate-parameters-v0.16.4.csv (function_name,parameter); reviewed semantics retained"
write(parameters[c("function_name", "parameter", "default_expression", "source", "current_rd",
  "rd_changed", "current_rd_text", "semantic_source", "current_semantic_override", "baseline", "review")], "parameters")
returns <- read("returns")
returns$current_return <- vapply(returns$function_name, function(n)
  docs[[n]]$value %||% "Upstream original object", character(1))
returns$rd_changed <- returns$current_return != returns$documented_return
returns$review <- "KEEP_DOCUMENTED_MEMBERS; current Rd plus detailed data/engine contracts; finite runtime only"
write(returns[c("function_name", "documentation", "detailed_contract", "rd_changed",
  "current_return", "review")], "returns")
s3 <- read("s3")
stopifnot(setequal(s3$registration, grep("^S3method", readLines("NAMESPACE"), value = TRUE)))
old_s3 <- read.csv(file.path(out, "api_s3.csv"))
s3$source <- vapply(s3$registration, function(r)
  sources[[old_s3$method[match(r, old_s3$registration)]]], character(1))
s3$review <- "KEEP_DISPATCH; ggproto callbacks remain internal, original generic protocol retained"
write(s3[c("registration", "source", "policy", "review")], "s3")
aesthetics <- read("aesthetics")
all_source <- paste(vapply(list.files("R", "[.]R$", full.names = TRUE),
  function(f) paste(readLines(f), collapse = "\n"), character(1)), collapse = "\n")
aesthetics$literal_in_current_source <- vapply(aesthetics$aesthetic_or_role,
  grepl, logical(1), x = all_source, fixed = TRUE)
stopifnot(all(aesthetics$literal_in_current_source))
aesthetics$review <- "KEEP_ROLE; R/zzz.R normalization + family geoms/scales; literal presence is not execution proof"
write(aesthetics, "aesthetics")
migrations <- read("migrations")
stopifnot(!anyDuplicated(paste(migrations$function_name, migrations$old_entry,
  migrations$entry_kind)))
migrations$review <- "RETAIN_MIGRATION; historical upper bounds distinct from first removal; representative guards only"
write(migrations[c("function_name", "old_entry", "entry_kind", "current_status", "replacement",
  "removed_by_version", "first_removed_version", "rejection_change_version", "guard", "decision", "review")], "migrations")
compat <- read("compatibility")
compat$rc_evidence <- "NOT_RUN_IN_RC; inherited source contract / finite evidence retained"
write(compat, "compatibility")
inputs <- sort(unique(c(list.files("R", "[.]R$", full.names = TRUE),
  list.files("man", "[.]Rd$", full.names = TRUE), "DESCRIPTION", "NAMESPACE", "NEWS.md",
  "tools/audit-rc.R", list.files(out, "^candidate-.*[.]csv$", full.names = TRUE))))
write(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs)),
  revision = revision, state = "HEAD plus file hashes; working bytes, not behavior evidence"), "source-hashes")
cat(length(exports), "exports;", nrow(parameters), "parameters;", nrow(aesthetics),
  "roles;", nrow(s3), "S3;", nrow(migrations), "migrations.\n")
cat(sum(parameters$rd_changed), "parameter prose differences;", sum(returns$rd_changed),
  "return prose differences from candidate; review rc tables.\n")
