# Focused audit check: introspect exported formals without calling the APIs.
# Run after the static inventory, from the repository root. No package tests,
# build/check, graphics output, vignette rebuilding or dependency installation.
if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload is required")
generated <- "design-roadmap/audit-generated"
out <- "design-roadmap/audit-validation"
read_inventory <- function(file) read.csv(file.path(generated, file),
  check.names = FALSE, na.strings = character())
exports <- read_inventory("api_exports.csv")
arguments <- read_inventory("api_arguments.csv")
hashes <- read_inventory("source_hashes.csv")
if (!identical(unname(tools::md5sum(hashes$path)), hashes$md5))
  stop("Source inputs differ from the static inventory; regenerate it first")
manifest <- read_inventory("manifest.csv")
revision <- manifest$value[manifest$key == "git_head_at_generation"]
fingerprint <- manifest$value[manifest$key == "source_fingerprint_md5"]
pkgload::load_all(quiet = TRUE)

format_default <- function(formals, i, parameter) {
  if (parameter == "...") return("<ellipsis>")
  text <- tryCatch(paste(deparse(formals[[i]], width.cutoff = 400L), collapse = " "),
    error = function(e) "")
  if (!nzchar(trimws(text))) "<required>" else text
}
results <- lapply(exports$symbol[!exports$is_third_party_reexport], function(name) {
  value <- getExportedValue("ggchord", name)
  expected <- arguments[arguments$is_exported & arguments[["function"]] == name, ]
  expected <- expected[order(expected$parameter_index), ]
  if (!is.function(value)) stop("Expected a local function export: ", name)
  actual <- formals(value)
  defaults <- vapply(seq_along(actual), function(i)
    format_default(actual, i, names(actual)[[i]]), character(1))
  expected_defaults <- expected$default_expr
  expected_defaults[expected$parameter == "..."] <- "<ellipsis>"
  names_match <- identical(unname(names(actual)), unname(expected$parameter))
  defaults_match <- identical(unname(defaults), unname(expected_defaults))
  data.frame(symbol = name, parameter_count = length(actual),
    names_and_order_match = names_match, defaults_match = defaults_match,
    status = if (names_match && defaults_match) "PASS_SIGNATURE_ONLY" else "MISMATCH",
    source_revision = revision, source_fingerprint_md5 = fingerprint,
    R_version = R.version.string,
    ggplot2_version = as.character(utils::packageVersion("ggplot2")),
    stringsAsFactors = FALSE)
})
results <- do.call(rbind, results)
if (!identical(unname(tools::md5sum(hashes$path)), hashes$md5))
  stop("Source inputs changed during signature introspection; rerun after regeneration")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(out, "api_export_signatures.csv"), row.names = FALSE)
writeLines(c("# Public signature introspection", "",
  "Rebuild with `Rscript tools/validate-api-signatures.R` after the static inventory.",
  "This loads the local development namespace and compares local exported formals only.",
  "The APIs themselves are not called. No data, geometry, parameter meaning, validation,",
  "return value, side effect or coordinate compatibility is certified by this check.",
  "Upstream ggplot2/grid reexports remain outside this repository's signature check.",
  "Source revision, input fingerprint and runtime versions are stored on each row."),
  file.path(out, "README.md"))
if (any(results$status != "PASS_SIGNATURE_ONLY")) {
  print(results[results$status != "PASS_SIGNATURE_ONLY", c("symbol", "names_and_order_match", "defaults_match")])
  stop("Export signature mismatch; do not accept the static coverage yet")
}
cat(nrow(results), "local exported signatures match names, order and defaults.\n")
