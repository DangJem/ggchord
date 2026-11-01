# Environment observations and deliberately untested combinations, no installs.
pkgload::load_all(quiet = TRUE)
row <- function(component, declaration, actual, status, scope) data.frame(
  component = component, declared_support = declaration, observed = actual,
  evidence_status = status, scope = scope)
rows <- list(
  row("R", ">=4.1.0", R.version.string, "LOCAL_FINITE", "rc-behavior-v0.16.6.csv; lower bound NOT_RUN"),
  row("ggplot2", ">=4.0.0", as.character(utils::packageVersion("ggplot2")), "LOCAL_FINITE", "candidate identity/signature and six behavior cases; lower bound NOT_RUN"),
  row("OS", "No OS exclusion in DESCRIPTION; ggplot2/grid and readers impose requirements",
    paste(Sys.info()[["sysname"]], Sys.info()[["release"]], R.version$platform), "LOCAL_FINITE", "This checkout only; no other OS executed"),
  row("PNG", "grDevices PNG via ggsave / target render", "9x8in 120dpi sans; 6x5 and 9x7 target measurements", "LOCAL_FINITE", "two visual scenes; two diagnostics; one preview; not all fonts/devices"),
  row("PDF", "R graphics device", "private pdf(NULL) 8x8in Helvetica metrics only", "METRICS_ONLY", "No v0.16.6 exported PDF acceptance or font embedding check"),
  row("SVG/Cairo/Quartz", "Device-specific optional support", "not invoked for RC export", "NOT_RUN", "v0.16.5 device facts are historical; no current PDF/SVG matrix"),
  row("Windows/Linux/other macOS", "Same declared R/dependency requirements", "no v0.16.6 run", "NOT_RUN", "Four-platform v0.16.5 CI recorded in RELEASE-ACCEPTANCE-v0.16.5.md, not inherited PASS"),
  row("R4.1 / ggplot2 4.0.0", "Declared lower bounds", "not installed/exercised in this RC", "NOT_RUN", "Needs separately authorized release compatibility evidence"))
for (p in c("rlang", "scales", "readxl", "digest", "GenomicRanges", "InteractionSet",
  "IRanges", "S4Vectors", "jsonlite", "Biostrings", "svglite", "qpdf", "knitr", "rmarkdown")) {
  found <- requireNamespace(p, quietly = TRUE)
  rows[[length(rows) + 1L]] <- row(p,
    if (p %in% c("rlang", "scales")) "Imports" else "Suggests; required only for documented optional path",
    if (found) as.character(utils::packageVersion(p)) else "NOT_INSTALLED",
    if (p %in% c("rlang", "scales")) "LOCAL_TRANSITIVE" else if (found) "INSTALLED_NOT_EXERCISED" else "NOT_RUN",
    "Presence is not behavior evidence for XLSX/Bioconductor/device/vignette paths")
}
rows[[length(rows) + 1L]] <- row("GenBank parser", "External Python + Biopython selected explicitly by reader",
  "not probed or installed by RC tools", "NOT_RUN", "No native parser or optional environment guarantee inferred")
write.csv(do.call(rbind, rows), "design-roadmap/audit-generated/rc-support-v0.16.6.csv", row.names = FALSE)
cat("Support declaration, installed versions and unknown combinations separated.\n")
