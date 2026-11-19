# Draw each published script in a fresh session with a single local package installation.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1L]])
out <- normalizePath(args[[2L]], mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
specs <- jsonlite::fromJSON(file.path(root, "website/handbook-examples.json"), simplifyVector = FALSE)
if (length(args) < 3L || args[[3L]] != "--example") {
  ids <- if (length(args) > 2L) args[-c(1L, 2L)] else names(specs)
  lib <- tempfile("handbook-library-")
  dir.create(lib)
  install <- system2(file.path(R.home("bin"), "R"),
                     c("CMD", "INSTALL", "--no-docs", "--no-html", "-l", shQuote(lib), shQuote(root)))
  if (install != 0L) stop("Could not install the current package for the website examples")
  script <- file.path(root, "tools/website/draw-handbook.R")
  run <- function(id) {
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(script), shQuote(root), shQuote(out),
                        "--example", shQuote(lib), shQuote(id)))
    if (status == 0L) cat("Rendered", id, "\n")
    status
  }
  status <- if (.Platform$OS.type == "unix") {
    parallel::mclapply(ids, run, mc.cores = 2L, mc.preschedule = FALSE)
  } else {
    lapply(ids, run)
  }
  if (!all(vapply(status, identical, logical(1), 0L))) stop("Handbook example failed")
} else {
  .libPaths(c(args[[4L]], .libPaths()))
  id <- args[[5L]]
  spec <- specs[[id]]
  work <- tempfile("handbook-example-")
  dir.create(work)
  setwd(work)
  value <- source(file.path(root, "website/examples", paste0(id, ".R")), local = .GlobalEnv)$value
  if (spec$kind == "data") {
    writeLines(capture.output(print(value)), file.path(out, paste0(id, ".txt")), useBytes = TRUE)
  } else {
    if (!inherits(value, c("ggplot", "ggchord_fusion", "grob", "gtable"))) stop(id, ": expected a plot")
    background <- if (is.null(spec$background)) "white" else spec$background
    ggplot2::ggsave(file.path(out, paste0(id, ".png")), value, width = spec$width, height = spec$height,
                    dpi = spec$dpi, units = "in", bg = background, limitsize = FALSE)
  }
}
