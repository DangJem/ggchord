# Draw each published script in a fresh session with a single local package installation.
check_fusion_strand <- function(composition) {
  find_grobs <- function(grob, match) {
    if (match(grob)) return(list(grob))
    children <- if (inherits(grob, "gtable")) grob$grobs else grob$children
    unlist(lapply(children, find_grobs, match = match), recursive = FALSE)
  }
  text_labels <- function(grob) {
    texts <- find_grobs(grob, function(grob) inherits(grob, "text"))
    unlist(lapply(texts, function(text) as.character(text$label)), use.names = FALSE)
  }
  guides <- find_grobs(composition, function(grob) {
    if (!inherits(grob, "gtable")) return(FALSE)
    titles <- grob$grobs[grob$layout$name == "title"]
    any(vapply(titles, function(title) "Strand" %in% text_labels(title), logical(1)))
  })
  if (length(guides) != 1L) stop("case-fusion: expected one visible Strand legend")
  guide <- guides[[1L]]
  labels <- unlist(lapply(guide$grobs[startsWith(guide$layout$name, "label-")], text_labels))
  keys <- guide$grobs[startsWith(guide$layout$name, "key-")]
  if (length(keys) != 2L || length(labels) != 2L || !setequal(labels, c("+", "-")))
    stop("case-fusion: Strand must have one positive and one negative key")
  for (i in seq_along(keys)) {
    marks <- find_grobs(keys[[i]], function(grob) {
      !inherits(grob, c("gTree", "null", "zeroGrob")) &&
        !startsWith(grob$name, "legend.key.")
    })
    if (length(marks) != 1L || !inherits(marks[[1L]], "polygon"))
      stop("case-fusion: Strand key ", labels[i], " must draw only its gene arrow")
    mark <- marks[[1L]]
    fill <- mark$gp$fill
    if (!length(fill) || !any(grDevices::col2rgb(fill, alpha = TRUE)[4L, ] > 0L))
      stop("case-fusion: Strand key ", labels[i], " must show a filled gene arrow")
    x <- unique(data.frame(x = as.numeric(mark$x), y = as.numeric(mark$y)))$x
    left <- sum(x == min(x))
    right <- sum(x == max(x))
    points <- if (labels[i] == "+") right == 1L && left > 1L else left == 1L && right > 1L
    if (!points) stop("case-fusion: Strand key ", labels[i], " points in the wrong direction")
  }
}

checked_fusion_grob <- function(plot) {
  # Build and check the same grob on the actual export device before drawing.
  grid::recordGrob({
    grob <- build_fusion(fusion)$grob
    check_strand(grob)
    grid::grid.draw(grob, recording = FALSE)
  }, list(fusion = plot, build_fusion = getFromNamespace("ggchord_build_fusion", "ggchord"),
    check_strand = check_fusion_strand))
}

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
    if (identical(id, "case-fusion")) value <- checked_fusion_grob(value)
    background <- if (is.null(spec$background)) "white" else spec$background
    ggplot2::ggsave(file.path(out, paste0(id, ".png")), value, width = spec$width, height = spec$height,
                    dpi = spec$dpi, units = "in", bg = background, limitsize = FALSE)
  }
}
