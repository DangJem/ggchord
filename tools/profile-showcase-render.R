# From the repository root:
# Rscript --vanilla tools/profile-showcase-render.R OUTPUT [timing|profile] [case ...]
# The output directory holds the ordinary PNGs and timings for later comparison.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply an output directory")
output <- args[1L]
mode <- if (length(args) >= 2L) match.arg(args[2L], c("timing", "profile")) else "timing"
sizes <- list(base = c(8, 5), dotplot = c(8, 5), chord = c(12, 10),
              collinear = c(12, 8), circular = c(16.5, 8))
cases <- if (length(args) > 2L) args[-(1:2)] else names(sizes)
stopifnot(all(cases %in% names(sizes)), !anyDuplicated(cases))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
  library(ggplot2)
  library(dplyr)
})

lines <- readLines("tools/documentation-examples.R", warn = FALSE)
markers <- grep("^# ---- [a-z]+ ----$", lines)
starts <- head(markers, -1L)
keys <- sub("^# ---- ([a-z]+) ----$", "\\1", lines[starts])
recipes <- setNames(lapply(seq_along(starts), function(i) {
  parse(text = lines[seq.int(starts[i] + 1L, markers[i + 1L] - 1L)])
}), keys)
for (key in c("sequencepalette", "docstyle", "showcasedata", "chordshowcasegenes",
              "showcaseplasmiddata", "dotplotdata")) invisible(eval(recipes[[key]]))
plots <- c(dotplot = "dotplotview", chord = "showcasechord",
           collinear = "showcasecollinear", circular = "showcaseplasmid")

elapsed <- function() unname(proc.time()[["elapsed"]])
render <- function(case, name, profile = FALSE) {
  size <- sizes[[case]]
  png_path <- file.path(output, paste0(name, ".png"))
  ragg::agg_png(png_path, width = size[1], height = size[2], units = "in",
                res = 144, background = "white")
  device <- grDevices::dev.cur()
  on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
  if (profile) {
    profile_path <- file.path(output, paste0(name, ".prof"))
    Rprof(profile_path, interval = .01)
    on.exit(Rprof(NULL), add = TRUE)
  }
  t0 <- elapsed()
  plot <- if (case == "base") ggchord(seq_data_example) + geom_seq() + coord_collinear() else
    eval(recipes[[plots[[case]]]], envir = globalenv())
  t1 <- elapsed()
  built <- ggplot_build(plot)
  t2 <- elapsed()
  table <- ggplot_gtable(built)
  t3 <- elapsed()
  grid::grid.newpage()
  grid::grid.draw(table)
  grDevices::dev.off(device)
  t4 <- elapsed()
  if (profile) {
    Rprof(NULL)
    samples <- summaryRprof(profile_path)
    write.csv(samples$by.total, file.path(output, paste0(name, "-cpu.csv")))
    write.csv(samples$by.self, file.path(output, paste0(name, "-self.csv")))
  }
  data.frame(case, png = basename(png_path), setup = t1 - t0, build = t2 - t1,
    gtable = t3 - t2, draw_close = t4 - t3, total = t4 - t0,
    refinements = attr(table, "ggchord_render_context")$refinement_passes,
    md5 = unname(tools::md5sum(png_path)))
}

invisible(render("base", paste0(mode, "-warm")))
results <- list()
for (i in seq_len(if (mode == "timing") 3L else 1L)) {
  for (case in if (i %% 2L == 0L) rev(cases) else cases) {
    row <- render(case, paste(case, mode, i, sep = "-"), mode == "profile")
    results[[length(results) + 1L]] <- row
    message(sprintf("%s %d  %.3fs", case, i, row$total))
  }
}
write.csv(do.call(rbind, results), file.path(output, paste0(mode, ".csv")), row.names = FALSE)
