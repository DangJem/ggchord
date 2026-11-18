# From the repository root:
# Rscript --vanilla tools/long-term-audit/v0.50.0-render.R OUTPUT timing|profile|views
# Rscript --vanilla tools/long-term-audit/v0.50.0-render.R BASELINE compare|compare-migrated CANDIDATE
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L,
  args[2] %in% c("timing", "profile", "views", "compare", "compare-migrated"))
if (args[2] %in% c("compare", "compare-migrated")) {
  stopifnot(length(args) == 3L)
  folders <- args[c(1L, 3L)]
  results <- lapply(folders, function(folder) rbind(
    read.csv(file.path(folder, "timing.csv")),
    read.csv(file.path(folder, "views.csv"))
  ))
  metadata <- c("case", "png", "width", "height", "dpi")
  stopifnot(identical(results[[1L]][metadata], results[[2L]][metadata]))
  checksums <- lapply(folders, function(folder)
    unname(tools::md5sum(file.path(folder, results[[1L]]$png))))
  stopifnot(!anyNA(checksums[[1L]]), identical(checksums[[1L]], checksums[[2L]]))
  check_prepared <- function(old, new, path) {
    if (identical(old, new)) return(invisible(NULL))
    if (is.data.frame(old) && is.data.frame(new)) {
      added <- setdiff(names(new), names(old))
      stopifnot(nrow(old) == nrow(new), all(names(old) %in% names(new)),
        all(added %in% c(".source_row", "metadata_.source_row", "source_strand", "directionality",
          "anchor_layer")))
      if ("anchor_layer" %in% added) {
        stopifnot(identical(path, "circular$gene_labels"),
          all(new$anchor_layer %in% c(NA_character_, "layer-0002")))
      }
      # Preparation assigns a local row index and preserves the previous
      # source index as metadata when that column already exists.
      if (".source_row" %in% names(old) && "metadata_.source_row" %in% added) {
        new$.source_row <- new$metadata_.source_row
      }
      for (field in intersect(names(old), c("position_name", ".position_name"))) {
        if (identical(old[[field]], new[[field]])) next
        previous_names <- c(plasmid = "offset", feature_stack_plasmid = "feature_stack_offset",
          feature_stack_legacy = "feature_stack_strand")
        renamed <- old[[field]] %in% names(previous_names)
        old[[field]][renamed] <- unname(previous_names[old[[field]][renamed]])
      }
      if (path %in% c("circular$gene_polys", "circular$gene_labels")) {
        inputs <- new.env(parent = emptyenv())
        load("data/flagship_circular_example.rda", envir = inputs)
        features <- inputs$flagship_circular_example$dense$features
        source <- match(new$match_id, features$match_id)
        at <- !is.na(source)
        for (field in intersect(names(old), c("type", "repeat_type"))) {
          if (identical(old[[field]], new[[field]])) next
          stopifnot(all(old[[field]][at] == features$anno[source[at]]))
          old[[field]][at] <- features$type[source[at]]
        }
      }
      for (field in names(old)) {
        if (!identical(old[[field]], new[[field]])) stop(path, "$", field, " differed")
      }
      return(invisible(NULL))
    }
    stopifnot(is.list(old), is.list(new), identical(names(old), names(new)))
    for (field in names(old)) {
      check_prepared(old[[field]], new[[field]], paste(path, field, sep = "$"))
    }
  }
  for (case in unique(results[[1L]]$case)) {
    layouts <- lapply(file.path(folders, paste0(case, ".rds")), readRDS)
    if (args[2] == "compare") {
      stopifnot(identical(layouts[[1L]], layouts[[2L]]))
    } else {
      check_prepared(layouts[[1L]], layouts[[2L]], case)
    }
    message(case, " — checked device PNGs and layout snapshot")
  }
  quit(save = "no")
}
output <- args[1]
mode <- args[2]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
  library(ggplot2)
  library(dplyr)
})
r_files <- sort(list.files("R", full.names = TRUE))
inputs <- c(r_files, "tools/documentation-examples.R", "tools/site-flagships.R",
            sort(list.files("data", full.names = TRUE)))
write.csv(data.frame(file = inputs, md5 = unname(tools::md5sum(inputs))),
          file.path(output, paste0(mode, "-inputs.csv")), row.names = FALSE)
writeLines(c(system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
             paste("ragg", utils::packageVersion("ragg")),
             capture.output(sessionInfo())),
           file.path(output, paste0(mode, "-environment.txt")))
lines <- readLines("tools/documentation-examples.R", warn = FALSE)
markers <- grep("^# ---- [a-z]+ ----$", lines)
keys <- sub("^# ---- ([a-z]+) ----$", "\\1", head(lines[markers], -1L))
recipes <- setNames(lapply(seq_along(keys), function(i) {
  parse(text = lines[seq.int(markers[i] + 1L, markers[i + 1L] - 1L)])
}), keys)
for (key in c("sequencepalette", "docstyle", "showcasedata", "chordshowcasegenes",
              "showcaseplasmiddata", "dotplotdata")) invisible(eval(recipes[[key]]))
if (mode == "views") source("tools/site-flagships.R")

elapsed <- function() unname(proc.time()[["elapsed"]])
render <- function(case, name, capture = FALSE, profile = FALSE) {
  size <- switch(case, base = c(8, 5), circular = c(16.5, 8), radial = c(12, 9),
                 sankey = c(12, 9), fusion = c(16, 14))
  background <- if (case %in% c("radial", "sankey", "fusion")) flagship_paper else "white"
  dpi <- if (case %in% c("radial", "sankey", "fusion")) 120 else 144
  png_path <- file.path(output, paste0(name, ".png"))
  ragg::agg_png(png_path, width = size[1], height = size[2], units = "in",
                res = dpi, background = background)
  device <- grDevices::dev.cur()
  on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
  if (profile) {
    path <- tempfile(fileext = ".prof")
    Rprof(path, interval = .01)
    on.exit(Rprof(NULL), add = TRUE)
  }
  t0 <- elapsed()
  p <- switch(case, base = ggchord(seq_data_example) + geom_seq() + coord_collinear(),
    circular = eval(recipes[["showcaseplasmid"]], envir = globalenv()),
    radial = flag_radial, sankey = flag_sankey, fusion = flag_fusion)
  t1 <- elapsed()
  if (case == "fusion") {
    result <- ggchord:::ggchord_build_fusion(p)
    t2 <- t3 <- elapsed()
    table <- result$grob
  } else {
    built <- ggplot_build(p)
    t2 <- elapsed()
    table <- ggplot_gtable(built)
    t3 <- elapsed()
  }
  grid::grid.newpage()
  grid::grid.draw(table)
  grDevices::dev.off(device)
  t4 <- elapsed()
  if (profile) {
    Rprof(NULL)
    samples <- summaryRprof(path)
    write.csv(samples$by.total, file.path(output, paste0(name, "-cpu.csv")))
    write.csv(samples$by.self, file.path(output, paste0(name, "-self.csv")))
  }
  context <- attr(table, "ggchord_render_context")
  if (capture) {
    if (case == "fusion") {
      # Fusion's public export accompanies this exact grob and draw.
      snapshot <- result$export[c("panels", "instances", "connections", "fronts",
        "paths", "focus_regions", "focus_fronts", "source_data", "relationship_data",
        "source_observations")]
    } else {
      # Read measured geometry from the draw above, before it can be replaced
      # by a layout-only solve using estimated panel dimensions.
      layout <- attr(table, "ggchord_render_layout")
      snapshot <- layout[vapply(layout, is.data.frame, logical(1))]
      snapshot$render <- context[c("panel_mm", "refinement_passes")]
    }
    saveRDS(snapshot, file.path(output, paste0(case, ".rds")))
  }
  data.frame(case, png = basename(png_path), width = size[1], height = size[2], dpi,
    setup = t1 - t0, build = t2 - t1, gtable = t3 - t2, draw_close = t4 - t3,
    total = t4 - t0, refinements = if (is.null(context)) NA else context$refinement_passes,
    md5 = unname(tools::md5sum(png_path)))
}

invisible(render("base", paste0(mode, "-warm")))
cases <- switch(mode, timing = rep("circular", 3L), profile = "circular",
                views = c("radial", "sankey", "fusion"))
results <- lapply(seq_along(cases), function(i) {
  case <- cases[i]
  result <- render(case, paste(case, mode, i, sep = "-"),
    capture = mode == "views" || mode == "timing" && i == 1L,
    profile = mode == "profile")
  message(sprintf("%s %d  %.3fs", case, i, result$total))
  result
})
write.csv(do.call(rbind, results), file.path(output, paste0(mode, ".csv")), row.names = FALSE)
