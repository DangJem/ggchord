# Render website display files from the documented records and selected layers.
# Historical audit images stay in their original directories.
pkgload::load_all(quiet = TRUE)
source("tools/site-flagships.R")
source("tools/home-previews.R")
home <- jsonlite::read_json("website/home-previews.json")
for (row in home) {
  ggplot2::ggsave(row$file, get(row$object), width = row$width, height = row$height,
    units = "in", dpi = row$dpi, bg = home_paper)
}
rows <- jsonlite::read_json("website/reference-catalog.json")
source_code <- readLines("tools/site-course-examples.R", warn = FALSE)
starts <- grep("^# ---- [a-z]+ ----$", source_code)
names(starts) <- sub("^# ---- ([a-z]+) ----$", "\\1", source_code[starts])
page <- readLines("website/zh/plasmid-catalog.md", warn = FALSE)
for (row in rows) {
  figure <- match(paste0("{{figure:", row$figure, "}}"), page)
  snippets <- grep("^\\{\\{code:", page[seq_len(figure - 1L)], value = TRUE)
  keys <- strsplit(sub("^\\{\\{code:([a-z,]+)\\}\\}$", "\\1", tail(snippets, 1L)), ",")[[1L]]
  code <- unlist(lapply(keys, function(key) {
    start <- unname(starts[[key]])
    end <- min(starts[starts > start]) - 1L
    source_code[(start + 1L):end]
  }), use.names = FALSE)
  exprs <- parse(text = code)
  env <- new.env(parent = globalenv())
  # The last expression is the page's sized preview. Use its plot in the PNG device.
  stopifnot(identical(exprs[[length(exprs)]][[1L]], as.name("view_ggchord")))
  eval(exprs[-length(exprs)], envir = env)
  id <- sub("^catalog", "complete", row$figure)
  ggplot2::ggsave(file.path("man/figures", paste0(id, ".png")), env$p,
    width = 10, height = 8, units = "in", dpi = 110, bg = "white")
  cat("Rendered", row$accver, "\n")
}
