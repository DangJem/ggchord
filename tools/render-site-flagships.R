# Run from the repository root. Optional arguments select home preview IDs.
pkgload::load_all(quiet = TRUE)
source("tools/site-flagships.R")
source("tools/home-previews.R")
previews <- jsonlite::read_json("website/home-previews.json")
selected <- commandArgs(trailingOnly = TRUE)
known <- vapply(previews, "[[", character(1), "id")
if (length(selected) && !all(selected %in% known))
  stop("Select homepage previews by ID: ", paste(known, collapse = ", "))
for (row in previews) {
  if (length(selected) && !row$id %in% selected) next
  ggplot2::ggsave(row$file, get(row$object), width = row$width, height = row$height,
    units = "in", dpi = row$dpi, bg = home_paper)
  message("Rendered ", row$id)
}
