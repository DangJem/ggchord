# Local Fusion maintenance rendering; not a package or release check.
pkgload::load_all(quiet = TRUE)
source("tools/site-flagships.R")
source("tools/home-previews.R")
inputs <- c("R/geom-residue.R", "R/scale-residue.R", "R/alignment-panel.R", "R/magnify.R",
  "R/fusion-routing.R", "R/fusion.R", "R/fusion-render.R", "DESCRIPTION", "tools/site-flagships.R", "tools/home-previews.R")
snapshot <- function() setNames(lapply(inputs, function(path)
  digest::digest(file = path, algo = "sha256")), inputs)
before <- snapshot()
rows <- list(
  flagfusion = list(plot = flag_fusion, width = 16, height = 14, dpi = 120),
  homefusion = list(plot = home_fusion, width = 10, height = 10, dpi = 160)
)
record <- list(status = "rendering", generator = "tools/site-audit/render-residue-fusion.R",
  scope = "Two normal Fusion PNGs and source/leader checks; no package or browser checks", figures = list())
for (id in names(rows)) {
  row <- rows[[id]]
  path <- paste0("man/figures/", id, ".png")
  ggplot2::ggsave(path, row$plot, width = row$width, height = row$height,
    dpi = row$dpi, bg = flagship_paper)
  report <- export_ggchord_fusion(row$plot, width = row$width, height = row$height)
  stopifnot(nrow(report$panels) == 6L, nrow(report$connections) == 3L,
    nrow(report$focus_regions) == 3L, all(report$leaders$obstacle_collisions == 0),
    all(fusion_base_cells$source_verified),
    all(report$focus_regions$frame_start == report$focus_regions$start),
    all(report$focus_regions$frame_end == report$focus_regions$end), identical(before, snapshot()))
  record$figures[[id]] <- list(file = path, width = row$width, height = row$height, dpi = row$dpi,
    sha256 = digest::digest(file = path, algo = "sha256"), panels = nrow(report$panels),
    connections = nrow(report$connections), focus_regions = report$focus_regions, leaders = report$leaders)
  cat("Rendered", id, "\n")
}
record$status <- "complete"
record$base_selection <- fusion_base_cds
record$source_positions <- range(fusion_base_cells$source_position)
record$source_sha256 <- before
jsonlite::write_json(record, "tools/site-audit/residue-fusion-v0.28.0.json", auto_unbox = TRUE, pretty = TRUE)
