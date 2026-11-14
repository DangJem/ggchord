pkgload::load_all(quiet = TRUE)
render_source_sha <- digest::digest(file='tools/site-flagships.R',algo='sha256',serialize=FALSE)
source('tools/site-flagships.R')
items <- list(flagradial = flag_radial, flagsankey = flag_sankey,
  flagdotplot = flag_dotplot, flagfusion = flag_fusion)
selected <- commandArgs(trailingOnly = TRUE)
if (length(selected)) {
  stopifnot(all(selected %in% names(items)))
  items <- items[selected]
}
record <- 'tools/site-audit/flagship-rendering.json'
prior <- if (file.exists(record)) jsonlite::read_json(record) else list(figures = list())
results <- list()
for (id in names(items)) {
  cat('Rendering:', id, '\n'); flush.console()
  path <- paste0('man/figures/', id, '.png')
  size <- if (id == 'flagfusion') c(16, 12) else c(12, 9)
  ggplot2::ggsave(path, items[[id]], width = size[1], height = size[2],
    units = 'in', dpi = 120, bg = flagship_paper)
  stopifnot(identical(render_source_sha,digest::digest(file='tools/site-flagships.R',algo='sha256',serialize=FALSE)))
  results[[id]] <- list(file = path, sha256 = digest::digest(file = path, algo = 'sha256', serialize = FALSE),
    width = size[1], height = size[2], dpi = 120, background = flagship_paper,
    rendered_source_sha256 = render_source_sha)
  prior$figures[id] <- results[id]
  jsonlite::write_json(prior, record, auto_unbox = TRUE, pretty = TRUE)
}
record <- 'tools/site-audit/flagship-rendering.json'
prior <- if (file.exists(record)) jsonlite::read_json(record) else list(figures = list())
prior$figures[names(results)] <- results
prior$source_sha256 <- digest::digest(file = 'tools/site-flagships.R', algo = 'sha256', serialize = FALSE)
prior$inputs <- list(sequences = nrow(flagship_sequences), hsp = nrow(flagship_hsp), cds = nrow(flagship_features), retained_flagship_parent_data = nrow(show_links), retained_chord_features=nrow(chord_genes), retained_linear_features=nrow(show_linear_features), linear_displayed_hsp=nrow(show_adjacent)+nrow(show_cross))
prior$scope <- 'Existing Chord/Collinear design and biological data retained; Chord PNG separately refreshed for the default-style bug fix. Radial/Sankey derive from these parents, Dotplot/Fusion use expanded source tables; normal PNGs only, no browser or publication checks'
jsonlite::write_json(prior, record, auto_unbox = TRUE, pretty = TRUE)
