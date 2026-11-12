# Targeted source-backed size profiles on one PNG device, not release checks.
pkgload::load_all(quiet = TRUE)
library(dplyr)
arguments <- commandArgs(trailingOnly = TRUE)
after <- '--after' %in% arguments
suffix <- if (after) '-after' else ''
plant <- read_delta('examples/genome-alignments/arabidopsis-first-pair.delta',
  executable = Sys.getenv('DELTA_DECODER', 'show-coords'))
base <- ggchord(plant$seq_data, ribbon_data = plant$ribbon_data) +
  geom_seq() + geom_seq_label(size = 3) + geom_link_ribbon()
data(plasmid_example_pUC19)
annotations <- find_common_features(plasmid_example_pUC19)
circle <- ggchord(plasmid_example_pUC19[c('accver', 'length')], feature_data = annotations) +
  geom_seq(seq_style = 'double') + geom_feature_plasmid() + geom_feature_label_repel() +
  coord_circular() + theme_ggchord_plasmid()
plots <- list(chord = base + coord_chord(), collinear = base + coord_collinear(),
  vertical = base + coord_collinear(seq_axis = 'vertical'),
  radial = base + coord_radial_synteny(), sankey = base + coord_sankey(),
  dotplot = ggchord(plant$seq_data, ribbon_data = plant$ribbon_data) +
    geom_alignment_segment() + coord_dotplot(), circular = circle)
measure <- function(plot, id, width, height, keep) {
  path <- if (keep) paste0('tools/long-term-audit/v0.27.0-adaptive-quality-', id, suffix, '.png') else tempfile(fileext = '.png')
  ggchord_diagnostic_device(path, 'png', width, height, 100)
  on.exit({grDevices::dev.off(); if (!keep) unlink(path)})
  original <- plot$ggchord$data
  table <- ggplot2::ggplotGrob(plot)
  grid::grid.newpage(); grid::grid.draw(table)
  built <- attr(table, 'ggchord_render_built'); layout <- attr(table, 'ggchord_render_layout')
  context <- attr(table, 'ggchord_render_context')
  text <- ggchord_rbind_fill(lapply(seq_along(built$data), function(i) {
    data <- built$data[[i]]; layer <- built$plot$layers[[i]]
    if (is.data.frame(layer$data) && nrow(layer$data) == nrow(data) && 'source_row' %in% names(layer$data))
      data$source_row <- layer$data$source_row
    txt <- ggchord_diagnostic_text(ggchord_diagnostic_text_styles(data, layer),
      built$plot$coordinates, built$layout$panel_params[[1L]], context,
      parse = layer$geom_params$parse %||% FALSE, size_unit = layer$geom_params$size.unit %||% 'mm')
    if (!nrow(txt)) return(NULL)
    txt$layer_id <- layer$ggchord_layer_id %||% paste0('native-', i)
    txt
  })) %||% data.frame()
  quality <- ggchord_label_quality(layout, text, built, context)
  stopifnot(identical(original, plot$ggchord$data), identical(context$status, 'measured'))
  if (id %in% c('chord', 'collinear', 'vertical', 'radial', 'sankey')) {
    labels <- text |> filter(label %in% plant$seq_data$accver)
    stopifnot(nrow(labels) > 0L,
      isTRUE(all.equal(labels$size_mm, rep(3, nrow(labels)), check.attributes = FALSE)))
  }
  list(width = width, height = height, measurement = context$status,
    guide_layout = context$guide_layout, adaptive_layout = context$adaptive_layout,
    quality = quality$summary, contacts = quality$conflicts, text = quality$text,
    source_preserved = TRUE, explicit_geom_text_preserved = id != 'dotplot' && id != 'circular',
    file = if (keep) path else NULL, not_checked = quality$not_checked)
}
profiles <- list()
for (id in names(plots)) {
  cat('Adaptive quality:', id, '\n'); flush.console()
  p <- plots[[id]] + adapt_ggchord('adaptive')
  profiles[[id]] <- list(landscape = measure(p, id, 9, 6, TRUE),
    portrait = measure(p, id, 6, 9, FALSE))
}
jsonlite::write_json(list(scope = 'Real unchanged plant intervals and pUC19 annotations; six coordinate modes plus vertical Collinear; landscape/portrait single-device profiles',
  profiles = profiles, runtime_sha256 = as.list(stats::setNames(vapply(list.files('R', '[.]R$', full.names = TRUE),
    function(path) digest::digest(file = path, algo = 'sha256'), character(1)), list.files('R', '[.]R$', full.names = TRUE))),
  unverified = c('complete publication quality for arbitrary view/size combinations',
    'all geometry-constant scaling, exact ink and font portability', 'release device/platform matrix')),
  paste0('tools/long-term-audit/v0.27.0-adaptive-quality', suffix, '-results.json'), pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Adaptive quality profiles recorded; contact counts are observations, not a blanket pass\n')
