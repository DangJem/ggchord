# Single-device normal-size drawing-order regression, not package checks.
source('tools/audit-source-snapshot.R')
audit <- audit_source_snapshot(additional = 'tools/long-term-audit/check-axis-draw-order.R')
pkgload::load_all(quiet = TRUE)
before <- '--before' %in% commandArgs(trailingOnly = TRUE)
if (before) {
  # Isolate the previous automatic-axis insertion, with all geometry unchanged.
  previous <- new.env(parent = asNamespace('ggchord'))
  eval(parse(text = system2('git', c('show', 'e81a357:R/ggchord-frame.R'),
    stdout = TRUE)), previous)
  assignInNamespace('prepare_ggchord_plot', previous$prepare_ggchord_plot, 'ggchord')
}
data(plasmid_example_pUC19)
sequence <- plasmid_example_pUC19
features <- find_common_features(sequence)
plot <- ggchord(sequence[c('accver', 'length')], feature_data = features) +
  geom_seq(seq_style = 'double') + geom_feature_plasmid() +
  geom_feature_label_repel() + coord_circular() + theme_ggchord_plasmid()
original_layers <- vapply(plot$layers, `[[`, character(1), 'ggchord_type')
original_data <- plot$ggchord$data
folder <- 'tools/long-term-audit'
prefix <- paste0('v0.27.0-axis-draw-order-', if (before) 'before' else 'after')
file <- file.path(folder, paste0(prefix, '.png'))
ggchord_diagnostic_device(file, 'png', 10, 8, 110)
table <- ggplot2::ggplotGrob(plot)
grid::grid.newpage(); grid::grid.draw(table)
grDevices::dev.off()
built <- attr(table, 'ggchord_render_built')
layout <- attr(table, 'ggchord_render_layout')
layers <- built$plot$layers
axis <- which(vapply(layers, function(layer)
  inherits(layer$geom, 'GeomChordAxis'), logical(1)))
marks <- axis[vapply(layers[axis], function(layer)
  any(layer$data$.component %in% c('line', 'major_tick', 'minor_tick')), logical(1))]
text <- axis[vapply(layers[axis], function(layer)
  any(layer$data$.component %in% 'text'), logical(1))]
feature <- which(vapply(layers, function(layer)
  identical(layer$ggchord_type, 'gene_poly'), logical(1)))
types <- vapply(layers, function(layer) layer$ggchord_type %||% '', character(1))
stopifnot(identical(unname(types[nzchar(types)]), unname(original_layers)),
  identical(plot$ggchord$data, original_data))
if (before) stopifnot(min(marks) > max(feature)) else
  stopifnot(max(marks) < min(feature), min(text) > max(feature))
saveRDS(list(axis = ggchord_axis_geometry(layout), inputs = layout$layer_inputs,
  geometry = layout$layer_geometry, source = plot$ggchord$data),
  file.path('/tmp', paste0(prefix, '.rds')))
audit_assert_source_snapshot(audit)
jsonlite::write_json(list(scope = if (before)
  'Previous insertion function from e81a357 evaluated in isolation; current shared geometry' else
  'Current native implementation; actual 10x8in/110dpi PNG',
  file = file, source_features = nrow(features),
  drawing_order = list(axis_marks = marks, feature = feature, axis_text = text),
  user_layer_order_preserved = TRUE,
  physical_converged = attr(table, 'ggchord_render_context')$converged,
  runtime_sha256 = as.list(audit$engine)),
  file.path(folder, paste0(prefix, '-results.json')), pretty = TRUE, auto_unbox = TRUE)
cat(prefix, ': drawing order and source controls passed\n')
