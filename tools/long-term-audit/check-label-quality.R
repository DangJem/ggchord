# One PNG per requested case; no release device matrix or package test suite.
pkgload::load_all(quiet = TRUE)
library(dplyr)
folder <- 'tools/long-term-audit'
arguments <- commandArgs(trailingOnly = TRUE)
features_check <- '--features' %in% arguments
identity_check <- '--identity' %in% arguments
clearance_check <- '--clearance' %in% arguments
after <- '--after' %in% arguments
prefix <- if (clearance_check) 'v0.27.0-quality-clearance-' else if (identity_check) 'v0.27.0-quality-identity-' else if (features_check) 'v0.27.0-quality-features-' else if (after) 'v0.27.0-quality-after-' else 'v0.27.0-quality-'
measure <- function(plot, id, width = 10, height = 8) {
  path <- file.path(folder, paste0(prefix, id, '.png'))
  ggchord_diagnostic_device(path, 'png', width, height, 110)
  on.exit(grDevices::dev.off())
  table <- ggplot2::ggplotGrob(plot)
  grid::grid.newpage(); grid::grid.draw(table)
  layout <- attr(table, 'ggchord_render_layout')
  context <- attr(table, 'ggchord_render_context')
  built <- attr(table, 'ggchord_render_built')
  texts <- lapply(seq_along(built$data), function(i) {
    d <- built$data[[i]]; layer <- built$plot$layers[[i]]
    if (is.data.frame(layer$data) && nrow(layer$data) == nrow(d) && 'source_row' %in% names(layer$data))
      d$source_row <- layer$data$source_row
    d <- ggchord_diagnostic_text_styles(d, layer)
    text <- ggchord_diagnostic_text(d, built$plot$coordinates, built$layout$panel_params[[1L]], context,
      parse = layer$geom_params$parse %||% FALSE, size_unit = layer$geom_params$size.unit %||% 'mm')
    if (!nrow(text)) return(NULL)
    text$layer_id <- layer$ggchord_layer_id %||% paste0('native-', i)
    text
  })
  text <- ggchord_rbind_fill(texts) %||% data.frame()
  result <- ggchord_label_quality(layout, text, built, context)
  result$label_visibility <- ggchord_label_visibility(layout)
  axis <- ggchord_axis_geometry(layout)
  if ('.axis_clearance_status' %in% names(axis)) result$axis_clearance <- unique(axis[
    axis$.component %in% 'text' & !is.na(axis$.axis_clearance_status),
    intersect(c('accver', '.axis_label', '.axis_offset_explicit', '.axis_clearance_status', '.axis_clearance_mm'), names(axis)), drop = FALSE])
  jsonlite::write_json(result, file.path(folder, paste0(prefix, id, '.json')),
    auto_unbox = TRUE, pretty = TRUE, na = 'null')
  c(result$summary, list(width = width, height = height, file = path,
    measurement = context$status, visibility_reasons = as.list(table(ggchord_label_visibility(layout)$reason))))
}
# Explicit geometry controls test detectors rather than invented observations.
identity_projection <- list(transform = function(data, panel) data)
context <- list(status = 'measured', panel_bounds_mm = c(xmin = 0, ymin = 0),
  panel_mm = c(width = 100, height = 100), canvas_mm = c(width = 100, height = 100))
mock <- list(plot = list(coordinates = identity_projection), layout = list(panel_params = list(list())))
empty <- list(layer_geometry = list(), layer_inputs = list())
boxes <- data.frame(layer_id = 'format', source_row = 1:2, text_object = c('one', 'two'),
  centre_x_mm = c(40, 43), centre_y_mm = 50, width_mm = 10, height_mm = 10,
  angle_radians = 0, xmin_mm = c(35, 38), xmax_mm = c(45, 48), ymin_mm = 45, ymax_mm = 55)
result <- ggchord_label_quality(empty, boxes, mock, context)
stopifnot(result$summary$n_text_overlaps == 1L, result$summary$n_feature_polygons == 0L)
empty$layer_geometry <- list(format = list(gene_label_repel = data.frame(
  .component = 'segment', source_row = 1:2, accver = 'format',
  x = c(.1, .1), y = c(.2, .8), xend = c(.8, .8), yend = c(.8, .2))))
result <- ggchord_label_quality(empty, boxes, mock, context)
stopifnot(result$summary$n_leader_crossings == 1L, result$summary$n_leaders_through_text > 0L)
empty$layer_geometry$format$gene_label_repel$xend <- empty$layer_geometry$format$gene_label_repel$x
empty$layer_geometry$format$gene_label_repel$yend <- empty$layer_geometry$format$gene_label_repel$y
result <- ggchord_label_quality(empty, boxes, mock, context)
stopifnot(result$summary$n_leaders == 0L, result$summary$n_leaders_through_text == 0L)
# Synthetic polygon controls: same constructor source is intentional, a different
# constructor row is foreign, and absent association stays unresolved.
empty$layer_geometry <- list(format = list(gene_poly = data.frame(
  x = c(.3, .5, .5, .3), y = c(.4, .4, .6, .6), group = 1L,
  source_row = 1L, .component = 'polygon')))
feature_input <- data.frame(accver = 'format', anno = 'control')
attr(feature_input, 'ggchord_source') <- list(table = 'feature', rows = 1L)
empty$layer_inputs <- list(format = list(gene_poly = feature_input,
  gene_label_repel = feature_input))
result <- ggchord_label_quality(empty, boxes[1L, ], mock, context)
stopifnot(result$summary$n_text_on_own_feature == 1L, result$summary$n_text_over_feature == 0L)
attr(empty$layer_inputs$format$gene_label_repel, 'ggchord_source')$rows <- 2L
result <- ggchord_label_quality(empty, boxes[1L, ], mock, context)
stopifnot(result$summary$n_text_over_feature == 1L)
empty$layer_inputs <- list()
result <- ggchord_label_quality(empty, boxes[1L, ], mock, context)
stopifnot(result$summary$n_text_feature_unresolved == 1L, result$summary$n_text_over_feature == 0L)
empty$layer_geometry$format$gene_poly$.component <- 'boundary'
result <- ggchord_label_quality(empty, boxes[1L, ], mock, context)
stopifnot(result$summary$n_feature_polygons == 0L, result$summary$n_feature_parts_unchecked == 1L)
objects <- c('pSB1C3', 'pUC19', 'pBR322', 'pBluescript_II_SK_plus', 'pET_28a_plus',
  'pETDuet_1', 'pDONR221', 'pcDNA3_1_plus', 'pTRE_Tight_BI', 'pTRIPZ', 'pCAMBIA1300',
  'pEarleyGate_201', 'pSpCas9_BB_2A_GFP_PX458')
requested <- setdiff(arguments, c('--after', '--features', '--identity', '--clearance'))
if (length(requested)) objects <- intersect(objects, requested)
results <- list()
for (id in objects) {
  cat('Quality:', id, '\n'); flush.console()
  environment <- new.env()
  data(list = paste0('plasmid_example_', id), package = 'ggchord', envir = environment)
  sequence <- environment[[paste0('plasmid_example_', id)]]
  features <- find_common_features(sequence)
  sites <- find_restriction_sites(sequence)
  nearest <- if (nrow(sites)) unique(vapply(sequence$length*c(1,3,5)/6,
    function(x) which.min(abs(sites$position-x)), integer(1))) else integer()
  position <- position_feature_stack(spacing = .085, base_position = position_plasmid())
  plot <- ggchord(sequence[c('accver', 'length')], feature_data = features) + geom_seq(seq_style = 'double') +
    geom_feature_plasmid(position = position) + geom_feature_label_repel(position = position,
      external = TRUE, label_fit = 'shrink', label_min_size = 2.5) +
    geom_restriction_site(data = sites[nearest, ], label_size = 2.5) + coord_circular() + theme_ggchord_plasmid()
  original <- plot$ggchord$data
  results[[id]] <- measure(plot, id)
  stopifnot(identical(original, plot$ggchord$data))
}
if (clearance_check) {
  before <- jsonlite::fromJSON(file.path(folder, 'v0.27.0-label-quality-feature-results.json'), simplifyVector = FALSE)$cases
  for (id in names(results)) {
    current <- results[[id]]
    stopifnot(current$n_source_labels == before[[id]]$n_source_labels,
      current$n_visible_source_labels == current$n_source_labels,
      current$n_text_overlaps == 0L, current$n_leader_crossings == 0L,
      current$n_leaders_through_text == 0L, current$n_label_outside_canvas == 0L,
      current$n_text_over_feature == 0L, current$n_text_feature_unresolved == 0L)
  }
}
source_files <- list.files('R', '[.]R$', full.names = TRUE)
jsonlite::write_json(list(scope = 'Original catalogue inputs; 10x8 inches, one PNG device. Contacts reported, not assumed zero; not a cross-application comparison.',
  source_sha256 = as.list(stats::setNames(vapply(source_files,
    function(path) digest::digest(file = path, algo = 'sha256'), character(1)), source_files)),
  detector_controls = 'Explicit synthetic geometry distinguishes overlap, crossing, through-text and own/foreign/unresolved polygon contacts',
  cases = results, comparison_scope = 'SnapGene is a visual reference only; matched external exports are not required under the clarified user instruction.',
  unverified = c('non-polygon strokes/markers, visibility/occlusion and exact ink', 'release device/platform matrix')),
  file.path(folder, if (clearance_check) 'v0.27.0-label-quality-clearance-results.json' else if (identity_check) 'v0.27.0-label-quality-identity-results.json' else if (features_check) 'v0.27.0-label-quality-feature-results.json' else if (after) 'v0.27.0-label-quality-after-results.json' else 'v0.27.0-label-quality-results.json'),
  auto_unbox = TRUE, pretty = TRUE, na = 'null')
cat('Label quality observations recorded\n')
