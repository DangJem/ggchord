# Explicit synthetic detector controls, not biological observations.
source('tools/audit-source-snapshot.R')
audit_session <- audit_source_snapshot(additional =
  'tools/long-term-audit/check-capped-contour-quality.R')
pkgload::load_all(quiet = TRUE)
projection <- list(transform = function(data, panel) data)
built <- list(plot = list(coordinates = projection), layout = list(panel_params = list(list())))
context <- list(status = 'measured', panel_bounds_mm = c(xmin = 0, ymin = 0),
  panel_mm = c(width = 100, height = 100), canvas_mm = c(width = 100, height = 100))
text <- data.frame(layer_id = 'control', source_row = 1L, text_object = 'label',
  centre_x_mm = 40, centre_y_mm = 50, width_mm = 10, height_mm = 10,
  angle_radians = 0, xmin_mm = 35, xmax_mm = 45, ymin_mm = 45, ymax_mm = 55)
source_input <- data.frame(accver = 'schema_control', anno = 'label')
attr(source_input, 'ggchord_source') <- list(table = 'feature', rows = 1L)
layout <- list(layer_geometry = list(control = list(gene_poly = data.frame(
  x = c(.3, .5, .5, .3), y = c(.4, .4, .6, .6), group = 1L, source_row = 1L,
  .component = 'polygon', feature_shape = 'capped_line'))),
  layer_inputs = list(control = list(gene_poly = source_input,
    gene_label_repel = source_input)))
own <- ggchord_label_quality(layout, text, built, context)
stopifnot(own$summary$n_feature_polygons == 1L,
  own$summary$n_feature_parts_unchecked == 0L, own$summary$n_text_on_own_feature == 1L,
  own$summary$n_text_over_feature == 0L, own$summary$n_text_feature_unresolved == 0L)
attr(layout$layer_inputs$control$gene_label_repel, 'ggchord_source')$rows <- 2L
foreign <- ggchord_label_quality(layout, text, built, context)
stopifnot(foreign$summary$n_text_over_feature == 1L,
  foreign$summary$n_text_on_own_feature == 0L)
layout$layer_inputs <- list()
unresolved <- ggchord_label_quality(layout, text, built, context)
stopifnot(unresolved$summary$n_text_feature_unresolved == 1L)
layout$layer_geometry$control$gene_poly$.component <- 'boundary'
boundary <- ggchord_label_quality(layout, text, built, context)
stopifnot(boundary$summary$n_feature_polygons == 0L,
  boundary$summary$n_feature_parts_unchecked == 1L)
audit_assert_source_snapshot(audit_session)
jsonlite::write_json(list(status = 'passed',
  scope = 'Synthetic capped contour/source-identity controls, not a plot quality certification',
  own = own$summary, foreign = foreign$summary, unresolved = unresolved$summary,
  non_polygon_boundary = boundary$summary, runtime_sha256 = as.list(audit_session$engine),
  input_sha256 = as.list(audit_session$inputs)),
  'tools/long-term-audit/v0.27.0-capped-contour-controls.json',
  pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Capped contour own/foreign/unresolved and non-polygon boundary controls: passed\n')
