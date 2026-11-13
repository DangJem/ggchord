# Tools-only observation: production diagnostics/layout remain unchanged.
source('tools/audit-source-snapshot.R')
point_mode <- '--points' %in% commandArgs(trailingOnly = TRUE)
point_source <- if (point_mode) system.file('extdata', 'gridss.vcf',
  package = 'StructuralVariantAnnotation') else character()
if (point_mode && !nzchar(point_source)) stop('The existing optional point source is unavailable')
audit_session <- audit_source_snapshot(additional =
  c('tools/long-term-audit/observe-feature-contours.R', point_source))
pkgload::load_all(quiet = TRUE)
audit_assert_source_snapshot(audit_session)
data(plasmid_example_pUC19)
sequence <- plasmid_example_pUC19
features <- find_common_features(sequence)
if (point_mode) {
  imported <- read_vcf(point_source, structural = TRUE)
  features <- imported$features$variant
  stopifnot(nrow(features) > 0L, all(features$start == features$end))
  # Display IDs verbatim; no variant or biological interpretation is performed.
  features$anno <- features$feature_id
  sequence <- imported$report$records[[1L]]$sequence_info
  sequence <- sequence[sequence$accver %in% features$accver, c('accver', 'length')]
  stopifnot(all(is.finite(sequence$length)), all(sequence$length > 0))
}
# Inspect the same detector without its capped-line exclusion, in an isolated
# function. Closed transparent contours are conservative envelopes, not ink.
detector <- ggchord_label_feature_quality
lines <- deparse(body(detector), width.cutoff = 500L)
at <- grep('is_polygon <- is_polygon & !part\\$feature_shape', lines)
stopifnot(length(at) == 1L, grepl('if \\("feature_shape"', lines[at - 1L]))
body(detector) <- parse(text = paste(lines[-c(at - 1L, at)], collapse = '\n'))[[1L]]

observe <- function(shape) {
  position <- if (point_mode) position_identity() else
    position_feature_stack(spacing = .085, base_position = position_plasmid())
  feature_layer <- if (point_mode) geom_feature(position = position) else if (shape == 'preset') geom_feature_plasmid(position = position) else
    geom_feature(shape = shape, position = position)
  plot <- ggchord(sequence[c('accver', 'length')], feature_data = features) +
    geom_seq(seq_style = 'double') + feature_layer +
    geom_feature_label_repel(position = position, external = TRUE, label_fit = 'shrink',
      label_min_size = 2.5)
  plot <- if (point_mode) plot + coord_collinear() + theme_ggchord() else
    plot + coord_circular() + theme_ggchord_plasmid()
  original <- plot$ggchord$data
  file <- if (point_mode) 'tools/long-term-audit/v0.27.0-variant-point-contours.png' else
    tempfile(fileext = '.png')
  ggchord_diagnostic_device(file, 'png', 10, 8, 110)
  on.exit({grDevices::dev.off(); if (!point_mode) unlink(file)})
  table <- ggplot2::ggplotGrob(plot)
  grid::grid.newpage(); grid::grid.draw(table)
  layout <- attr(table, 'ggchord_render_layout')
  context <- attr(table, 'ggchord_render_context')
  built <- attr(table, 'ggchord_render_built')
  text <- ggchord_rbind_fill(lapply(seq_along(built$data), function(i) {
    data <- built$data[[i]]; layer <- built$plot$layers[[i]]
    if (is.data.frame(layer$data) && nrow(layer$data) == nrow(data) &&
        'source_row' %in% names(layer$data)) data$source_row <- layer$data$source_row
    txt <- ggchord_diagnostic_text(ggchord_diagnostic_text_styles(data, layer),
      built$plot$coordinates, built$layout$panel_params[[1L]], context,
      parse = layer$geom_params$parse %||% FALSE,
      size_unit = layer$geom_params$size.unit %||% 'mm')
    if (!nrow(txt)) return(NULL)
    txt$layer_id <- layer$ggchord_layer_id %||% paste0('native-', i)
    txt
  }))
  quality <- ggchord_label_quality(layout, text, built, context)
  text <- quality$text
  boxes <- data.frame(cx = text$centre_x_mm, cy = text$centre_y_mm,
    ow = text$width_mm, oh = text$height_mm, angle = text$angle_radians,
    xmin = text$xmin_mm, xmax = text$xmax_mm, ymin = text$ymin_mm, ymax = text$ymax_mm)
  to_mm <- function(x, y) {
    xy <- built$plot$coordinates$transform(data.frame(x = x, y = y),
      built$layout$panel_params[[1L]])
    data.frame(x = context$panel_bounds_mm[['xmin']] + xy$x * context$panel_mm[['width']],
      y = context$panel_bounds_mm[['ymin']] + xy$y * context$panel_mm[['height']])
  }
  extended <- detector(layout, text, boxes, to_mm)
  paths <- list(); contacts <- list()
  for (id in names(layout$layer_geometry)) {
    data <- layout$layer_geometry[[id]]$gene_poly
    if (!is.data.frame(data) || !nrow(data) || !'.component' %in% names(data)) next
    for (group in unique(data$group)) {
      part <- data[data$group == group & data$.component %in% c('point', 'boundary'), , drop = FALSE]
      if (!nrow(part)) next
      path_id <- paste(id, group, sep = '/')
      paths[[length(paths) + 1L]] <- data.frame(path_id = path_id,
        component = part$.component[1L], vertices = nrow(part))
      if (nrow(part) < 2L) next
      xy <- to_mm(part$x, part$y)
      for (j in seq_len(nrow(xy) - 1L)) {
        segment <- data.frame(x0 = xy$x[j], y0 = xy$y[j],
          x1 = xy$x[j + 1L], y1 = xy$y[j + 1L])
        for (i in seq_len(nrow(boxes))) if (ggchord_inner_leader_hits(segment, boxes[i, , drop = FALSE]))
          contacts[[length(contacts) + 1L]] <- data.frame(path_id = path_id,
            label_id = text$label_id[i], source_row = part$source_row[1L])
      }
    }
  }
  stopifnot(identical(original, plot$ggchord$data), identical(context$status, 'measured'))
  list(production = quality$summary, measurement = quality$measurement,
    retained_png = if (point_mode) file else NULL,
    production_skipped = quality$feature_parts_unchecked,
    extended_polygon_count = length(unique(extended$polygons$polygon_id)),
    extended_contacts = extended$contacts, remaining_skipped = extended$skipped,
    observed_paths = ggchord_rbind_fill(paths) %||% data.frame(),
    path_centreline_contacts = unique(ggchord_rbind_fill(contacts) %||% data.frame()))
}
shapes <- if (point_mode) 'point_records' else c('preset', 'capped_line', 'lollipop')
cases <- setNames(lapply(shapes, observe), shapes)
audit_assert_source_snapshot(audit_session)
jsonlite::write_json(list(status = 'observation completed; no production override',
  scope = if (point_mode) 'Existing optional package gridss.vcf point records and original header lengths; Collinear 10x8 inch PNG; IDs displayed verbatim; no clinical interpretation' else
    'Unchanged pUC19 annotations, three explicitly selected display geometries, 10x8 inch PNG',
  runtime_sha256 = as.list(audit_session$engine), input_sha256 = as.list(audit_session$inputs),
  source_records = if (point_mode) features[c('accver', 'start', 'end', 'feature_id')] else NULL,
  source_sequence_lengths = if (point_mode) sequence else NULL,
  cases = cases, not_checked = c('stroke thickness, dash gaps, caps and joins',
    'visibility/occlusion, glyph outlines, arbitrary custom grobs',
    'path contact source association and release device/platform matrix')),
  if (point_mode) 'tools/long-term-audit/v0.27.0-variant-point-contour-observations.json' else
    'tools/long-term-audit/v0.27.0-feature-contour-observations.json',
  pretty = TRUE, auto_unbox = TRUE, na = 'null')
cat('Feature contour observations completed; no production change\n')
