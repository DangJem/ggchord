# Conservative envelope contacts, with source evidence kept separate from shape.
ggchord_label_feature_quality <- function(layout, text, boxes, to_mm) {
  sources <- ggchord_annotation_source_maps(layout)
  identity <- function(layer, row, components) {
    if (!nrow(sources) || is.na(row)) return(NA_character_)
    at <- sources$layer_id == layer & sources$input_row == row & sources$component %in% components
    ids <- unique(sources$source_id[at])
    if (length(ids) == 1L) ids else NA_character_
  }
  polygons <- list(); skipped <- list(); contacts <- list()
  for (id in names(layout$layer_geometry)) {
    data <- layout$layer_geometry[[id]]$gene_poly
    if (!is.data.frame(data) || !nrow(data)) next
    for (group in unique(data$group)) {
      part <- data[data$group == group, , drop = FALSE]
      source <- unique(part$source_row)
      source <- if (length(source) == 1L) source else NA_integer_
      key <- paste(id, group, sep = '/')
      is_polygon <- if ('.component' %in% names(part)) is.na(part$.component) | part$.component %in% 'polygon' else rep(TRUE, nrow(part))
      if (any(!is_polygon)) skipped[[length(skipped)+1L]] <- data.frame(polygon_id = key,
        layer_id = id, source_row = source, reason = 'non_polygon_components')
      part <- part[is_polygon, , drop = FALSE]
      if (!nrow(part)) next
      if (nrow(part) < 3L || any(!is.finite(part$x) | !is.finite(part$y))) {
        skipped[[length(skipped)+1L]] <- data.frame(polygon_id = key, layer_id = id,
          source_row = source, reason = 'not_a_finite_polygon')
        next
      }
      xy <- to_mm(part$x, part$y)
      source_id <- identity(id, source, 'gene_poly')
      polygons[[length(polygons)+1L]] <- data.frame(polygon_id = key, layer_id = id,
        source_row = source, source_id = source_id, x_mm = xy$x, y_mm = xy$y)
      if (!nrow(text)) next
      at <- which(text$xmin_mm <= max(xy$x) & text$xmax_mm >= min(xy$x) &
        text$ymin_mm <= max(xy$y) & text$ymax_mm >= min(xy$y))
      for (i in at) if (ggchord_box_overlaps_feature(boxes[i, , drop = FALSE], xy, tol = 1e-8)) {
        label_source <- identity(text$layer_id[i], text$source_row[i],
          c('gene_text', 'gene_text_repel', 'gene_label_repel'))
        relation <- if (is.na(label_source) || is.na(source_id)) 'unresolved' else
          if (identical(label_source, source_id)) 'own_source' else
            if (grepl('^(feature|gene):', label_source) && grepl('^(feature|gene):', source_id))
              'other_source' else 'unresolved'
        contacts[[length(contacts)+1L]] <- data.frame(label_id = text$label_id[i],
          feature_id = paste(id, source, sep = '/'), label_source_id = label_source,
          feature_source_id = source_id, relation = relation)
      }
    }
  }
  list(contacts = unique(ggchord_rbind_fill(contacts) %||% data.frame(label_id = character(),
    feature_id = character(), label_source_id = character(), feature_source_id = character(), relation = character())),
    polygons = ggchord_rbind_fill(polygons) %||% data.frame(),
    skipped = ggchord_rbind_fill(skipped) %||% data.frame())
}
