# Evidence from one measured panel. Detect contacts without moving/hiding data.
ggchord_label_quality <- function(layout, text, built, context) {
  segments <- list()
  coordinate <- built$plot$coordinates
  panel <- built$layout$panel_params[[1L]]
  to_mm <- function(x, y) {
    xy <- coordinate$transform(data.frame(x = x, y = y), panel)
    data.frame(x = context$panel_bounds_mm[['xmin']] + xy$x * context$panel_mm[['width']],
      y = context$panel_bounds_mm[['ymin']] + xy$y * context$panel_mm[['height']])
  }
  add <- function(id, component, source, part, x, y, xend, yend) {
    a <- to_mm(x, y); b <- to_mm(xend, yend)
    if (!nrow(a)) return()
    out <- data.frame(leader_id = paste(id, component, source, part, sep = '/'),
      layer_id = id, source_row = source, x0 = a$x, y0 = a$y, x1 = b$x, y1 = b$y)
    out <- out[stats::complete.cases(out[c('x0', 'y0', 'x1', 'y1')]), , drop = FALSE]
    out <- out[out$x0 != out$x1 | out$y0 != out$y1, , drop = FALSE]
    if (nrow(out)) segments[[length(segments)+1L]] <<- out
  }
  for (id in names(layout$layer_geometry)) {
    labels <- layout$layer_geometry[[id]]$gene_label_repel
    if (is.data.frame(labels) && nrow(labels)) {
      at <- which(labels$.component %in% 'segment')
      if (length(at)) {
        part <- labels$accver %||% rep('', nrow(labels))
        for (key in unique(paste(labels$source_row[at], part[at], sep = '\r'))) {
          index <- at[paste(labels$source_row[at], part[at], sep = '\r') == key]
          add(id, 'gene_label_repel', labels$source_row[index[1L]], part[index[1L]],
            labels$x[index], labels$y[index], labels$xend[index], labels$yend[index])
        }
      }
    }
    sites <- layout$layer_geometry[[id]]$restriction_site
    if (is.data.frame(sites) && nrow(sites)) {
      at <- which(sites$restriction_component %in% 'leader')
      for (source in unique(sites$source_row[at])) {
        index <- at[sites$source_row[at] == source]
        if (length(index) < 2L) next
        # Consecutive draw-order vertices, never a radial sort.
        add(id, 'restriction_site', source, '', sites$x[head(index, -1L)],
          sites$y[head(index, -1L)], sites$x[tail(index, -1L)], sites$y[tail(index, -1L)])
      }
    }
  }
  segments <- ggchord_rbind_fill(segments) %||% data.frame(leader_id = character(),
    layer_id = character(), source_row = integer(), x0 = numeric(), y0 = numeric(), x1 = numeric(), y1 = numeric())
  conflicts <- list()
  contact <- function(type, first, second) {
    conflicts[[length(conflicts)+1L]] <<- data.frame(type = type, first = first, second = second)
  }
  boxes <- data.frame()
  if (nrow(text)) {
    text$text_id <- paste0('text-', seq_len(nrow(text)))
    text$label_id <- paste(text$layer_id, text$text_object, sep = '/')
    boxes <- data.frame(cx = text$centre_x_mm, cy = text$centre_y_mm,
      ow = text$width_mm, oh = text$height_mm, angle = text$angle_radians,
      xmin = text$xmin_mm, xmax = text$xmax_mm, ymin = text$ymin_mm, ymax = text$ymax_mm)
    if (nrow(text) > 1L) for (i in seq_len(nrow(text)-1L)) {
      at <- seq.int(i+1L, nrow(text))
      # Glyphs or parts of one registered source label are intentional peers.
      same <- text$label_id[at] == text$label_id[i]
      at <- at[!same & text$xmin_mm[at] < text$xmax_mm[i] & text$xmax_mm[at] > text$xmin_mm[i] &
        text$ymin_mm[at] < text$ymax_mm[i] & text$ymax_mm[at] > text$ymin_mm[i]]
      if (length(at)) for (j in at[ggchord_oriented_box_overlaps(boxes[i, , drop = FALSE], boxes[at, , drop = FALSE])])
        contact('text_overlap', text$label_id[i], text$label_id[j])
    }
    for (i in seq_len(nrow(segments))) {
      at <- which(text$xmin_mm < max(segments$x0[i], segments$x1[i]) &
        text$xmax_mm > min(segments$x0[i], segments$x1[i]) &
        text$ymin_mm < max(segments$y0[i], segments$y1[i]) &
        text$ymax_mm > min(segments$y0[i], segments$y1[i]))
      for (j in at) if (ggchord_inner_leader_hits(segments[i, , drop = FALSE], boxes[j, , drop = FALSE]))
        contact('leader_through_text', segments$leader_id[i], text$label_id[j])
    }
  }
  features <- ggchord_label_feature_quality(layout, text, boxes, to_mm)
  if (nrow(features$contacts)) for (i in which(features$contacts$relation != 'own_source'))
    contact(if (features$contacts$relation[i] == 'other_source') 'text_over_feature' else
      'text_feature_unresolved', features$contacts$label_id[i], features$contacts$feature_id[i])
  ids <- unique(segments$leader_id)
  if (length(ids) > 1L) for (i in seq_len(length(ids)-1L)) for (j in seq.int(i+1L, length(ids))) {
    a <- segments[segments$leader_id == ids[i], , drop = FALSE]
    b <- segments[segments$leader_id == ids[j], , drop = FALSE]
    crosses <- any(vapply(seq_len(nrow(a)), function(k) any(vapply(seq_len(nrow(b)), function(m)
      ggchord_external_leader_crosses(
        data.frame(x1 = a$x0[k], y1 = a$y0[k], x2 = a$x1[k], y2 = a$y1[k]),
        data.frame(x1 = b$x0[m], y1 = b$y0[m], x2 = b$x1[m], y2 = b$y1[m])), logical(1))), logical(1)))
    if (crosses) contact('leader_crossing', ids[i], ids[j])
  }
  conflicts <- unique(ggchord_rbind_fill(conflicts) %||%
    data.frame(type = character(), first = character(), second = character()))
  lengths <- if (nrow(segments)) tapply(sqrt((segments$x1-segments$x0)^2 + (segments$y1-segments$y0)^2),
    segments$leader_id, sum) else numeric()
  visibility <- ggchord_label_visibility(layout)
  measured <- identical(context$status, 'measured')
  summary <- list(measured = measured, n_source_labels = nrow(visibility),
    n_visible_source_labels = if (nrow(visibility)) sum(visibility$visible) else 0L,
    n_hidden_source_labels = if (nrow(visibility)) sum(!visibility$visible) else 0L,
    n_text_boxes = nrow(text), n_text_objects = if (nrow(text)) length(unique(text$label_id)) else 0L,
    n_leaders = length(ids),
    n_text_overlaps = sum(conflicts$type == 'text_overlap'),
    n_leader_crossings = sum(conflicts$type == 'leader_crossing'),
    n_leaders_through_text = sum(conflicts$type == 'leader_through_text'),
    n_label_outside_canvas = if (nrow(text)) length(unique(text$label_id[text$xmin_mm < 0 |
      text$xmax_mm > context$canvas_mm[['width']] | text$ymin_mm < 0 |
      text$ymax_mm > context$canvas_mm[['height']]])) else 0L,
    mean_leader_length_mm = if (length(lengths)) mean(lengths) else NA_real_,
    max_leader_length_mm = if (length(lengths)) max(lengths) else NA_real_,
    n_feature_polygons = if (nrow(features$polygons)) length(unique(features$polygons$polygon_id)) else 0L,
    n_feature_parts_unchecked = nrow(features$skipped),
    n_text_over_feature = sum(features$contacts$relation == 'other_source'),
    n_text_feature_unresolved = sum(features$contacts$relation == 'unresolved'),
    n_text_on_own_feature = sum(features$contacts$relation == 'own_source'))
  if (!measured) for (key in setdiff(names(summary), 'measured')) summary[[key]] <- NA
  list(summary = summary, conflicts = conflicts, text = text, leaders = segments,
    feature_contacts = features$contacts, feature_polygons = features$polygons,
    feature_parts_unchecked = features$skipped,
    scope = 'Measured registered panel text boxes, annotation polylines and gene/feature polygon envelopes, including transparent polygons. Own-source label contacts are recorded separately; unproven source association remains unresolved. Strict leader interior contacts, not shared endpoints. No automatic layout change.',
    not_checked = c('non-polygon feature strokes/markers, stroke width and visibility/occlusion', 'unregistered custom grobs',
      'outer plot decoration/guide text', 'glyph outlines/ligatures or cross-application rendering'))
}
