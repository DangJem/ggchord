# Matching text owners may request different sizes, but their shared feature
# geometry must use one measurement recipe, including in dependency subplots.
ggchord_resolve_feature_geometry_context <- function(plot) {
  layers <- plot$layers
  collinear <- isTRUE(plot$coordinates$ggchord_collinear)
  circular <- isTRUE(plot$coordinates$ggchord_circular)
  features <- which(vapply(layers, function(l)
    identical(l$ggchord_params$type, "gene"), logical(1)))
  labels <- which(vapply(layers, function(l)
    isTRUE(l$ggchord_params$type %in% c("gene_label", "gene_label_repel")) &&
      (collinear || isTRUE(l$ggchord_params$is_feature_label)), logical(1)))
  if (!length(labels)) return(plot)
  fallback <- ggchord_annotation_fallback(plot$ggchord$data)
  sizes <- vector("list", length(layers))
  for (label_row in labels) {
    layer <- layers[[label_row]]
    p <- layer$ggchord_params
    anchor <- p$anchor_layer
    candidates <- if (isTRUE(p$is_feature_label)) features else
      Filter(function(i) !isTRUE(layers[[i]]$ggchord_params$is_feature), features)
    if (!is.null(anchor)) {
      if (!is.character(anchor) || length(anchor) != 1L || is.na(anchor) || !nzchar(anchor))
        ggchord_stop("anchor_layer must name one feature layer")
      candidates <- Filter(function(i) identical(layers[[i]]$ggchord_layer_id, anchor), candidates)
      if (!length(candidates)) ggchord_stop("Unknown feature anchor_layer '", anchor, "'")
    }
    independent <- isTRUE(p$position_explicit %||% p$position_supplied)
    if (independent && !circular) next
    input <- ggchord_resolve_layer_input(layer, fallback)
    if (!NROW(input)) next
    found <- list()
    for (i in candidates) {
      feature <- layers[[i]]
      if (independent && !identical(p$feature_position, feature$ggchord_params$feature_position)) next
      ref <- ggchord_resolve_layer_input(feature, fallback)
      label_source <- attr(input, "ggchord_source")
      geom_source <- attr(ref, "ggchord_source")
      if (!is.null(label_source) && !is.null(geom_source) &&
          !identical(label_source$table, geom_source$table)) next
      shared <- identical(layer$ggchord_input_data, feature$ggchord_input_data) ||
        (!is.null(label_source) && identical(label_source, geom_source))
      if (independent && (!shared || NROW(input) != NROW(ref))) next
      rows <- ggchord_match_feature_rows(input, ref, shared, strict = !is.null(anchor))
      if (!is.null(rows)) found[[feature$ggchord_layer_id]] <- list(index = i, rows = rows)
    }
    if (independent && length(found) != 1L) next
    if (length(found) > 1L) ggchord_stop("Multiple geometry layers match these labels: ",
      paste(names(found), collapse = ", "), ". Select one with anchor_layer.")
    if (!length(found)) {
      if (!is.null(anchor)) ggchord_stop("Could not match labels to anchor_layer '", anchor,
        "'; check feature_id, source interval and displayed accession")
      next
    }
    owner <- found[[1L]]
    copy <- reconstruct_layer(layer, layer$data, mapping = layer$mapping)
    copy$ggchord_params$matched_feature <- list(layer_id = names(found)[1L], rows = owner$rows)
    plot$layers[[label_row]] <- copy
    if (circular && isTRUE(p$is_feature_label)) {
      size <- p$gene_label_size %||% ggchord_theme_text_size(plot, "ggchord.feature.label", 2.5)
      sizes[[owner$index]] <- c(sizes[[owner$index]], size)
    }
  }
  for (row in features) if (length(sizes[[row]])) {
    # Clone before setting internal context; caller layer objects stay owned
    # by their plots and no style or mapping is copied to another layer.
    layer <- reconstruct_layer(layers[[row]], layers[[row]]$data,
      mapping = layers[[row]]$mapping)
    layer$ggchord_params <- layers[[row]]$ggchord_params
    layer$ggchord_params$feature_layout_label_size <- max(sizes[[row]])
    layer$ggchord_params$feature_layout_has_labels <- TRUE
    plot$layers[[row]] <- layer
  }
  plot
}
