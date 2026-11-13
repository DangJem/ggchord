# Matching text owners may request different sizes, but their shared feature
# geometry must use one measurement recipe, including in dependency subplots.
ggchord_resolve_feature_geometry_context <- function(plot) {
  if (!isTRUE(plot$coordinates$ggchord_circular)) return(plot)
  layers <- plot$layers
  features <- which(vapply(layers, function(layer)
    identical(layer$ggchord_params$type, "gene") &&
      isTRUE(layer$ggchord_params$is_feature), logical(1)))
  labels <- which(vapply(layers, function(layer)
    layer$ggchord_params$type %in% c("gene_label", "gene_label_repel") &&
      isTRUE(layer$ggchord_params$is_feature_label), logical(1)))
  if (!length(features) || !length(labels)) return(plot)
  fallback <- ggchord_annotation_fallback(plot$ggchord$data)
  sizes <- vector("list", length(layers))
  for (label_row in labels) {
    label <- layers[[label_row]]
    input <- ggchord_resolve_layer_input(label, fallback)
    for (feature_row in features) {
      feature <- layers[[feature_row]]
      reference <- ggchord_resolve_layer_input(feature, fallback)
      compare <- intersect(c("accver", "start", "end", "strand", "anno", "label"),
        intersect(names(input), names(reference)))
      same_input <- nrow(input) == nrow(reference) &&
        all(c("accver", "start", "end", "strand") %in% compare) &&
        identical(input[compare], reference[compare])
      same_position <- !isTRUE(label$ggchord_params$position_supplied) ||
        identical(label$ggchord_params$feature_position,
          feature$ggchord_params$feature_position)
      if (!same_input || !same_position) next
      size <- label$ggchord_params$gene_label_size %||%
        ggchord_theme_text_size(plot, "ggchord.feature.label", 2.5)
      sizes[[feature_row]] <- c(sizes[[feature_row]], size)
      break
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
