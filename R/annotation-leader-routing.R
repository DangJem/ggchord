# Protect static inner/axis text after the shared exterior queue has settled.
# Positions and real feature anchors are unchanged; only collided leader pieces
# can acquire a short deterministic polyline detour.
ggchord_annotation_text_obstacles <- function(registry, layout, include_axis = TRUE) {
  all <- list()
  unit <- layout$text_units_per_inch %||% .35
  add <- function(text, owner) {
    if (!nrow(text)) return()
    text <- text[!is.na(text$label) & nzchar(text$label), , drop = FALSE]
    if (!nrow(text)) return()
    text$text <- text$label; text$text_x <- text$x; text$text_y <- text$y
    text$text_angle <- text$angle %||% 0
    boxes <- ggchord_text_boxes(text, units_per_inch = unit, box_padding = 0)
    boxes$owner <- owner; boxes$source_row <- text$source_row %||% rep(NA_integer_, nrow(text))
    object <- if ('.arc_parent_group' %in% names(text)) as.character(text$.arc_parent_group) else
      if ('source_row' %in% names(text)) as.character(text$source_row) else as.character(seq_len(nrow(text)))
    boxes$text_object <- paste(owner, object, sep = '/')
    all[[length(all)+1L]] <<- boxes
  }
  for (owner in names(registry)) {
    d <- registry[[owner]]$gene_label_repel
    if (is.data.frame(d) && nrow(d)) {
      keep <- d$.component %in% c('text', 'arc_text')
      if ('.draw_as_arc' %in% names(d)) keep <- keep &
        !(d$.component %in% 'text' & d$.draw_as_arc %in% TRUE)
      add(d[keep, , drop = FALSE], owner)
    }
    d <- registry[[owner]]$restriction_site
    if (is.data.frame(d) && nrow(d)) add(d[d$restriction_component %in% 'label', , drop = FALSE], owner)
    d <- registry[[owner]]$seq_center_label
    if (is.data.frame(d) && nrow(d)) add(d, owner)
  }
  if (include_axis) {
    axis_layout <- layout
    axis_layout$layer_geometry <- registry
    axis <- ggchord_axis_geometry(axis_layout)
    if (nrow(axis)) add(axis[axis$.component %in% 'text', , drop = FALSE], '.axis')
  }
  ggchord_rbind_fill(all) %||% data.frame()
}

# The exterior queue already provides direct or one-shoulder routes.
# Preserve those routes instead of turning collisions into multi-bend graph
# detours. Text positions and fixed source roots remain coordinate-owned.
ggchord_detour_annotation_leaders <- function(registry, layout) {
  if (!isTRUE(layout$circular)) return(registry)
  ggchord_measure_external_leader_crossings(registry)
}
