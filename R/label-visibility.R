# A row is reported even when it produced no visible text. Geometry is retained.
ggchord_label_visibility <- function(layout) {
  rows <- list()
  for (owner in names(layout$layer_inputs)) for (component in intersect(
      names(layout$layer_inputs[[owner]]),c('gene_text','gene_text_repel','gene_label_repel','restriction_site'))) {
    input <- layout$layer_inputs[[owner]][[component]]
    if (is.null(input) || !nrow(input)) next
    d <- layout$layer_geometry[[owner]][[component]]
    if (is.null(d)) d <- data.frame()
    text <- if (nrow(d)) {
      keep <- if (component=='restriction_site') d$restriction_component %in% 'label' else
        if ('.component' %in% names(d)) d$.component %in% c('text','arc_text') else rep(TRUE,nrow(d))
      d[keep & !is.na(d$label) & nzchar(d$label),,drop=FALSE]
    } else d
    for (i in seq_len(nrow(input))) {
      at <- if(nrow(text) && "source_rows" %in% names(text)) which(vapply(text$source_rows,function(rows) i %in% rows,logical(1))) else
        if(nrow(text)) which(text$source_row==i) else integer()
      visible <- length(at)>0
      source_text <- as.character((input$label %||% input$anno %||% input$enzyme %||% rep(NA_character_,nrow(input)))[i])
      selected <- if (visible) text[at[1],,drop=FALSE] else data.frame()
      scalar <- function(field,default=NA_real_) if (visible && field %in% names(selected)) selected[[field]][1] else default
      rows[[length(rows)+1L]] <- data.frame(layer_id=owner,component=component,source_row=i,
        accver=as.character(input$accver[i]),source_text=source_text,visible=visible,
        reason=if(visible && 'source_rows' %in% names(selected) && length(selected$source_rows[[1]])>1)'represented_by_shared_label' else if(visible)'visible' else if(!is.na(source_text) && !nzchar(source_text))'empty_label' else 'not_placed_by_layer_policy',
        actual_label=as.character(scalar('label',NA_character_)),
        displayed_source_row=scalar('source_row'),
        text_x=scalar('text_x',scalar('x')),text_y=scalar('text_y',scalar('y')),
        anchor_x=scalar('anchor_x'),anchor_y=scalar('anchor_y'),
        size=scalar('size'),side=scalar('outer_annotation_region',NA_character_),
        track=scalar('outer_track'),leader_crossings=scalar('leader_crossing_count'),stringsAsFactors=FALSE)
    }
  }
  ggchord_restore_instance_data(ggchord_rbind_fill(rows) %||% data.frame(),layout$display_instances)
}

# Shared diagnostics for the final measured layout, without rebuilding it.
ggchord_label_fit_diagnostics <- function(layout) {
  fits <- lapply(names(layout$layer_layouts), function(id) {
    # Dependency sub-layouts may carry compatibility label diagnostics even
    # when their owner only draws a backbone or feature. Export only the
    # actual text owner's diagnostics, including owners whose text was hidden.
    if (!any(names(layout$layer_geometry[[id]]) %in%
        c("gene_text", "gene_text_repel", "gene_label_repel"))) return(NULL)
    d <- layout$layer_layouts[[id]]$label_fit_diagnostics
    if (!is.data.frame(d) || !nrow(d)) return(NULL)
    resolved <- lapply(layout$layer_geometry[[id]], function(g) {
      if (!is.data.frame(g) || !nrow(g) || !"label_size_requested" %in% names(g)) return(NULL)
      if (".component" %in% names(g)) g <- g[g$.component %in% c("text", "arc_text"), , drop = FALSE]
      if (!nrow(g)) return(NULL)
      unique(g[intersect(c("source_row", "label_size_final", "label_fit_reason"), names(g))])
    })
    resolved <- ggchord_rbind_fill(resolved)
    if (is.data.frame(resolved) && nrow(resolved)) {
      at <- match(d$source_row, resolved$source_row)
      for (field in intersect(c("label_size_final", "label_fit_reason"), names(resolved)))
        d[[field]][!is.na(at)] <- resolved[[field]][at[!is.na(at)]]
      d$label_fit_reason[is.na(at)] <- "hidden_after_layout"
    } else d$label_fit_reason <- "hidden_after_layout"
    d$layer_id <- id
    d
  })
  ggchord_rbind_fill(fits) %||% data.frame()
}
