# Axis-label pruning, plot bounds and the final layout contract.
ggchord_layout_finalize_step <- function(context) evalq({
  # ====================================================================
  # Step 8c: optionally hide axis labels that overlap other elements
  # ====================================================================
  if (isTRUE(axis_label_hide_overlaps) && nrow(axis_ticks) > 0 &&
      any(!is.na(axis_ticks$label))) {
    content <- ggchord_repel_points(seq_arcs, gene_polys,
                                    data.frame(x = numeric(0), y = numeric(0)),
                                    data.frame(), show_axis = FALSE)
    axis_ticks <- ggchord_hide_text_overlaps(
      axis_ticks, content,
      units_per_inch = text_units_per_inch
    )
  }

  # ====================================================================
  # Step 9: compute plot extremes
  # ====================================================================
  extremes <- get_plot_extremes(
    allRibbon = ribbon_polys,
    seqArcs = seq_arcs,
    axisLines = axis_lines,
    axisTicks = axis_ticks,
    gene_polys = gene_polys,
    gene_arrows = gene_labels,
    seq_labels = seq_labels_df,
    show_axis = show_axis
  )

  # ====================================================================
  # Step 10: assemble and return the layout object
  # ====================================================================
  layout <- list(
    sequence_reference = list(refs = seq_refs, starts = starts, ends = ends,
      lens = lens, orientation = orientation, radius = seqRadius),
    # Geometric data
    seq_arcs       = seq_arcs,
    ribbon_polys   = ribbon_polys,
    region_polys   = region_polys,
    gene_polys     = gene_polys,
    gene_labels    = gene_labels,
    gene_label_segments = gene_label_segments,
    gene_label_clip_units = gene_label_clip_units,
    text_units_per_inch = text_units_per_inch,
    gene_label_layout = gene_label_layout,
    feature_shape_pal = feature_shape_pal,
    feature_shape_order = feature_shape_order,
    seq_labels_df  = seq_labels_df,
    axis_lines     = axis_lines,
    axis_ticks     = axis_ticks,

    # Extremes
    extremes       = extremes,

    # Colors and labels
    seq_colors     = seq_colors,
    seq_labels     = seq_labels,
    seqs           = seqs,
    seqRadius      = seqRadius,

    # Ribbon-related
    ribbon_color_scheme = ribbon_color_scheme,
    ribbon_colors  = ribbon_colors,
    ribbon_alpha   = ribbon_alpha,
    ribbon_color_by = ribbon_color_by,
    ribbon_color_limits = ribbon_color_limits,
    ribbon_color_breaks = ribbon_color_breaks,
    ribbon_color_name = ribbon_color_name,
    ribbon_use_outline = ribbon_use_outline,
    ribbon_use_linetype = ribbon_use_linetype,

    # Gene-related
    gene_pal           = gene_pal,
    gene_color_scheme  = gene_color_scheme,
    final_gene_order   = final_gene_order,
    gene_label_show    = gene_label_show,
    gene_label_size    = gene_label_size,

    # Axis-related
    show_axis           = show_axis,
    axisLabelOrientation = axisLabelOrientation,

    # Metadata
    rotation        = rotation,
    n_sequences     = n
  )

  class(layout) <- "chord_layout"
  layout
}, envir = context)
