# layout.R - core chord layout computation
# Pre-computes Cartesian (x, y) coordinates for all geometric elements from standardized parameters,
# for direct use by the geom_* layers

#' Compute the chord layout
#'
#' Pre-computes the coordinates of all geometric elements (sequence arcs, ribbons, gene arrows, axes, etc.)
#' into Cartesian (x, y) coordinates and stores them in a layout list.
#'
#' @param seqs Vector of sequence IDs (order already processed)
#' @param lens Named vector of sequence lengths (names = accver)
#' @param seq_labels Named vector of sequence labels
#' @param seqRadius Named vector of sequence radii
#' @param seq_curvature Named vector of sequence curvatures
#' @param orientation Named vector of sequence orientations (1 or -1)
#' @param seq_gap Named vector of sequence gap proportions
#' @param ribbonGap Named vector of ribbon gaps
#' @param ribbon_gap_auto Whether ribbon endpoints may move closer to sequence
#'   arcs when no gene or feature geometry occupies the local interval.
#' @param link_avoid Automatic link spacing mode.
#' @param ribbon_obstacles Optional normalized gene/feature obstacle table.
#' @param ribbon_data Alignment data (already validated)
#' @param gene_data Gene data (already validated)
#' @param draw_gene_geometry Whether gene/feature polygons should be generated.
#' @param gene_label_layout Character, default "radial". Deterministic
#'   automatic layout: "radial", "auto", or "arc".
#' @param gene_label_side Character, default "auto". Which side of the arc the
#'   labels sit on: "auto" (strand-based placement), "inside" (toward the chord
#'   center) or "outside" (away from the center, avoiding ribbon overlap).
#' @param gene_label_segment_linetype Character or numeric, default "auto".
#'   Leader-line linetype; "auto" uses solid lines except for labels moved to
#'   the other side of their arc, which use dashed lines.
#' @param gene_label_segment_overlap Character: "fade", "clip", or "show".
#' @param gene_label_segment_overlap_alpha Relative opacity for covered leader
#'   portions when `gene_label_segment_overlap = "fade"`.
#' @param gene_label_orientation Fixed-label text orientation: "radial",
#'   "tangent", or "horizontal". Automatic repel layouts manage their own
#'   text orientation.
#' @param gene_label_overlap Fixed-label collision policy: "hide", "nudge",
#'   or "allow".
#' @param seq_label_orientation Character, default "arc". Sequence label text
#'   orientation: "arc" (rotated along the arc, kept readable) or "horizontal"
#'   (all labels horizontal, extending away from the chord center).
#' @param seq_label_hjust Optional named vector or NULL, default NULL. Per-seq
#'   horizontal justification; NULL uses 0.5 (arc mode) or a side-based value
#'   (horizontal mode).
#' @param seq_label_vjust Optional named vector or NULL, default NULL. Per-seq
#'   vertical justification; NULL uses 0.5.
#' @param rotation Global rotation angle (degrees)
#' @param debug Whether to output debug information
#'
#' @return A chord layout list
#' @noRd
compute_chord_layout <- function(
    seqs, lens, seq_labels, seq_colors,
    seqRadius, seq_curvature, orientation, seq_gap,
    # Ribbon parameters
    ribbon_data = NULL, ribbonGap,
    ribbon_gap_auto = FALSE, ribbon_obstacles = NULL, link_avoid = "none",
    ribbon_color_scheme, ribbon_colors, ribbon_alpha,
    ribbon_color_by = NULL,
    ribbon_color_limits = NULL,
    ribbon_color_breaks = NULL,
    ribbon_color_name = NULL,
    ribbon_alpha_by = NULL,
    ribbon_alpha_range = c(0.15, 0.9),
    ribbon_outline_by = NULL,
    ribbon_outline_colors = NULL,
    ribbon_linetype_by = NULL,
    ribbon_linetypes = NULL,
    ribbon_direction = "none",
    ribbon_direction_colors = c(same = "black", reverse = "grey50"),
    ribbon_direction_linetypes = c(same = "solid", reverse = "dashed"),
    ribbon_direction_alpha = c(same = 1, reverse = 0.45),
    ribbon_ctrl_point,
    region_data = NULL,
    region_fill = "#F59E0B",
    region_color = "#B45309",
    region_alpha = 0.25,
    region_width = 0.08,
    region_offset = 0,
    region_side = "inside",
    # Gene parameters
    gene_data = NULL, draw_gene_geometry = TRUE,
    geneWidth,
    arrow_head_length = 0.04, arrow_head_width = 1,
    arrow_head_style = "shouldered",
    short_feature = "auto", circular = FALSE,
    geneLabelRadialOffset, geneLabelCircumOffset,
    geneLabelCircumLimit, geneLabelRotation,
    gene_label_show, gene_label_size,
    gene_label_family = "", gene_label_fontface = 1,
    gene_label_lineheight = 1.2,
    gene_label_wrap = NULL,
    gene_label_fit = "wrap",
    gene_label_max_lines = 2L,
    gene_label_orientation = "horizontal",
    gene_label_overlap = "hide",
    gene_label_repel_layer = FALSE,
    gene_label_repel_max_overlaps = Inf,
    gene_label_layout = "radial",
    feature_label_external = TRUE,
    gene_label_side = "auto",
    gene_label_segment_overlap = "fade",
    gene_label_segment_overlap_alpha = 0.18,
    gene_label_segment_linetype = "auto",
    gene_color_scheme, gene_colors, gene_order,
    feature_shape_pal = NULL, feature_shape_order = NULL,
    # Sequence label parameters
    seq_label_text = NULL, seq_label_radius = NULL,
    seq_label_rotation = NULL, seq_label_size = NULL,
    seq_label_orientation = "arc",
    seq_label_hjust = NULL, seq_label_vjust = NULL,
    # Axis parameters
    axisGap, axisMaj, axisMajLen, axisMin, axisMinLen,
    labelSize, labelOffset, axisLabelOrientation,
    axis_breaks = NULL, axis_minor_breaks = NULL, axis_labels = NULL,
    axis_label_hide_overlaps = FALSE,
    show_axis,
    # Global parameters
    rotation, debug = FALSE,
    geometry_cache = NULL
) {
  # Each focused module contributes one delayed expression. Evaluating them
  # in this frame preserves the original local state and layout contract.
  context <- environment()
  eval(ggchord_layout_sequence_step, envir = context)
  eval(ggchord_layout_axis_step, envir = context)
  eval(ggchord_layout_ribbon_step, envir = context)
  eval(ggchord_layout_annotation_step, envir = context)
  eval(ggchord_layout_transform_step, envir = context)
  eval(ggchord_layout_label_step, envir = context)
  eval(ggchord_layout_finalize_step, envir = context)
}
