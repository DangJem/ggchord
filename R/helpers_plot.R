# helpers_plot.R - core plotting functions (used in v0.2.0, deprecated since v0.3.0)
# After the layered architecture refactor in v0.3.0, rendering logic was distributed across the geom_* layers
# This file is kept for reference and may be removed in a future version

# [Deprecated] Core plotting function from v0.2.0, called directly by ggchord()
#
# chordPlotFunc <- function(allRibbon, ribbon_alpha, ribbon_color_scheme, ribbon_colors,
#     show_legend, gene_polys, gene_pal, gene_color_scheme, final_gene_order,
#     seqArcs, gene_arrows, gene_label_show, gene_label_size, show_axis,
#     axisLines, axisTicks, axisLabelOrientation, seq_colors, seq_labels,
#     seqs, extremes, panel_margin, title) {
#   ... (no longer called in the new architecture; rendering logic is distributed across geom-seq.R / geom-ribbon.R /
#       geom-gene.R / geom-axis.R)
# }
