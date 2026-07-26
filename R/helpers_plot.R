# helpers_plot.R — 核心绘图函数（v0.2.0 使用，v0.3.0 后已废弃）
# v0.3.0 重构为分层架构后，渲染逻辑已分散到各 geom_* 图层中
# 此文件保留作为参考，未来版本可能删除

# [已废弃] v0.2.0 的核心绘图函数，由 ggchord() 直接调用
#
# chordPlotFunc <- function(allRibbon, ribbon_alpha, ribbon_color_scheme, ribbon_colors,
#     show_legend, gene_polys, gene_pal, gene_color_scheme, final_gene_order,
#     seqArcs, gene_arrows, gene_label_show, gene_label_size, show_axis,
#     axisLines, axisTicks, axisLabelOrientation, seq_colors, seq_labels,
#     seqs, extremes, panel_margin, title) {
#   ...（在新架构中不再被调用，渲染逻辑已分散到 geom-seq.R / geom-ribbon.R /
#       geom-gene.R / geom-axis.R 中）
# }
