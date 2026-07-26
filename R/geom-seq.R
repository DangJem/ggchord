# geom-seq.R — 序列弧线图层
# 从包级环境获取预计算的序列弧线数据，用 geom_path 渲染
# 序列布局参数在此层指定，存入环境供 print 时使用

#' 添加序列弧线图层
#'
#' 绘制弦图中表示序列的弧线（或直线，取决于曲率设置）。
#' 序列布局参数（顺序、方向、半径、曲率、颜色等）在此指定。
#'
#' @param mapping 默认 NULL（使用预计算数据）
#' @param data 默认 NULL（从布局中自动获取）
#' @param seq_order 字符型向量，可选。指定序列的绘制顺序
#' @param seq_labels 字符型向量/命名向量，可选。序列标签
#' @param seq_orientation 数值型 (1 或 -1)，可选。序列方向，默认 1
#' @param seq_gap 数值型，可选。序列间空白比例，默认 0.03
#' @param seq_radius 数值型 (> 0)，可选。序列弧线半径，默认 1.0
#' @param seq_curvature 数值型，可选。弧线曲率 (0=直线, 1=标准弧, >1=更弯)，默认 1.0
#' @param seq_colors 颜色向量/命名向量，可选。序列颜色
#' @param linewidth 弧线宽度，默认 1.2
#' @param show_legend 是否显示该图层的图例，默认 TRUE
#' @param ... 传递给 \code{geom_path()} 的其他参数
#'
#' @return ggplot2 layer 列表
#' @export
geom_seq <- function(mapping = NULL, data = NULL,
                     seq_order = NULL,
                     seq_labels = NULL,
                     seq_orientation = NULL,
                     seq_gap = NULL,
                     seq_radius = NULL,
                     seq_curvature = NULL,
                     seq_colors = NULL,
                     linewidth = 1.2,
                     show_legend = TRUE,
                     ...) {
  set_seq_params(list(
    seq_order      = seq_order,
    seq_labels     = seq_labels,
    seq_orientation = seq_orientation,
    seq_gap        = seq_gap,
    seq_radius     = seq_radius,
    seq_curvature  = seq_curvature,
    seq_colors     = seq_colors
  ))

  # 仅返回 layer（scale 由 print.ggchord 统一设置）
  list(
    ggplot2::layer(
      data        = data.frame(x = numeric(0), y = numeric(0),
                               seq_id = character(0)),
      mapping     = aes(x = x, y = y, group = seq_id, color = seq_id),
      stat        = "identity",
      geom        = GeomPath,
      position    = "identity",
      show.legend = show_legend,
      inherit.aes = FALSE,
      params      = list(
        linewidth = linewidth,
        arrow = grid::arrow(type = "closed", length = grid::unit(3, "mm")),
        ...
      )
    )
  )
}
