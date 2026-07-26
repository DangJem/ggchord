# geom-axis.R — 坐标轴图层
# 从包级环境获取预计算的轴线、刻度线和标签数据
# 坐标轴参数在此层指定，存入环境供 print 时使用

#' 添加坐标轴图层
#'
#' 绘制弦图中每条序列的坐标轴（含轴线、主/次刻度线和标签）。
#' 坐标轴参数（间距、刻度数量/长度、标签大小/方向等）在此指定。
#'
#' @param mapping 默认 NULL（使用预计算数据）
#' @param data 默认 NULL（从布局中自动获取）
#' @param show_axis 逻辑型。是否显示坐标轴，默认 TRUE
#' @param axis_gap 数值型/向量，可选。序列与坐标轴间距，默认 0.04
#' @param axis_tick_major_number 整数型/向量，可选。主刻度数量，默认 5
#' @param axis_tick_major_length 数值型/向量，可选。主刻度长度比例，默认 0.02
#' @param axis_tick_minor_number 整数型/向量，可选。次刻度数量，默认 4
#' @param axis_tick_minor_length 数值型/向量，可选。次刻度长度比例，默认 0.01
#' @param axis_label_size 数值型/向量，可选。刻度标签字号，默认 3
#' @param axis_label_offset 数值型/向量，可选。标签偏移比例，默认 1.5
#' @param axis_label_orientation 字符/数值/向量，可选。标签方向，默认 "horizontal"
#' @param show_legend 是否显示图例，默认 FALSE（坐标轴不参与图例）
#' @param ... 传递给 geom_path/geom_segment/geom_text 的其他参数
#'
#' @return ggplot2 layer 列表
#' @export
geom_axis <- function(mapping = NULL, data = NULL,
                      show_axis = NULL,
                      axis_gap = NULL,
                      axis_tick_major_number = NULL,
                      axis_tick_major_length = NULL,
                      axis_tick_minor_number = NULL,
                      axis_tick_minor_length = NULL,
                      axis_label_size = NULL,
                      axis_label_offset = NULL,
                      axis_label_orientation = NULL,
                      show_legend = FALSE,
                      ...) {
  set_axis_params(list(
    show_axis                = show_axis,
    axis_gap                 = axis_gap,
    axis_tick_major_number   = axis_tick_major_number,
    axis_tick_major_length   = axis_tick_major_length,
    axis_tick_minor_number   = axis_tick_minor_number,
    axis_tick_minor_length   = axis_tick_minor_length,
    axis_label_size          = axis_label_size,
    axis_label_offset        = axis_label_offset,
    axis_label_orientation   = axis_label_orientation
  ))

  empty_id <- data.frame(x = numeric(0), y = numeric(0),
                         seq_id = character(0))
  empty_seg <- data.frame(x0 = numeric(0), y0 = numeric(0),
                          x1 = numeric(0), y1 = numeric(0),
                          label = character(0), label_x = numeric(0),
                          label_y = numeric(0), size = numeric(0),
                          seq_id = character(0))

  list(
    geom_path(data = empty_id,
              mapping = aes(x = x, y = y, group = seq_id),
              color = "black", linewidth = 0.3,
              inherit.aes = FALSE, show.legend = show_legend, ...),
    geom_segment(data = empty_seg,
                 mapping = aes(x = x0, y = y0,
                               xend = x1, yend = y1),
                 color = "black", linewidth = 0.3,
                 inherit.aes = FALSE, show.legend = show_legend, ...),
    geom_text(data = empty_seg[integer(0), ],
              mapping = aes(x = label_x, y = label_y,
                            label = label, size = size),
              inherit.aes = FALSE, color = "black",
              angle = 0, show.legend = show_legend, ...)
  )
}
