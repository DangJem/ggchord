# geom-ribbon.R — 比对 ribbon 图层
# 从包级环境获取预计算的 ribbon 多边形数据，用 geom_polygon 渲染
# Ribbon 参数在此层指定，存入环境供 print 时使用

#' 添加比对 ribbon 图层
#'
#' 绘制 BLAST 比对结果对应的彩色 ribbon。配色方案和间距参数在此指定。
#'
#' @param mapping 默认 NULL（使用预计算数据）
#' @param data 默认 NULL（从布局中自动获取）
#' @param ribbon_color_scheme 字符型。配色方案 "pident"、"query"、"single"，默认 "pident"
#' @param ribbon_colors 颜色向量，可选。ribbon 颜色参数
#' @param ribbon_alpha 数值型 (0-1)。ribbon 透明度，默认 0.35
#' @param ribbon_ctrl_point 向量/列表，可选。贝塞尔控制点，默认 c(0,0)
#' @param ribbon_gap 数值型/向量，可选。序列与 ribbon 间距，默认 0.15
#' @param alpha ribbon 透明度（可覆盖 ribbon_alpha），默认使用布局中的值
#' @param show_legend 是否显示图例，默认 TRUE
#' @param ... 传递给 \code{geom_polygon()} 的其他参数
#'
#' @return ggplot2 layer 列表
#' @export
geom_ribbon <- function(mapping = NULL, data = NULL,
                        ribbon_color_scheme = NULL,
                        ribbon_colors = NULL,
                        ribbon_alpha = NULL,
                        ribbon_ctrl_point = NULL,
                        ribbon_gap = NULL,
                        alpha = NULL,
                        show_legend = TRUE,
                        ...) {
  set_ribbon_params(list(
    ribbon_color_scheme = ribbon_color_scheme,
    ribbon_colors       = ribbon_colors,
    ribbon_alpha        = alpha %||% ribbon_alpha,
    ribbon_ctrl_point   = ribbon_ctrl_point,
    ribbon_gap          = ribbon_gap
  ))

  # 占位数据
  empty_polys <- data.frame(
    x = numeric(0), y = numeric(0),
    group = integer(0),
    # ggnewscale builds the plot as soon as a new fill scale is added.
    # Keep both possible mapped columns on the placeholder so that this
    # preliminary build succeeds before print.ggchord() injects real data.
    pident = numeric(0),
    fill = character(0),
    alpha = numeric(0),
    stringsAsFactors = FALSE
  )

  # 检测 ribbon_color_scheme 来决定 fill 映射变量
  # 但 print 时会注入真实数据并替换 fill 列
  scheme <- ribbon_color_scheme %||% "pident"

  if (scheme == "pident") {
    fill_aes <- aes(x = x, y = y, group = group, fill = pident, alpha = alpha)
  } else {
    fill_aes <- aes(x = x, y = y, group = group, fill = fill, alpha = alpha)
  }

  list(
    ggplot2::layer(
      data        = empty_polys,
      mapping     = fill_aes,
      stat        = "identity",
      geom        = GeomPolygon,
      position    = "identity",
      show.legend = show_legend,
      inherit.aes = FALSE,
      params      = list(...)
    )
  )
}
