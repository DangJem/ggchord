# coord-chord.R — 轻量弦图坐标系统
# 使用 ggplot2 内置的 coord_fixed
# print.ggchord 阶段会动态替换 coord，此函数仅提供初始占位

#' 弦图坐标系统
#'
#' 轻量 Coord。在 \code{ggchord()} 中创建占位坐标，
#' 打印阶段根据布局计算的实际极值替换为合适的范围。
#'
#' @param layout 弦图布局对象（由 ggchord() 内部传入，可 NULL）
#'
#' @return Coord 对象，用于 ggplot2 的 \code{+} 叠加
#' @export
coord_chord <- function(layout = NULL) {
  # 初始占位范围，打印时会被 print.ggchord 替换
  coord_fixed(
    ratio = 1,
    xlim  = c(-5, 5),
    ylim  = c(-5, 5),
    clip  = "off"
  )
}
