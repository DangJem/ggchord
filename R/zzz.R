# zzz.R — 包级环境和基础设施
# ggchord 包级环境，用于在 ggchord() 和各 geom 之间传递数据、参数和布局

#' 包级环境
#'
#' 内部环境，存储弦图的原始数据、各 geom 的参数、以及延迟计算的布局结果。
#'
#' @keywords internal
.chord_env <- new.env(parent = emptyenv())

# ====================================================================
# 数据存储
# ====================================================================

#' 设置原始数据到包级环境
#' @keywords internal
set_chord_data <- function(seq_data, ribbon_data, gene_data) {
  .chord_env$seq_data    <- seq_data
  .chord_env$ribbon_data <- ribbon_data
  .chord_env$gene_data   <- gene_data
}

#' 获取原始数据
#' @keywords internal
get_chord_data <- function() {
  list(
    seq_data    = .chord_env$seq_data,
    ribbon_data = .chord_env$ribbon_data,
    gene_data   = .chord_env$gene_data
  )
}

# ====================================================================
# 参数存储（各 geom 在叠加阶段将参数存入此环境）
# ====================================================================

# 全局参数
set_global_params <- function(params) {
  .chord_env$global_params <- params
}
get_global_params <- function() {
  .chord_env$global_params
}

# 序列参数（来自 geom_seq）
set_seq_params <- function(params) {
  .chord_env$seq_params <- params
}
get_seq_params <- function() {
  .chord_env$seq_params
}

# Ribbon 参数（来自 geom_ribbon）
set_ribbon_params <- function(params) {
  .chord_env$ribbon_params <- params
}
get_ribbon_params <- function() {
  .chord_env$ribbon_params
}

# 基因参数（来自 geom_gene）
set_gene_params <- function(params) {
  .chord_env$gene_params <- params
}
get_gene_params <- function() {
  .chord_env$gene_params
}

# 坐标轴参数（来自 geom_axis）
set_axis_params <- function(params) {
  .chord_env$axis_params <- params
}
get_axis_params <- function() {
  .chord_env$axis_params
}

# ====================================================================
# 布局存储（延迟计算，print 时执行）
# ====================================================================

#' 设置弦图布局到包级环境
#' @keywords internal
set_chord_layout <- function(layout) {
  .chord_env$layout <- layout
}

#' 从包级环境获取弦图布局
#' @keywords internal
get_chord_layout <- function() {
  layout <- .chord_env$layout
  if (is.null(layout)) {
    stop(
      "未找到弦图布局数据。请先调用 ggchord() 创建布局。",
      call. = FALSE
    )
  }
  layout
}

# ====================================================================
# 环境清理
# ====================================================================

#' 清空包级环境（用于重置状态）
#' @keywords internal
clear_chord_env <- function() {
  rm(list = ls(.chord_env, all.names = TRUE), envir = .chord_env)
}

# ====================================================================
# 辅助运算符
# ====================================================================

#' NULL 合并运算符
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
