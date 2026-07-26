# geom-gene.R — 基因箭头图层
# 从包级环境获取预计算的基因箭头多边形和标签数据
# 基因参数在此层指定，存入环境供 print 时使用

#' 添加基因箭头图层
#'
#' 在弦图上绘制基因注释箭头。基因布局参数（偏移、宽度、配色等）在此指定。
#' 当 \code{ggnewscale} 包可用时，自动使用 \code{new_scale_fill()} 确保
#' 与 ribbon 的 fill scale 独立。
#'
#' @param mapping 默认 NULL（使用预计算数据）
#' @param data 默认 NULL（从布局中自动获取）
#' @param gene_offset 数值型/向量/列表，可选。基因箭头径向偏移，默认 0.1
#' @param gene_width 数值型/向量/列表，可选。基因箭头宽度，默认 0.05
#' @param gene_color_scheme 字符型。"strand" 或 "manual"，默认 "strand"
#' @param gene_colors 颜色向量，可选。基因箭头填充色
#' @param gene_order 字符型向量，可选。基因在图例中的显示顺序
#' @param gene_label_show 逻辑型。是否显示基因标签，默认 FALSE
#' @param gene_label_size 数值型。标签字号，默认 2.5
#' @param gene_label_rotation 数值型/向量/列表，可选。标签旋转角度，默认 0
#' @param gene_label_radial_offset 数值型/向量/列表，可选。标签径向偏移，默认 0
#' @param gene_label_circum_offset 数值型/向量/列表，可选。标签周向偏移，默认 0
#' @param gene_label_circum_limit 逻辑型/向量/列表，可选。限制周向偏移，默认 TRUE
#' @param show_legend 是否显示图例，默认 TRUE
#' @param show_label 是否显示基因标签（覆盖 gene_label_show），默认 NULL
#' @param label_size 标签字号（覆盖 gene_label_size），默认 NULL
#' @param ... 传递给 \code{geom_polygon()} 的其他参数
#'
#' @return ggplot2 layer 列表
#' @export
geom_gene <- function(mapping = NULL, data = NULL,
                      gene_offset = NULL,
                      gene_width = NULL,
                      gene_color_scheme = NULL,
                      gene_colors = NULL,
                      gene_order = NULL,
                      gene_label_show = NULL,
                      gene_label_size = NULL,
                      gene_label_rotation = NULL,
                      gene_label_radial_offset = NULL,
                      gene_label_circum_offset = NULL,
                      gene_label_circum_limit = NULL,
                      show_legend = TRUE,
                      show_label = NULL,
                      label_size = NULL,
                      ...) {
  set_gene_params(list(
    gene_offset               = gene_offset,
    gene_width                = gene_width,
    gene_color_scheme         = gene_color_scheme,
    gene_colors               = gene_colors,
    gene_order                = gene_order,
    gene_label_show           = gene_label_show,
    gene_label_size           = gene_label_size,
    gene_label_rotation       = gene_label_rotation,
    gene_label_radial_offset  = gene_label_radial_offset,
    gene_label_circum_offset  = gene_label_circum_offset,
    gene_label_circum_limit   = gene_label_circum_limit,
    show_label_override       = show_label,
    label_size_override       = label_size
  ))

  layers <- list()

  # new_scale_fill 确保与 ribbon 的 fill scale 独立
  if (requireNamespace("ggnewscale", quietly = TRUE)) {
    layers[[1]] <- ggnewscale::new_scale_fill()
  }

  # 手动配色按注释映射；默认配色按链方向映射。
  gene_scheme <- gene_color_scheme %||% "strand"
  fill_mapping <- if (identical(gene_scheme, "manual")) {
    aes(x = x, y = y, group = group, fill = anno)
  } else {
    aes(x = x, y = y, group = group, fill = strand)
  }

  # 仅返回 layer，不创建 scale（scale 由 print.ggchord 统一设置）
  layers[[length(layers) + 1]] <- ggplot2::layer(
    data        = data.frame(x = numeric(0), y = numeric(0),
                             group = integer(0),
                             strand = factor(character(0), levels = c("+", "-")),
                             anno = character(0), ord = integer(0)),
    mapping     = fill_mapping,
    stat        = "identity",
    geom        = GeomPolygon,
    position    = "identity",
    show.legend = show_legend,
    inherit.aes = FALSE,
    params      = list(color = "black", key_glyph = draw_key_gene_arrow, ...)
  )

  # 占位 Text layer（供基因标签注入；实际数据在 print 时注入）
  layers[[length(layers) + 1]] <- geom_text(
    data        = data.frame(x = numeric(0), y = numeric(0),
                             text_x = numeric(0), text_y = numeric(0),
                             text = character(0), text_angle = numeric(0),
                             hjust = numeric(0), vjust = numeric(0)),
    mapping     = aes(x = text_x, y = text_y, label = text,
                      angle = text_angle, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    show.legend = FALSE
  )

  layers
}
