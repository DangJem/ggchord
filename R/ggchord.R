# ggchord.R — 构造函数
# 弦图 ggplot2 分层 API 的入口点
# v0.4.0: ggchord() 仅保留数据参数和全局样式参数，
# 布局参数移至各 geom_* 图层，坐标计算延迟到 print 阶段

# 全局变量声明（避免 R CMD check 的 NOTE）
globalVariables(c(
  "x", "y", "group", "pident", "fill", "strand", "anno", "seq_id",
  "text_x", "text_y", "text", "text_angle", "hjust", "vjust",
  "x0", "y0", "x1", "y1", "label", "label_x", "label_y", "size"
))

#' ggchord: ggplot2 分层的多序列比对弦图
#'
#' ggchord 采用 ggplot2 分层语法可视化多序列 BLAST 比对结果。
#' \code{ggchord()} 构造函数负责数据校验和全局设置；
#' 各 \code{geom_*} 图层按需叠加，各自负责布局参数和视觉渲染。
#' 布局在 \code{print()} 时延迟计算。
#'
#' @param seq_data data.frame/tibble，必需。序列基本信息
#' @param ribbon_data data.frame/tibble，可选。BLAST 比对结果
#' @param gene_data data.frame/tibble，可选。基因注释数据
#' @param title 字符型。图表主标题，默认 NULL
#' @param rotation 数值型。全局旋转角度（度），默认 45
#' @param panel_margin 数值型/列表，可选。面板边距，默认 0
#' @param show_legend 逻辑型。是否显示图例，默认 TRUE
#' @param debug 逻辑型。是否输出调试信息，默认 FALSE
#'
#' @return ggchord 对象（继承 ggplot），可通过 + 叠加 geom_* 图层
#' @export
#'
#' @examples
#' library(ggchord)
#' data(seq_data_example)
#' data(ribbon_data_example)
#' data(gene_data_example)
#'
#' p <- ggchord(
#'   seq_data = seq_data_example,
#'   ribbon_data = ribbon_data_example,
#'   gene_data = gene_data_example
#' ) +
#'   geom_seq() +
#'   geom_ribbon() +
#'   geom_gene() +
#'   geom_axis()
#' print(p)
#'
#' @import ggplot2
#' @import ggnewscale
#' @import RColorBrewer
#' @import grDevices
#' @import grid
ggchord <- function(
    seq_data,
    ribbon_data = NULL,
    gene_data = NULL,
    title = NULL,
    rotation = 45,
    panel_margin = 0,
    show_legend = TRUE,
    debug = FALSE
) {
  # 清空上一轮的参数
  clear_chord_env()

  # ====================================================================
  # 1. 校验数据
  # ====================================================================
  required_seq_cols <- c("seq_id", "length")
  if (!all(required_seq_cols %in% colnames(seq_data))) {
    stop("seq_data 必须包含以下列：", paste(required_seq_cols, collapse = ", "))
  }
  if (any(seq_data$length <= 0)) {
    stop("seq_data 中的 'length' 值必须为正数")
  }
  if (anyDuplicated(seq_data$seq_id)) {
    stop("seq_data 中的 'seq_id' 必须唯一")
  }

  if (!is.null(ribbon_data)) {
    required_ribbon_cols <- c("qaccver", "saccver", "length", "pident",
                              "qstart", "qend", "sstart", "send")
    if (!all(required_ribbon_cols %in% colnames(ribbon_data))) {
      stop("ribbon_data 必须包含以下列：",
           paste(required_ribbon_cols, collapse = ", "))
    }
    if (nrow(ribbon_data) == 0) warning("ribbon_data 中没有有效比对数据")
    if (debug) cat("比对数据行数：", nrow(ribbon_data), "\n")
  }

  if (!is.null(gene_data)) {
    required_gene_cols <- c("seq_id", "start", "end", "strand", "anno")
    if (!all(required_gene_cols %in% colnames(gene_data))) {
      stop("gene_data 必须包含以下列：",
           paste(required_gene_cols, collapse = ", "))
    }
    if (nrow(gene_data) == 0) warning("gene_data 中没有有效基因注释数据")
    if (any(!gene_data$strand %in% c("+", "-"))) {
      stop("gene_data 中的 'strand' 只能是 '+' 或 '-'")
    }
    if (debug) cat("基因注释数据行数：", nrow(gene_data), "\n")
  }

  # ====================================================================
  # 2. 存储原始数据 + 全局参数
  # ====================================================================
  set_chord_data(seq_data, ribbon_data, gene_data)

  set_global_params(list(
    rotation     = rotation,
    panel_margin = panel_margin,
    show_legend  = show_legend,
    debug        = debug
  ))

  # ====================================================================
  # 3. 构建基础 ggplot 对象（coord 在 print 时被替换为正确范围）
  # ====================================================================
  margin_vals <- process_panel_margin(panel_margin)

  p <- ggplot() +
    coord_chord() +
    labs(title = title) +
    theme(
      plot.title       = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.margin      = margin(t = margin_vals$t, r = margin_vals$r,
                                b = margin_vals$b, l = margin_vals$l),
      legend.background = element_blank(),
      legend.box.spacing  = unit(10, "mm"),
      legend.spacing      = unit(5, "mm"),
      legend.text         = element_text(size = 8),
      legend.title        = element_text(size = 10, face = "bold"),
      axis.title          = element_blank(),
      axis.line           = element_blank(),
      axis.ticks          = element_blank(),
      axis.text           = element_blank(),
      panel.background    = element_blank(),
      legend.position     = if (isTRUE(show_legend)) "right" else "none"
    )

  class(p) <- c("ggchord", class(p))
  p
}

# ====================================================================
# +.ggchord 方法
# ====================================================================

#' @export
`+.ggchord` <- function(e1, e2) {
  if (is.list(e2) && !inherits(e2, "ggplot") && !inherits(e2, "LayerInstance")) {
    p <- e1
    for (elem in e2) {
      p <- p + elem
    }
  } else {
    p <- NextMethod()
  }
  class(p) <- unique(c("ggchord", class(p)))
  p
}

# ====================================================================
# print.ggchord：延迟计算布局 + 注入数据到 layer + 渲染
# ====================================================================

#' @export
print.ggchord <- function(x, ...) {
  # 第一步：收集所有参数
  global <- get_global_params()
  data_list <- get_chord_data()

  seq_params    <- get_seq_params()
  ribbon_params <- get_ribbon_params()
  gene_params   <- get_gene_params()
  axis_params   <- get_axis_params()

  if (is.null(seq_params)) seq_params <- list()

  # --- 处理序列 ---
  seq_data <- data_list$seq_data
  seqs     <- seq_data$seq_id
  lens     <- setNames(seq_data$length, seqs)

  if (!is.null(seq_params$seq_order)) {
    if (!all(seq_params$seq_order %in% seqs)) {
      stop("seq_order 包含未知序列 ID")
    }
    seqs <- seq_params$seq_order
    lens <- lens[seqs]
  }
  n <- length(seqs)

  seq_labels    <- process_sequence_param(seq_params$seq_labels, seqs,
                                          "seq_labels", default_value = seqs)
  seqRadius     <- process_sequence_param(seq_params$seq_radius, seqs,
                                          "seq_radius", 1.0)
  orientation   <- process_sequence_param(seq_params$seq_orientation, seqs,
                                          "seq_orientation", 1)
  seq_gap       <- process_sequence_param(seq_params$seq_gap, seqs,
                                          "seq_gap", 0.03)
  seq_curvature <- process_sequence_param(seq_params$seq_curvature, seqs,
                                          "seq_curvature", 1.0)

  if (any(seqRadius <= 0)) stop("seq_radius 必须为正数")
  if (any(!orientation %in% c(-1, 1))) {
    stop("seq_orientation 只能取 1 或 -1")
  }
  if (any(seq_gap < 0 | seq_gap >= 0.5)) {
    stop("seq_gap 必须位于 [0, 0.5) 区间")
  }

  if (!is.null(seq_params$seq_colors)) {
    seq_colors <- process_sequence_param(seq_params$seq_colors, seqs, "seq_colors")
  } else {
    pal <- if (n <= 9) brewer.pal(n, "Set1")
           else colorRampPalette(brewer.pal(9, "Set1"))(n)
    seq_colors <- setNames(pal, seqs)
  }

  # --- 处理 Ribbon ---
  ribbonGap  <- process_sequence_param(ribbon_params$ribbon_gap %||% 0.15,
                                       seqs, "ribbon_gap", 0.15)
  ribbon_color_scheme <- ribbon_params$ribbon_color_scheme %||% "pident"
  ribbon_alpha    <- ribbon_params$ribbon_alpha %||% 0.35
  ribbon_ctrl_pt  <- ribbon_params$ribbon_ctrl_point %||% c(0, 0)

  ribbon_colors <- ribbon_params$ribbon_colors
  if (!ribbon_color_scheme %in% c("pident", "query", "single")) {
    stop("ribbon_color_scheme 必须是 'pident'、'query' 或 'single'")
  }
  if (!is.numeric(ribbon_alpha) || length(ribbon_alpha) != 1 ||
      ribbon_alpha < 0 || ribbon_alpha > 1) {
    stop("ribbon_alpha 必须位于 [0, 1] 区间")
  }

  # ribbon_colors 校验只在实际有 ribbon_data 时执行
  has_ribbon_data <- !is.null(data_list$ribbon_data) &&
                     nrow(data_list$ribbon_data) > 0

  if (has_ribbon_data) {
    if (is.null(ribbon_colors)) {
      ribbon_colors <- switch(ribbon_color_scheme,
        single = "steelblue",
        query  = {
          mix <- 0.5
          sapply(seq_colors, function(col) {
            cols <- col2rgb(col)
            light_cols <- cols + (255 - cols) * mix
            rgb(light_cols[1,], light_cols[2,], light_cols[3,],
                maxColorValue = 255)
          })
        },
        pident = c("#440154FF","#482878FF","#3E4A89FF","#31688EFF",
                    "#26828EFF","#1F9E89FF","#35B779FF","#6DCD59FF",
                    "#B4DE2CFF","#FDE725FF"))
    }
    if (ribbon_color_scheme == "query") {
      ribbon_colors <- process_sequence_param(ribbon_colors, seqs,
                                              "ribbon_colors")
    } else if (ribbon_color_scheme == "pident" && length(ribbon_colors) < 2) {
      stop("'pident' 配色至少需要两个 ribbon_colors")
    } else if (ribbon_color_scheme == "single" && length(ribbon_colors) < 1) {
      stop("'single' 配色需要至少一个 ribbon_colors")
    }
  }

  # --- 处理基因 ---
  gene_off  <- gene_params$gene_offset %||% 0.1
  gene_w    <- gene_params$gene_width %||% 0.05
  gene_cs   <- gene_params$gene_color_scheme %||% "strand"
  gene_cols <- gene_params$gene_colors
  gene_ord  <- gene_params$gene_order
  gene_ls   <- gene_params$show_label_override %||%
    gene_params$gene_label_show %||% FALSE
  gene_lsz  <- gene_params$label_size_override %||%
    gene_params$gene_label_size %||% 2.5
  gene_lr   <- gene_params$gene_label_rotation %||% 0
  gene_lro  <- gene_params$gene_label_radial_offset %||% 0
  gene_lco  <- gene_params$gene_label_circum_offset %||% 0
  gene_lcl  <- if (is.null(gene_params$gene_label_circum_limit)) TRUE
               else gene_params$gene_label_circum_limit

  if (!gene_cs %in% c("strand", "manual")) {
    stop("gene_color_scheme 必须是 'strand' 或 'manual'")
  }

  geneGap    <- process_gene_param(gene_off, seqs, "gene_offset", 0.1, FALSE)
  geneWidth  <- process_gene_param(gene_w, seqs, "gene_width", 0.05, FALSE)
  geneLabelRadialOffset <- process_gene_param(gene_lro, seqs,
                                              "gene_label_radial_offset", 0, FALSE)
  geneLabelCircumOffset <- process_gene_param(gene_lco, seqs,
                                              "gene_label_circum_offset", 0, FALSE)
  geneLabelCircumLimit  <- process_gene_param(gene_lcl, seqs,
                                              "gene_label_circum_limit", TRUE, TRUE)
  geneLabelRotation     <- process_gene_param(gene_lr, seqs,
                                              "gene_label_rotation", 0, FALSE)

  # --- 处理坐标轴 ---
  show_axis  <- axis_params$show_axis %||% TRUE
  axisGap    <- process_sequence_param(axis_params$axis_gap %||% 0.04,
                                       seqs, "axis_gap", 0.04)
  axisMaj    <- process_sequence_param(axis_params$axis_tick_major_number %||% 5,
                                       seqs, "axis_tick_major_number", 5)
  axisMajLen <- process_sequence_param(axis_params$axis_tick_major_length %||% 0.02,
                                       seqs, "axis_tick_major_length", 0.02)
  axisMin    <- process_sequence_param(axis_params$axis_tick_minor_number %||% 4,
                                       seqs, "axis_tick_minor_number", 4)
  axisMinLen <- process_sequence_param(axis_params$axis_tick_minor_length %||% 0.01,
                                       seqs, "axis_tick_minor_length", 0.01)
  labelSize  <- process_sequence_param(axis_params$axis_label_size %||% 3,
                                       seqs, "axis_label_size", 3)
  labelOffset <- process_sequence_param(axis_params$axis_label_offset %||% 1.5,
                                        seqs, "axis_label_offset", 1.5)
  axisLabelOri <- process_axis_orientation(
    axis_params$axis_label_orientation %||% "horizontal", seqs
  )

  # ====================================================================
  # 第二步：计算布局
  # ====================================================================
  layout <- compute_chord_layout(
    seqs = seqs, lens = lens, seq_labels = seq_labels,
    seq_colors = seq_colors, seqRadius = seqRadius,
    seq_curvature = seq_curvature, orientation = orientation,
    seq_gap = seq_gap,
    ribbon_data = data_list$ribbon_data, ribbonGap = ribbonGap,
    ribbon_color_scheme = ribbon_color_scheme,
    ribbon_colors = ribbon_colors, ribbon_alpha = ribbon_alpha,
    ribbon_ctrl_point = ribbon_ctrl_pt,
    gene_data = data_list$gene_data,
    geneGap = geneGap, geneWidth = geneWidth,
    geneLabelRadialOffset = geneLabelRadialOffset,
    geneLabelCircumOffset = geneLabelCircumOffset,
    geneLabelCircumLimit = geneLabelCircumLimit,
    geneLabelRotation = geneLabelRotation,
    gene_label_show = gene_ls, gene_label_size = gene_lsz,
    gene_color_scheme = gene_cs, gene_colors = gene_cols,
    gene_order = gene_ord,
    axisGap = axisGap, axisMaj = axisMaj, axisMajLen = axisMajLen,
    axisMin = axisMin, axisMinLen = axisMinLen,
    labelSize = labelSize, labelOffset = labelOffset,
    axisLabelOrientation = axisLabelOri,
    show_axis = show_axis,
    rotation = global$rotation, debug = global$debug
  )

  set_chord_layout(layout)

  # ====================================================================
  # 第三步：注入数据到 layers
  # ====================================================================
  # 按 geom 类型区分 layer 归属
  seq_indices    <- integer(0)
  ribbon_indices <- integer(0)
  gene_poly_indices <- integer(0)
  gene_text_indices <- integer(0)
  axis_line_indices <- integer(0)
  axis_seg_indices  <- integer(0)
  axis_text_indices <- integer(0)

  seq_path_assigned <- FALSE
  for (i in seq_along(x$layers)) {
    lyr <- x$layers[[i]]
    gname <- class(lyr$geom)[1]
    data_names <- names(lyr$data)

    if (gname == "GeomPath") {
      # Sequence and axis paths share a geom; their placeholder columns do not.
      if ("seq_id" %in% data_names) {
        # geom_seq(), when present, is the first such path.  A plot containing
        # only geom_axis() must retain its axis path as an axis path.
        if (!seq_path_assigned && !is.null(seq_params)) {
          seq_indices <- c(seq_indices, i)
          seq_path_assigned <- TRUE
        } else {
          axis_line_indices <- c(axis_line_indices, i)
        }
      }
    } else if (gname %in% c("GeomPolygon", "NewGeomPolygon")) {
      # ggnewscale renames the ribbon geom to NewGeomPolygon.  The gene
      # placeholder is distinguished reliably by its strand column.
      if ("strand" %in% data_names) {
        gene_poly_indices <- c(gene_poly_indices, i)
      } else {
        ribbon_indices <- c(ribbon_indices, i)
      }
    } else if (gname == "GeomSegment") {
      axis_seg_indices <- c(axis_seg_indices, i)
    } else if (gname == "GeomText") {
      # Gene labels use text_x/text_y; axis labels use label_x/label_y.
      if ("text_x" %in% data_names) {
        gene_text_indices <- c(gene_text_indices, i)
      } else {
        axis_text_indices <- c(axis_text_indices, i)
      }
    }
  }

  # 注入序列弧线数据
  if (length(seq_indices) > 0 && length(layout$seq_arcs) > 0) {
    arc_df <- do.call(rbind, layout$seq_arcs)
    for (idx in seq_indices) x$layers[[idx]]$data <- arc_df
  }

  # 注入 ribbon 数据
  if (length(ribbon_indices) > 0 && !is.null(layout$ribbon_polys)) {
    for (idx in ribbon_indices) x$layers[[idx]]$data <- layout$ribbon_polys
  }

  # 注入基因数据
  if (length(gene_poly_indices) > 0 && nrow(layout$gene_polys) > 0) {
    for (idx in gene_poly_indices) x$layers[[idx]]$data <- layout$gene_polys
  }

  # 注入基因标签
  if (length(gene_text_indices) > 0 && nrow(layout$gene_labels) > 0) {
    for (idx in gene_text_indices) x$layers[[idx]]$data <- layout$gene_labels
  }

  # 注入轴线数据
  if (length(axis_line_indices) > 0 && nrow(layout$axis_lines) > 0) {
    for (idx in axis_line_indices) x$layers[[idx]]$data <- layout$axis_lines
  }

  # 注入刻度线数据
  if (length(axis_seg_indices) > 0 && nrow(layout$axis_ticks) > 0) {
    for (idx in axis_seg_indices) x$layers[[idx]]$data <- layout$axis_ticks
  }

  # 注入刻度标签数据
  label_data <- subset(layout$axis_ticks, !is.na(label))
  if (length(axis_text_indices) > 0 && nrow(label_data) > 0) {
    for (idx in axis_text_indices) x$layers[[idx]]$data <- label_data
  }

  # ====================================================================
  # 第四步：动态设置 scales
  # 注意：需要临时移除 ggchord class，否则 +.ggchord 和 +.gg 方法冲突
  # ====================================================================
  cls <- class(x)
  class(x) <- setdiff(cls, "ggchord")

  # 序列颜色 scale
  x <- x + scale_color_manual(
    name   = "Seq ID",
    values = layout$seq_colors,
    labels = layout$seq_labels,
    breaks = layout$seqs
  )

  # Ribbon fill scale
  if (!is.null(layout$ribbon_polys)) {
    if (layout$ribbon_color_scheme == "pident") {
      x <- x + scale_fill_stepsn(
        name    = "Identity(%)",
        colours = layout$ribbon_colors,
        limits  = c(0, 100),
        breaks  = c(0, 50, 80, 90, 95, 100),
        guide   = guide_colorbar(
          theme = theme(legend.title.position = "top",
                        legend.key.height = unit(8, "mm")),
          order = 1, position = "left"
        )
      )
    } else {
      x <- x + scale_fill_identity()
    }
  }

  # 基因 fill scale
  if (nrow(layout$gene_polys) > 0) {
    # 如果 ggnewscale 不可用且 ribbon 也用了 fill scale，跳过基因 fill scale
    # 避免 "Continuous value supplied to discrete scale" 冲突
    ggnewscale_ok <- requireNamespace("ggnewscale", quietly = TRUE)
    ribbon_uses_fill <- !is.null(layout$ribbon_polys)
    if (!ggnewscale_ok && ribbon_uses_fill) {
      # 无 ggnewscale 时无法隔离两个 fill scale，跳过基因着色
    } else {
      if (layout$gene_color_scheme == "strand") {
        x <- x + scale_fill_manual(
          name   = "Strand",
          breaks = c("+", "-"),
          values = layout$gene_pal
        )
      } else {
        x <- x + scale_fill_manual(
          name   = "Gene Annotation",
          breaks = layout$final_gene_order,
          values = layout$gene_pal
        )
      }
    }
  }

  # 坐标轴文本 size scale
  x <- x + scale_size_identity()

  # 恢复 ggchord class
  class(x) <- unique(c("ggchord", class(x)))

  # ====================================================================
  # 第五步：更新 coord 范围
  # ====================================================================
  ext <- layout$extremes
  pad <- 0.05 * max(ext$x_max - ext$x_min, ext$y_max - ext$y_min, 1)
  x$coordinates <- coord_fixed(
    ratio = 1,
    xlim  = c(ext$x_min - pad, ext$x_max + pad),
    ylim  = c(ext$y_min - pad, ext$y_max + pad),
    clip  = "off"
  )

  # ====================================================================
  # 第六步：渲染
  # ====================================================================
  cls <- class(x)
  class(x) <- setdiff(cls, c("ggchord"))
  print(x)
  invisible()
}

# ====================================================================
# ggplot_build.ggchord
# ====================================================================

#' @export
ggplot_build.ggchord <- function(plot) {
  try_layout <- tryCatch(get_chord_layout(), error = function(e) NULL)
  if (is.null(try_layout)) {
    warning("布局尚未计算，请通过 print() 渲染。")
  }
  NextMethod()
}
