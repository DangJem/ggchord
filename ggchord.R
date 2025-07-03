#' ggchord: 多序列弦图（带基因注释箭头）
#'
#' 该函数用于绘制包含多个序列的弦图，可展示序列间的比对关系和基因注释信息。
#' 支持自定义序列排列、颜色、连接带样式等多种参数，适合用于基因组比对可视化。
#'
#' @param seq_data data.frame/tibble，包含序列信息，必须包含列：
#'   - seq_id: 序列唯一标识
#'   - length: 序列长度
#' @param ribbon_data data.frame/tibble，包含BLAST比对结果，可选，若提供则必须包含列：
#'   - qaccver: 查询序列ID
#'   - saccver: 目标序列ID
#'   - length: 比对长度
#'   - pident: 序列相似度百分比
#'   - qstart: 查询序列起始位置
#'   - qend: 查询序列结束位置
#'   - sstart: 目标序列起始位置
#'   - send: 目标序列结束位置
#' @param gene_track data.frame/tibble，包含基因注释信息，可选，若提供则必须包含列：
#'   - seq_id: 序列唯一标识
#'   - start: 基因起始位置
#'   - end: 基因结束位置
#'   - strand: 链方向（+或-）
#'   - anno: 基因注释
#' @param title 字符串，图形主标题，默认 "Multi-sequence Chord Diagram with Gene Annotations"
#' @param seq_order 字符向量，可选，指定序列绘制顺序；若为NULL则使用 seq_data 中的顺序
#' @param seq_labels 字符向量或命名向量，可选，序列标签；若为NULL则使用 seq_id
#' @param seq_orientation 数值向量或单值，可选，每条序列的方向：1（正向）或 -1（反向）；默认正向
#' @param seq_gap 数值或向量，长度与序列数一致，定义每条序列头部到下一条序列尾部的弧度比例 [0,0.5)，默认0.05
#' @param seq_radius 数值或向量，序列圆弧半径，支持单值或与序列数相同的向量，默认1.0
#' @param seq_colors 颜色向量或命名向量，定义各序列圆弧颜色；若为NULL则基于 RColorBrewer Set1 自动生成
#' @param gene_offset 数值、向量或列表，基因箭头与序列圆弧之间的径向偏移距离。支持：
#'   - 单值：所有序列的所有链使用相同偏移
#'   - 向量：长度与序列数一致，每个序列的所有链使用相同偏移
#'   - 列表：命名列表，每个元素对应一个序列，元素可为单值（该序列所有链）或包含"+"和"-"的命名向量（区分链），默认0.03
#' @param gene_width 数值或向量，基因箭头宽度，默认0.1
#' @param gene_label_show 逻辑值，是否显示基因标签，默认FALSE
#' @param gene_label_size 数值，基因注释文字大小，默认2.5
#' @param gene_color_scheme 字符，指定基因颜色方案，可选"strand"（按链方向）或"manual"（手动指定），默认"strand"
#' @param gene_colors 颜色向量，用于指定基因箭头的填充色，具体行为取决于gene_color_scheme：
#'   - "strand"模式：支持命名向量（仅"+"/"-"）、非命名向量（先"+"后"-"）或单值（正负链同色），缺省则"+"为红色、"-"为蓝色
#'   - "manual"模式：支持命名向量（对应anno）、非命名向量（截断多余，补齐不足），缺省则使用原始anno颜色
#' @param gene_order 字符向量，可选，指定基因在图例中的显示顺序；若为NULL则使用基因在数据中出现的顺序
#' @param ribbon_color_scheme 字符，连接带配色方案，可选 "single"、"query" 或 "pident"，默认"single"
#' @param ribbon_colors 连接带颜色参数：
#'   - single: 单一颜色（单值或向量取第一个）
#'   - query: 按查询序列映射颜色（命名或非命名向量或单值）
#'   - pident: 渐变色阶向量，用于按相似度百分比生成渐变
#' @param ribbon_alpha 数值，连接带透明度 [0,1]，默认0.4
#' @param ribbon_ctrl_point 数值向量或列表，可选，贝塞尔控制点，长度2或list(c1,c2,..)，默认NULL（圆心）
#' @param ribbon_gap 数值或向量，序列圆弧与连接带之间的径向距离，默认0.1
#' @param axis_gap 数值或向量，坐标轴与序列圆弧之间径向距离，支持负值，默认0.05
#' @param axis_tick_major_number 整数或向量，每条序列主刻度数，默认5
#' @param axis_tick_major_length 数值或向量，主刻度线长度比例，默认0.02
#' @param axis_tick_minor_number 整数或向量，每两个主刻度之间的次刻度数，默认4
#' @param axis_tick_minor_length 数值或向量，次刻度线长度比例，默认0.01
#' @param axis_label_size 数值或向量，坐标轴刻度文字大小，默认3
#' @param axis_label_offset 数值或向量，基于默认标签位置（1.5倍刻度长度）的偏移量，0时为原位置，正值向外，负值向内，默认0
#' @param rotation 数值，整体图形绕原点旋转角度（度，逆时针为正），默认0
#' @param show_legend 逻辑值，是否显示图例，默认TRUE
#' @param show_axis 逻辑值，是否显示坐标轴，默认TRUE
#' @param debug 逻辑值，是否打印调试信息，默认FALSE
#' @return 返回一个 ggplot2 对象
#' @import ggplot2
#' @import RColorBrewer
#' @import grDevices
#' @import ggnewscale
#' @export

# 加载所需包
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(ggnewscale)

# 辅助函数：生成坐标轴主刻度断点
breakPointsFunc <- function(max_value, n = 5, tol = 0.5) {
  if (max_value <= 0) return(c(0, max_value))
  
  ticks <- pretty(c(0, max_value), n = n)
  ticks <- ticks[ticks >= 0 & ticks <= max_value]
  ticks <- sort(unique(c(0, ticks, max_value)))
  
  if (length(ticks) >= 3) {
    d <- diff(ticks)
    med <- median(d[-length(d)])
    if (d[length(d)] < med * tol) {
      ticks <- ticks[-(length(ticks) - 1)]
    }
  }
  
  return(ticks)
}

# 辅助函数：线性映射到角度索引
map_idx <- function(r, angles) {
  r <- pmin(1, pmax(0, ifelse(is.na(r), 0, r)))
  tgt <- angles[1] + r * (angles[length(angles)] - angles[1])
  which.min(abs(angles - tgt))
}

# 辅助函数：生成贝塞尔曲线点
bezier_pts <- function(p0, p3, c1, c2, n = 100) {
  t <- seq(0, 1, length.out = n)
  bx <- (1 - t)^3*p0[1] + 3*(1 - t)^2*t*c1[1] + 3*(1 - t)*t^2*c2[1] + t^3*p3[1]
  by <- (1 - t)^3*p0[2] + 3*(1 - t)^2*t*c1[2] + 3*(1 - t)*t^2*c2[2] + t^3*p3[2]
  data.frame(x=bx, y=by)
}

# 通用参数处理（序列参数）
process_sequence_param <- function(param, seqs, param_name, default_value = NULL, allow_null = FALSE) {
  n <- length(seqs)
  if (is.null(param)) {
    if (allow_null) return(NULL)
    if (!is.null(default_value)) return(setNames(rep(default_value, n), seqs))
    stop(paste0(param_name, " 不能为空且未指定默认值"))
  }
  if (length(param) == 1) return(setNames(rep(param, n), seqs))
  if (length(param) == n && is.null(names(param))) return(setNames(param, seqs))
  if (!is.null(names(param))) {
    if (!all(names(param) %in% seqs)) stop(paste0(param_name, " 包含未知序列ID: ",
                                                  paste(setdiff(names(param), seqs), collapse=", ")))
    return(param[seqs])
  }
  stop(paste0("无法处理的", param_name, "格式，请提供单值、命名向量或与序列数相同的非命名向量"))
}

# 处理gene_offset参数的辅助函数
process_gene_offset <- function(gene_offset, seqs, default = 0.03) {
  n <- length(seqs)
  result <- setNames(lapply(seqs, function(id) c("+" = default, "-" = default)), seqs)
  
  if (is.null(gene_offset)) {
    return(result)
  }
  
  if (length(gene_offset) == 1 && !is.list(gene_offset) && is.numeric(gene_offset)) {
    val <- gene_offset
    return(setNames(lapply(seqs, function(id) c("+" = val, "-" = val)), seqs))
  }
  
  if (is.vector(gene_offset) && !is.list(gene_offset) && length(gene_offset) == n) {
    if (!is.null(names(gene_offset))) {
      for (id in seqs) {
        if (id %in% names(gene_offset)) {
          val <- gene_offset[id]
          result[[id]] <- c("+" = val, "-" = val)
        }
      }
    } else {
      for (i in seq_along(seqs)) {
        val <- gene_offset[i]
        result[[seqs[i]]] <- c("+" = val, "-" = val)
      }
    }
    return(result)
  }
  
  if (is.list(gene_offset)) {
    list_names <- names(gene_offset)
    if (is.null(list_names)) {
      if (length(gene_offset) != n) {
        stop("gene_offset 列表长度与序列数不匹配，且未命名")
      }
      for (i in seq_along(seqs)) {
        elem <- gene_offset[[i]]
        if (length(elem) == 1 && is.numeric(elem)) {
          result[[seqs[i]]] <- c("+" = elem, "-" = elem)
        } else if (is.vector(elem) && all(names(elem) %in% c("+", "-"))) {
          current <- result[[seqs[i]]]
          if ("+" %in% names(elem)) current["+"] <- elem["+"]
          if ("-" %in% names(elem)) current["-"] <- elem["-"]
          result[[seqs[i]]] <- current
        } else {
          stop(paste("gene_offset 列表中第", i, "个元素格式错误，应为单值或包含'+', '-'的命名向量"))
        }
      }
    } else {
      for (id in list_names) {
        if (!id %in% seqs) {
          stop(paste("gene_offset 列表包含未知序列ID:", id))
        }
        elem <- gene_offset[[id]]
        current <- result[[id]]
        if (length(elem) == 1 && is.numeric(elem)) {
          current <- c("+" = elem, "-" = elem)
        } else if (is.vector(elem) && all(names(elem) %in% c("+", "-"))) {
          if ("+" %in% names(elem)) current["+"] <- elem["+"]
          if ("-" %in% names(elem)) current["-"] <- elem["-"]
        } else {
          stop(paste("gene_offset 列表中元素", id, "格式错误，应为单值或包含'+', '-'的命名向量"))
        }
        result[[id]] <- current
      }
    }
    return(result)
  }
  
  stop("gene_offset 格式错误，请提供单值、向量或列表")
}

# 处理gene_colors参数的辅助函数（strand模式）
process_strand_colors <- function(gene_colors) {
  default <- c("+" = "#E41A1C", "-" = "#377EB8")
  if (is.null(gene_colors)) {
    return(default)
  }
  if (!is.null(names(gene_colors))) {
    if (!all(names(gene_colors) %in% c("+", "-"))) {
      stop("'strand'模式下，gene_colors命名向量只能包含'+'和'-'")
    }
    res <- default
    res[names(gene_colors)] <- gene_colors
    return(res)
  } else {
    len <- length(gene_colors)
    if (len == 1) {
      return(c("+" = gene_colors[1], "-" = gene_colors[1]))
    } else if (len == 2) {
      return(c("+" = gene_colors[1], "-" = gene_colors[2]))
    } else {
      stop("'strand'模式下，非命名gene_colors长度必须为1或2")
    }
  }
}

# 处理gene_colors参数的辅助函数（manual模式）
process_manual_colors <- function(gene_colors, unique_anno, gene_order) {
  if (!is.null(gene_order)) {
    unknown <- setdiff(gene_order, unique_anno)
    if (length(unknown) > 0) {
      stop("'gene_order' 包含未知基因注释: ", paste(unknown, collapse = ","))
    }
    final_order <- c(gene_order, setdiff(unique_anno, gene_order))
  } else {
    final_order <- unique_anno
  }
  n_anno <- length(final_order)
  
  if (n_anno <= 9) {
    default_pal <- brewer.pal(n_anno, "Set1")
  } else {
    default_pal <- colorRampPalette(brewer.pal(9, "Set1"))(n_anno)
  }
  names(default_pal) <- final_order
  
  if (is.null(gene_colors)) {
    return(default_pal)
  }
  
  if (!is.null(names(gene_colors))) {
    unknown <- setdiff(names(gene_colors), unique_anno)
    if (length(unknown) > 0) {
      warning("'manual'模式下，gene_colors包含未知注释: ", paste(unknown, collapse = ","))
    }
    res <- default_pal
    common_names <- intersect(names(gene_colors), final_order)
    res[common_names] <- gene_colors[common_names]
    return(res)
  } else {
    len <- length(gene_colors)
    res <- character(n_anno)
    if (len >= 1) {
      res[1:min(len, n_anno)] <- gene_colors[1:min(len, n_anno)]
    }
    if (len < n_anno) {
      res[(len + 1):n_anno] <- default_pal[(len + 1):n_anno]
    }
    names(res) <- final_order
    return(res)
  }
}

# 主函数：绘制多序列弦图
ggchord <- function(
    seq_data,
    ribbon_data            = NULL,
    gene_track             = NULL,
    title                  = "Multi-sequence Chord Diagram with Gene Annotations",
    seq_order              = NULL,
    seq_labels             = NULL,
    seq_orientation        = NULL,
    seq_gap                = 0.05,
    seq_radius             = 1.0,
    seq_colors             = NULL,
    gene_offset            = 0.03,
    gene_width             = 0.1,
    gene_label_show        = FALSE,
    gene_label_size        = 2.5,
    gene_color_scheme      = c("strand", "manual"),
    gene_colors            = NULL,
    gene_order             = NULL,
    ribbon_color_scheme    = c("single","query","pident"),
    ribbon_colors          = NULL,
    ribbon_alpha           = 0.4,
    ribbon_ctrl_point      = NULL,
    ribbon_gap             = 0.1,
    axis_gap               = 0.05,
    axis_tick_major_number = 5,
    axis_tick_major_length = 0.02,
    axis_tick_minor_number = 4,
    axis_tick_minor_length = 0.01,
    axis_label_size        = 3,
    axis_label_offset      = 0.1,
    rotation               = 45,
    show_legend            = TRUE,
    show_axis              = TRUE,
    debug                  = FALSE
) {
  # 调试信息：函数开始
  if (debug) cat("[DEBUG] 开始执行ggchord函数\n")
  
  # 检查必要的包是否安装
  if (!"ggplot2" %in% installed.packages()) stop("需要安装 ggplot2 包")
  if (!"RColorBrewer" %in% installed.packages()) stop("需要安装 RColorBrewer 包")
  if (!"grDevices" %in% installed.packages()) stop("需要安装 grDevices 包")
  if (!"ggnewscale" %in% installed.packages()) stop("需要安装 ggnewscale 包以支持多填充色映射")
  
  ribbon_color_scheme <- match.arg(ribbon_color_scheme)
  gene_color_scheme <- match.arg(gene_color_scheme)
  
  # 1. 验证输入数据格式
  if (debug) cat("[DEBUG] 开始验证输入数据格式\n")
  
  # 验证序列数据（必选）
  required_seq_cols <- c("seq_id", "length")
  if (!all(required_seq_cols %in% colnames(seq_data))) {
    stop("seq_data 必须包含以下列: ", paste(required_seq_cols, collapse = ", "))
  }
  if (any(seq_data$length <= 0)) {
    stop("seq_data 中的 length 必须为正数")
  }
  if (debug) {
    cat("[DEBUG] 序列数据验证通过，包含", nrow(seq_data), "条序列\n")
    print(data.frame(seq_id = seq_data$seq_id, length = seq_data$length))
  }
  
  # 验证BLAST数据（可选）
  fb_all <- NULL
  if (!is.null(ribbon_data)) {
    required_ribbon_cols <- c("qaccver", "saccver", "length", "pident", 
                              "qstart", "qend", "sstart", "send")
    if (!all(required_ribbon_cols %in% colnames(ribbon_data))) {
      stop("ribbon_data 必须包含以下列: ", paste(required_ribbon_cols, collapse = ", "))
    }
    fb_all <- ribbon_data
    if (nrow(fb_all) == 0) warning("ribbon_data 中没有有效比对数据")
    if (debug) {
      cat("[DEBUG] 比对数据验证通过，共", nrow(fb_all), "行\n")
      cat("[DEBUG] 比对数据包含的查询序列ID:", paste(unique(fb_all$qaccver), collapse = ", "), "\n")
      cat("[DEBUG] 比对数据包含的目标序列ID:", paste(unique(fb_all$saccver), collapse = ", "), "\n")
    }
  } else if (debug) {
    cat("[DEBUG] 未提供ribbon_data，不绘制连接带\n")
  }
  
  # 验证基因注释数据（可选）
  if (!is.null(gene_track)) {
    required_gene_cols <- c("seq_id", "start", "end", "strand", "anno")
    if (!all(required_gene_cols %in% colnames(gene_track))) {
      stop("gene_track 必须包含以下列: ", paste(required_gene_cols, collapse = ", "))
    }
    if (nrow(gene_track) == 0) warning("gene_track 中没有有效基因注释数据")
    if (debug) {
      cat("[DEBUG] 基因数据验证通过，共", nrow(gene_track), "个基因\n")
      cat("[DEBUG] 基因涉及的序列ID:", paste(unique(gene_track$seq_id), collapse = ", "), "\n")
      cat("[DEBUG] 基因注释类别:", paste(unique(gene_track$anno), collapse = ", "), "\n")
    }
    
    if (!is.null(gene_order)) {
      unknown_genes <- setdiff(gene_order, unique(gene_track$anno))
      if (length(unknown_genes) > 0) {
        stop("gene_order 包含未知基因注释: ", paste(unknown_genes, collapse = ", "))
      } else if (debug) {
        cat("[DEBUG] gene_order验证通过，包含已知基因注释\n")
      }
    }
  } else if (debug) {
    cat("[DEBUG] 未提供gene_track，不绘制基因注释\n")
  }
  
  # 2. 处理序列信息
  if (debug) cat("[DEBUG] 开始处理序列信息\n")
  
  seqs <- seq_data$seq_id
  lens <- setNames(seq_data$length, seqs)
  if (length(seqs) < 1) stop("seq_data 中至少需要包含1条序列")
  
  # 处理序列顺序
  if (!is.null(seq_order)) {
    if (!all(seq_order %in% seqs)) {
      stop("seq_order 包含未知序列ID: ", paste(setdiff(seq_order, seqs), collapse = ", "))
    }
    seqs <- seq_order
    lens <- lens[seqs]
    if (debug) cat("[DEBUG] 已按seq_order调整序列顺序为:", paste(seqs, collapse = ", "), "\n")
  } else if (debug) {
    cat("[DEBUG] 未指定seq_order，使用seq_data中的顺序:", paste(seqs, collapse = ", "), "\n")
  }
  n <- length(seqs)
  
  # 3. 处理序列相关参数
  if (debug) cat("[DEBUG] 开始处理序列相关参数\n")
  
  seq_labels <- process_sequence_param(seq_labels, seqs, "seq_labels", default_value = seqs)
  if (debug) {
    cat("[DEBUG] 处理后的seq_labels:\n")
    print(seq_labels)
  }
  
  seqRadius <- process_sequence_param(seq_radius, seqs, "seq_radius", default_value = 1.0)
  if (debug) {
    cat("[DEBUG] 处理后的seq_radius:\n")
    print(seqRadius)
  }
  
  ribbonGap <- process_sequence_param(ribbon_gap, seqs, "ribbon_gap", default_value = 0.1)
  axisGap <- process_sequence_param(axis_gap, seqs, "axis_gap", default_value = 0.05)
  geneGap <- process_gene_offset(gene_offset, seqs, default = 0.03)
  if (debug) {
    cat("[DEBUG] 处理后的gene_offset（按序列和链）:\n")
    print(geneGap)
  }
  
  geneWidth <- process_sequence_param(gene_width, seqs, "gene_width", default_value = 0.1)
  axisMaj <- process_sequence_param(axis_tick_major_number, seqs, "axis_tick_major_number", default_value = 5)
  axisMajLen <- process_sequence_param(axis_tick_major_length, seqs, "axis_tick_major_length", default_value = 0.02)
  axisMin <- process_sequence_param(axis_tick_minor_number, seqs, "axis_tick_minor_number", default_value = 4)
  axisMinLen <- process_sequence_param(axis_tick_minor_length, seqs, "axis_tick_minor_length", default_value = 0.01)
  labelSize <- process_sequence_param(axis_label_size, seqs, "axis_label_size", default_value = 3)
  labelOffset <- process_sequence_param(axis_label_offset, seqs, "axis_label_offset", default_value = 0)
  
  orientation <- process_sequence_param(seq_orientation, seqs, "seq_orientation", default_value = 1)
  if (debug) {
    cat("[DEBUG] 处理后的seq_orientation（1=正向，-1=反向）:\n")
    print(orientation)
  }
  
  rot_rad <- rotation * pi / 180
  if (debug) cat("[DEBUG] 图形整体旋转角度:", rotation, "度（", rot_rad, "弧度）\n")
  
  # 处理序列间隙参数
  seq_gap <- process_sequence_param(seq_gap, seqs, "seq_gap")
  if (any(seq_gap < 0 | seq_gap >= 0.5)) {
    stop("seq_gap 必须在 [0, 0.5) 范围内")
  }
  if (debug) {
    cat("[DEBUG] 处理后的seq_gap（弧度比例）:\n")
    print(seq_gap)
  }
  
  # 处理序列颜色
  if (is.null(seq_colors)) {
    pal <- if (length(seqs) <= 9) {
      brewer.pal(length(seqs), "Set1")
    } else {
      colorRampPalette(brewer.pal(9, "Set1"))(length(seqs))
    }
    seq_colors <- setNames(pal, seqs)
  } else {
    seq_colors <- process_sequence_param(seq_colors, seqs, "seq_colors")
  }
  if (debug) {
    cat("[DEBUG] 序列颜色分配:\n")
    print(seq_colors)
  }
  
  # 4. 处理连接带颜色参数
  if (!is.null(ribbon_data)) {
    if (debug) cat("[DEBUG] 开始处理连接带颜色参数（scheme:", ribbon_color_scheme, "）\n")
    
    if (is.null(ribbon_colors)) {
      ribbon_colors <- switch(ribbon_color_scheme,
                              single = "steelblue",
                              query = {
                                mix <- 0.5  
                                sapply(seq_colors, function(col) {
                                  cols <- col2rgb(col)
                                  light_cols <- cols + (255 - cols) * mix
                                  rgb(light_cols[1,], light_cols[2,], light_cols[3,], maxColorValue = 255)
                                })
                              },
                              pident = c(
                                "#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF",
                                "#1F9E89FF","#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF"
                              )
      )
    }
    
    if (ribbon_color_scheme == "single") {
      singleCol <- if (length(ribbon_colors) > 1) ribbon_colors[[1]] else ribbon_colors
      if (debug) cat("[DEBUG] 连接带单一颜色:", singleCol, "\n")
    } else if (ribbon_color_scheme == "query") {
      queryCols <- process_sequence_param(ribbon_colors, seqs, "ribbon_colors")
      if (debug) {
        cat("[DEBUG] 连接带按查询序列着色:\n")
        print(queryCols)
      }
    } else if (ribbon_color_scheme == "pident") {
      if (length(ribbon_colors) < 2) {
        stop("pident 模式需指定至少两个颜色用作渐变色阶")
      }
      rampFunc <- colorRampPalette(ribbon_colors)
      if (debug) {
        cat("[DEBUG] 连接带按相似度着色，渐变色阶:", paste(ribbon_colors, collapse = " → "), "\n")
      }
    }
  }
  
  # 5. 计算基因颜色
  gene_pal <- NULL
  final_gene_order <- NULL
  if (!is.null(gene_track) && nrow(gene_track) > 0) {
    if (debug) cat("[DEBUG] 开始处理基因颜色参数（scheme:", gene_color_scheme, "）\n")
    
    valid_genes <- gene_track[gene_track$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      unique_anno <- unique(valid_genes$anno)
      if (!is.null(gene_order)) {
        final_order <- c(gene_order, setdiff(unique_anno, gene_order))
      } else {
        final_order <- unique_anno
      }
      final_gene_order <- final_order
      
      if (gene_color_scheme == "strand") {
        gene_pal <- process_strand_colors(gene_colors)
        if (debug) {
          cat("[DEBUG] 基因按链着色（+/-）:\n")
          print(gene_pal)
        }
      } else if (gene_color_scheme == "manual") {
        gene_pal <- process_manual_colors(gene_colors, unique_anno, gene_order)
        if (debug) {
          cat("[DEBUG] 基因按注释手动着色:\n")
          print(gene_pal)
        }
      }
    } else if (debug) {
      cat("[DEBUG] 基因数据中无有效序列ID匹配，不绘制基因箭头\n")
    }
  }
  
  # 6. 计算序列弧度和间隙弧度
  if (debug) cat("[DEBUG] 开始计算序列弧度和间隙\n")
  
  total_circ <- 2 * pi
  total_gap_prop <- sum(seq_gap)
  
  if (total_gap_prop >= 1) {
    stop("seq_gap 总和不能超过1（无法容纳所有序列）")
  }
  if (debug) {
    cat("[DEBUG] 总间隙比例:", total_gap_prop, "（小于1，有效）\n")
  }
  
  seq_total_prop <- 1 - total_gap_prop
  sum_lens <- sum(lens)
  theta <- (lens / sum_lens) * total_circ * seq_total_prop  # 每个序列的弧度
  gap_rads <- total_circ * seq_gap  # 每个间隙的实际弧度
  
  if (debug) {
    cat("[DEBUG] 各序列弧度（theta）:\n")
    print(theta)
    cat("[DEBUG] 各间隙弧度（gap_rads）:\n")
    print(gap_rads)
  }
  
  # 7. 计算序列起始和结束角度
  starts <- numeric(n)
  starts[1] <- 0
  
  if (n > 1) {
    for (i in 2:n) {
      starts[i] <- starts[i-1] + theta[i-1] + gap_rads[i-1]
    }
  }
  
  ends <- starts + theta
  names(starts) <- names(ends) <- seqs
  
  if (debug) {
    cat("[DEBUG] 序列起始角度（starts，弧度）:\n")
    print(starts)
    cat("[DEBUG] 序列结束角度（ends，弧度）:\n")
    print(ends)
  }
  
  # 8. 准备序列外层和内层坐标（略过详细输出，避免冗余）
  if (debug) cat("[DEBUG] 生成序列弧线坐标（外层和内层）\n")
  
  nSeg <- 500
  seqArcs <- lapply(seqs, function(id) {
    angs <- seq(starts[id], ends[id], length.out = nSeg)
    if (orientation[id] == -1) angs <- rev(angs)
    data.frame(x = seqRadius[id] * cos(angs), 
               y = seqRadius[id] * sin(angs), 
               seq_id = id)
  })
  
  innerArcs <- lapply(seqs, function(id) {
    angs <- seq(starts[id], ends[id], length.out = nSeg)
    if (orientation[id] == -1) angs <- rev(angs)
    r <- seqRadius[id] - ribbonGap[id]
    data.frame(x = r * cos(angs), y = r * sin(angs))
  })
  names(innerArcs) <- seqs
  
  # 9. 生成坐标轴刻度和线（输出刻度断点示例）
  if (debug) cat("[DEBUG] 生成坐标轴刻度数据\n")
  
  axisTicks <- do.call(rbind, lapply(seqs, function(id) {
    majors <- breakPointsFunc(lens[id], axisMaj[id])
    if (debug && id == seqs[1]) {  # 仅输出第一个序列的刻度示例
      cat("[DEBUG] 序列", id, "的主刻度断点:", paste(majors, collapse = ", "), "\n")
    }
    minors <- unlist(lapply(seq_len(length(majors)-1), function(i) {
      seq(majors[i], majors[i+1], length.out = axisMin[id] + 2)[-c(1, axisMin[id] + 2)]
    }))
    pts <- c(majors, minors)
    do.call(rbind, lapply(pts, function(v) {
      frac <- v / lens[id]
      ang <- if (orientation[id] == 1) {
        starts[id] + frac * (ends[id] - starts[id])
      } else {
        ends[id] - frac * (ends[id] - starts[id])
      }
      base_r <- seqRadius[id] + axisGap[id]
      tick_len <- if (v %in% majors) axisMajLen[id] else axisMinLen[id]
      dir <- if (axisGap[id] >= 0) 1 else -1
      end_r <- base_r + tick_len * dir
      lab_r <- base_r + (1.5 + labelOffset[id]) * tick_len * dir
      
      data.frame(
        x0 = base_r * cos(ang), y0 = base_r * sin(ang),
        x1 = end_r * cos(ang), y1 = end_r * sin(ang),
        label = if (v %in% majors) as.character(v) else NA,
        label_x = lab_r * cos(ang), 
        label_y = lab_r * sin(ang),
        size = labelSize[id], 
        seq_id = id
      )
    }))
  }))
  
  axisLines <- do.call(rbind, lapply(seqs, function(id) {
    angs <- seq(starts[id], ends[id], length.out = nSeg)
    if (orientation[id] == -1) angs <- rev(angs)
    r <- seqRadius[id] + axisGap[id]
    data.frame(x = r * cos(angs), y = r * sin(angs), seq_id = id)
  }))
  
  # 10. 生成连接带数据
  allRibbon <- NULL
  if (!is.null(ribbon_data) && !is.null(fb_all) && nrow(fb_all) > 0) {
    if (debug) cat("[DEBUG] 开始生成连接带数据\n")
    
    ribbons <- list()
    cntValid <- cntInvalid <- 0
    invalidReasons <- list(self = 0, noSeq = 0, tooShort = 0)  # 分类统计无效原因
    
    for (i in seq_len(nrow(fb_all))) {
      row <- fb_all[i,]
      q <- row$qaccver
      s <- row$saccver
      
      # 分类统计无效原因
      if (q == s) {
        invalidReasons$self <- invalidReasons$self + 1
        cntInvalid <- cntInvalid + 1
        next
      }
      if (!q %in% seqs || !s %in% seqs) {
        invalidReasons$noSeq <- invalidReasons$noSeq + 1
        cntInvalid <- cntInvalid + 1
        next
      }
      
      # 获取序列的内层弧线数据
      qr <- innerArcs[[q]]
      sr <- innerArcs[[s]]
      
      qi0 <- map_idx((row$qstart - 1)/lens[q], seq(starts[q], ends[q], length.out = nSeg))
      qi1 <- map_idx((row$qend - 1)/lens[q], seq(starts[q], ends[q], length.out = nSeg))
      si0 <- map_idx((row$sstart - 1)/lens[s], seq(starts[s], ends[s], length.out = nSeg))
      si1 <- map_idx((row$send - 1)/lens[s], seq(starts[s], ends[s], length.out = nSeg))
      
      if (abs(qi1 - qi0) < 1 || abs(si1 - si0) < 1) {
        invalidReasons$tooShort <- invalidReasons$tooShort + 1
        cntInvalid <- cntInvalid + 1
        next
      }
      
      segQ <- qr[qi0:qi1, , drop = FALSE]
      segS <- sr[si0:si1, , drop = FALSE]
      
      if (!is.null(ribbon_ctrl_point)) {
        if (is.list(ribbon_ctrl_point)) {
          cp <- ribbon_ctrl_point[[i]]
          c1 <- cp[[1]]
          c2 <- cp[[2]]
        } else {
          c1 <- c2 <- c(0, 0)
        }
      } else {
        c1 <- c2 <- c(0, 0)
      }
      
      b1 <- bezier_pts(as.numeric(segQ[1, ]), as.numeric(segS[1, ]), c1, c1)
      b2 <- bezier_pts(as.numeric(segQ[nrow(segQ), ]), as.numeric(segS[nrow(segS), ]), c2, c2)
      
      poly <- rbind(segQ, b2, segS[nrow(segS):1, ], b1[nrow(b1):1, ])
      if (ribbon_color_scheme == "pident") {
        poly$pident <- row$pident
      } else {
        poly$fill <- switch(ribbon_color_scheme, single = singleCol, query = queryCols[q])
      }
      poly$group <- i
      ribbons[[length(ribbons) + 1]] <- poly
      cntValid <- cntValid + 1
    }
    
    if (debug) {
      cat("[DEBUG] 连接带统计: 有效=", cntValid, ", 无效=", cntInvalid, "\n")
      cat("[DEBUG] 无效原因: 自身比对=", invalidReasons$self, 
          ", 序列不匹配=", invalidReasons$noSeq, 
          ", 比对过短=", invalidReasons$tooShort, "\n")
    }
    if (cntValid > 0) {
      allRibbon <- do.call(rbind, ribbons)
    } else {
      allRibbon <- NULL
      warning("没有有效连接带可绘制")
    }
  }
  
  # 11. 处理基因注释箭头
  gene_polys <- data.frame()
  if (!is.null(gene_track) && nrow(gene_track) > 0) {
    valid_genes <- gene_track[gene_track$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      if (debug) cat("[DEBUG] 开始生成基因箭头数据，共", nrow(valid_genes), "个有效基因\n")
      
      for (i in seq_len(nrow(valid_genes))) {
        gene   <- valid_genes[i, ]
        sid    <- gene$seq_id
        strand <- gene$strand
        anno   <- gene$anno
        
        st_rel <- min(gene$start, gene$end) / lens[sid]
        en_rel <- max(gene$start, gene$end) / lens[sid]
        if (orientation[sid] == 1) {
          a_s <- starts[sid] + st_rel*(ends[sid]-starts[sid])
          a_e <- starts[sid] + en_rel*(ends[sid]-starts[sid])
        } else {
          a_s <- ends[sid] - st_rel*(ends[sid]-starts[sid])
          a_e <- ends[sid] - en_rel*(ends[sid]-starts[sid])
        }
        if (strand == "-") {
          tmp <- a_s; a_s <- a_e; a_e <- tmp
        }
        
        if (debug && i <= 3) {  # 输出前3个基因的角度信息
          cat("[DEBUG] 基因", i, "（", anno, "）: 序列=", sid, ", 链=", strand, 
              ", 起始角度=", a_s, ", 结束角度=", a_e, "\n")
        }
        
        r0        <- seqRadius[sid] + geneGap[[sid]][strand]
        width     <- geneWidth[sid]
        r_out     <- r0 + width/2
        r_in      <- r0 - width/2
        head_frac <- 0.3
        a_break   <- a_s + (1-head_frac)*(a_e - a_s)
        
        n_body <- 40
        ang_body <- seq(a_s, a_break, length.out = n_body)
        body_out <- data.frame(
          x = r_out * cos(ang_body),
          y = r_out * sin(ang_body),
          group = i, ord = seq_len(n_body),
          anno = anno,
          strand = strand
        )
        body_in <- data.frame(
          x = r_in * cos(rev(ang_body)),
          y = r_in * sin(rev(ang_body)),
          group = i, ord = seq_len(n_body) + 3*n_body,
          anno = anno,
          strand = strand
        )
        
        n_head <- 40
        t_seq  <- seq(0, 1, length.out = n_head)
        ang_out <- a_break + t_seq*(a_e - a_break)
        r_out_edge <- r_out + t_seq*(r0 - r_out)
        edge_out <- data.frame(
          x = r_out_edge * cos(ang_out),
          y = r_out_edge * sin(ang_out),
          group = i, ord = seq_len(n_head) + n_body,
          anno = anno,
          strand = strand
        )
        ang_in <- a_break + (1 - t_seq)*(a_e - a_break)
        r_in_edge <- r_in + (1 - t_seq)*(r0 - r_in)
        edge_in <- data.frame(
          x = r_in_edge * cos(ang_in),
          y = r_in_edge * sin(ang_in),
          group = i, ord = seq_len(n_head) + 2*n_body,
          anno = anno,
          strand = strand
        )
        
        gene_polys <- rbind(gene_polys, body_out, edge_out, edge_in, body_in)
      }
    } else if (debug) {
      cat("[DEBUG] 无有效基因可绘制（基因序列ID与主序列不匹配）\n")
    }
  }
  
  # 12. 处理基因标签
  gene_arrows <- data.frame()
  if (!is.null(gene_track) && nrow(gene_track) > 0 && gene_label_show) {
    valid_genes <- gene_track[gene_track$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      gene_arrows <- do.call(rbind, lapply(1:nrow(valid_genes), function(i) {
        gene <- valid_genes[i, ]
        sid <- gene$seq_id
        mid_pos <- (gene$start + gene$end) / 2
        mid_rel <- mid_pos / lens[sid]
        if (orientation[sid] == 1) {
          mid_ang <- starts[sid] + mid_rel * (ends[sid] - starts[sid])
        } else {
          mid_ang <- ends[sid] - mid_rel * (ends[sid] - starts[sid])
        }
        if (debug && i <= 3) {
          cat("[DEBUG] 基因", i, "标签位置角度:", mid_ang, "弧度\n")
        }
        
        r0 <- seqRadius[sid] + geneGap[[sid]][gene$strand]
        label_r <- r0 + ifelse(gene$strand == "+", 0.1, -0.1)
        x_start <- r0 * cos(mid_ang)
        y_start <- r0 * sin(mid_ang)
        x_end <- (r0 + 0.05) * cos(mid_ang)
        y_end <- (r0 + 0.05) * sin(mid_ang)
        text_x <- label_r * cos(mid_ang)
        text_y <- label_r * sin(mid_ang)
        text_angle <- mid_ang * 180 / pi
        
        data.frame(
          x_start = x_start, y_start = y_start, x_end = x_end, y_end = y_end,
          text = gene$anno, text_x = text_x, text_y = text_y, 
          text_angle = text_angle, seq_id = gene$seq_id, group = i
        )
      }))
    }
  }
  
  # 应用旋转
  rotate_df <- function(df) {
    if (all(c("x","y") %in% names(df))) {
      x0 <- df$x; y0 <- df$y
      df$x <-  x0 * cos(rot_rad) - y0 * sin(rot_rad)
      df$y <-  x0 * sin(rot_rad) + y0 * cos(rot_rad)
    }
    if (all(c("x0","y0","x1","y1") %in% names(df))) {
      X <- df$x0; Y <- df$y0
      df$x0 <-  X * cos(rot_rad) - Y * sin(rot_rad)
      df$y0 <-  X * sin(rot_rad) + Y * cos(rot_rad)
      X1 <- df$x1; Y1 <- df$y1
      df$x1 <-  X1 * cos(rot_rad) - Y1 * sin(rot_rad)
      df$y1 <-  X1 * sin(rot_rad) + Y1 * cos(rot_rad)
    }
    if (all(c("label_x","label_y") %in% names(df))) {
      LX <- df$label_x; LY <- df$label_y
      df$label_x <-  LX * cos(rot_rad) - LY * sin(rot_rad)
      df$label_y <-  LX * sin(rot_rad) + LY * cos(rot_rad)
    }
    if (all(c("x_start","y_start","x_end","y_end") %in% names(df))) {
      Xs <- df$x_start; Ys <- df$y_start
      df$x_start <-  Xs * cos(rot_rad) - Ys * sin(rot_rad)
      df$y_start <-  Xs * sin(rot_rad) + Ys * cos(rot_rad)
      Xe <- df$x_end; Ye <- df$y_end
      df$x_end <-  Xe * cos(rot_rad) - Ye * sin(rot_rad)
      df$y_end <-  Xe * sin(rot_rad) + Ye * cos(rot_rad)
    }
    if (all(c("text_x","text_y") %in% names(df))) {
      TX <- df$text_x; TY <- df$text_y
      df$text_x <-  TX * cos(rot_rad) - TY * sin(rot_rad)
      df$text_y <-  TX * sin(rot_rad) + TY * cos(rot_rad)
      df$text_angle <- df$text_angle + rotation
    }
    df
  }
  
  if (debug) cat("[DEBUG] 应用图形旋转（", rotation, "度）\n")
  seqArcs   <- lapply(seqArcs,   rotate_df)
  innerArcs <- lapply(innerArcs, rotate_df)
  axisLines <- rotate_df(axisLines)
  axisTicks <- rotate_df(axisTicks)
  if (!is.null(allRibbon)) {
    allRibbon <- rotate_df(allRibbon)
  }
  if (nrow(gene_arrows) > 0) {
    gene_arrows <- rotate_df(gene_arrows)
  }
  if (nrow(gene_polys) > 0) {
    gene_polys <- rotate_df(gene_polys)
    gene_polys <- gene_polys[with(gene_polys, order(group, ord)), ]
  }
  
  # 13. 绘制图形
  if (debug) cat("[DEBUG] 开始构建ggplot图形对象\n")
  
  p <- ggplot() +
    { if (!is.null(allRibbon))
      geom_polygon(
        data = allRibbon,
        aes(
          x, y, group = group,
          fill = if (ribbon_color_scheme == "pident") pident else fill
        ),
        alpha = ribbon_alpha, 
        color = "grey40", 
        linewidth = 0.2
      )
    } +
    { if (!is.null(allRibbon)) {
      if (ribbon_color_scheme == "pident") {
        scale_fill_stepsn(
          name = "Identity(%)",
          colours = ribbon_colors,
          limits = c(0, 100), 
          breaks = c(0,50,80,90,95,100),
          guide = if (show_legend) guide_colorbar(theme = theme(legend.title.position = "top",legend.key.height = unit(90,"mm")),order = 1,position = "left") else "none"
        )
      } else {
        scale_fill_identity(guide = "none")
      }
    }
    } +
    { if (nrow(gene_polys) > 0) new_scale_fill() } +
    { if (nrow(gene_polys) > 0)
      geom_polygon(
        data  = gene_polys,
        aes(x = x, y = y, group = group, 
            fill = if (gene_color_scheme == "strand") strand else anno),
        color = "black"
      )
    } +
    { if (nrow(gene_polys) > 0)
      scale_fill_manual(
        name = if (gene_color_scheme == "strand") "Strand" else "Gene Annotation",
        breaks = if (gene_color_scheme == "strand") c("+","-") else final_gene_order,
        values = gene_pal,
        guide = if (show_legend) guide_legend(order = 3) else "none"
      )
    } +
    geom_path(
      data = do.call(rbind, seqArcs),
      aes(x, y, group = seq_id, color = seq_id),
      linewidth = 1.2,
      arrow = arrow(angle = 25, length = unit(2, "mm"), type = "closed")
    ) +
    { if (nrow(gene_arrows) > 0 && gene_label_show)
      geom_text(
        data = gene_arrows,
        aes(x = text_x, y = text_y, label = text, angle = text_angle),
        size = gene_label_size,
        color = "black",
        inherit.aes = FALSE
      )
    } +
    { if (show_axis)
      geom_path(
        data = axisLines,
        aes(x, y, group = seq_id),
        color = "black",
        linewidth = 0.3,
        inherit.aes = FALSE
      )
    } +
    { if (show_axis)
      geom_segment(
        data = axisTicks,
        aes(x = x0, y = y0, xend = x1, yend = y1),
        color = "black",
        linewidth = 0.3,
        inherit.aes = FALSE
      )
    } +
    { if (show_axis)
      geom_text(
        data = subset(axisTicks, !is.na(label)),
        aes(x = label_x, y = label_y, label = label, size = size),
        inherit.aes = FALSE,
        color = "black"
      )
    } +
    scale_size_identity() +
    scale_color_manual(
      name = "Seq ID",
      values = seq_colors,
      labels = seq_labels,
      guide = if (show_legend) guide_legend(order = 2) else "none"
    ) +
    theme_void() +
    coord_equal(clip = "off") +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.margin = margin(t = 0,r = 0,b = 0,l = 0),
      legend.box.spacing = unit(10,"mm"),
      legend.spacing = unit(5,"mm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 12, face = "bold"),  
    )
  
  if (debug) cat("[DEBUG] ggchord函数执行完成，返回ggplot对象\n")
  return(p)
}



# 示例使用方法
# 1. 读取数据
#读取序列长度数据
seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)

# 读取基因注释数据
gene_track <- read.delim("gene_track.tsv", sep = "\t", stringsAsFactors = FALSE) |> dplyr::slice_max(order_by = end-start, n = 5, by = seq_id)

# 读取并处理BLAST数据
read_blast <- function(file) {
  df <- read.delim(file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, comment.char = "#")
  colnames(df) <- c("qaccver","saccver","pident","length","mismatches","gapopen",
                    "qstart","qend","sstart","send","evalue","bitscore",
                    "qcovs","qlen","slen","sstrand","stitle")
  df
}
blast_files <- list.files(path = ".", pattern = "*.o7", full.names = TRUE)
all_blast <- do.call(rbind, lapply(blast_files, read_blast))
ribbon_data <- subset(all_blast, length >= 100)

## 2. 调用ggchord函数（示例）
p_final <- ggchord(
  seq_data = seq_data,
  ribbon_data = ribbon_data,
  gene_track = gene_track,  # 传入基因注释数据
  title = "Multi-sequence Chord Diagram with Gene Annotations",
  seq_gap = c(.1,.05,.03,.09),
  seq_radius = c(4,3,2,1),
  seq_orientation = c(-1, 1, -1, -1),
  gene_offset = list(.15, 
                     .25, 
                     .4, 
                     c("+"=.3,"-"=-.1)),
  gene_width  = .08,
  gene_label_show = F,
  ribbon_gap = c(.4,.2,.2,.2),
  ribbon_color_scheme = "pident",
  axis_gap = c(-.1,.1,.1,.1),
  axis_tick_major_number = 5,
  axis_tick_major_length = 0.03,
  axis_tick_minor_number = 5,
  axis_tick_minor_length = 0.01,
  axis_label_size = 2,
  axis_label_offset = .1,
  rotation = 15,
  debug = TRUE
)
print(p_final)
