#' ggchord: 基于ggplot2的多序列比对弦图工具
#'
#' 该函数用于绘制包含多个序列的弦图，可展示序列间的比对关系和基因注释信息。
#' 支持自定义序列排列、颜色、连接带样式、基因箭头样式等多种参数，适合用于基因组比对可视化。
#'
#' @param seq_data data.frame/tibble，包含序列信息，必须包含列：
#' - seq_id: 序列唯一标识
#' - length: 序列长度
#' @param ribbon_data data.frame/tibble，包含BLAST比对结果，可选，若提供则必须包含列：
#' - qaccver: 查询序列ID
#' - saccver: 目标序列ID
#' - length: 比对长度
#' - pident: 序列相似度百分比
#' - qstart: 查询序列起始位置
#' - qend: 查询序列结束位置
#' - sstart: 目标序列起始位置
#' - send: 目标序列结束位置
#' @param gene_track data.frame/tibble，包含基因注释信息，可选，若提供则必须包含列：
#' - seq_id: 序列唯一标识
#' - start: 基因起始位置
#' - end: 基因结束位置
#' - strand: 链方向（+或-）
#' - anno: 基因注释
#' @param title 字符串，图形主标题，默认NULL
#' @param seq_order 字符向量，可选，指定序列绘制顺序；若为NULL则使用 seq_data 中的顺序
#' @param seq_labels 字符向量或命名向量，可选，序列标签；若为NULL则使用 seq_id
#' @param seq_orientation 数值向量或单值，可选，每条序列的方向：1（正向）或 -1（反向）；默认正向
#' @param seq_gap 数值或向量，长度与序列数一致，定义每条序列头部到下一条序列尾部的弧度比例 [0,0.5)，默认0.03
#' @param seq_radius 数值或向量，序列圆弧半径，支持单值或与序列数相同的向量，默认1.0
#' @param seq_curvature 数值或向量，序列圆弧的弯曲程度：1为标准圆弧，0为直线，>1更弯曲，默认1.0
#' @param seq_colors 颜色向量或命名向量，定义各序列圆弧颜色；若为NULL则基于 RColorBrewer Set1 自动生成
#' @param gene_offset 数值、向量或列表，基因箭头与序列圆弧之间的径向偏移距离。支持：
#' - 单值：所有序列的所有链使用相同偏移
#' - 向量：长度与序列数一致，每个序列的所有链使用相同偏移
#' - 列表：命名列表，每个元素对应一个序列，元素可为单值（该序列所有链）或包含"+"和"-"的命名向量（区分链），默认0.03
#' @param gene_width 数值或向量，基因箭头宽度，默认0.1
#' @param gene_label_show 逻辑值，是否显示基因标签，默认FALSE
#' @param gene_label_rotation 数值、向量或列表，基因标签的旋转角度（度），支持与gene_offset相同的参数格式，默认0
#' @param gene_label_size 数值，基因注释文字大小，默认2.5
#' @param gene_label_radial_offset 数值、向量或列表，基因标签相对于箭头的径向偏移（正值向外，负值向内），支持与gene_offset相同的参数格式，默认0
#' @param gene_label_circum_offset 数值、向量或列表，基因标签沿序列周向的偏移比例（相对于基因长度），支持与gene_offset相同的参数格式，默认0
#' @param gene_label_circum_limit 逻辑值、向量或列表，是否限制周向偏移不超过基因长度的一半，支持与gene_offset相同的参数格式，默认TRUE
#' @param gene_color_scheme 字符，指定基因颜色方案，可选"strand"（按链方向）或"manual"（手动指定），默认"strand"
#' @param gene_colors 颜色向量，用于指定基因箭头的填充色，具体行为取决于gene_color_scheme：
#' - "strand"模式：支持命名向量（仅"+"/"-"）、非命名向量（先"+"后"-"）或单值（正负链同色），缺省则"+"为红色、"-"为蓝色
#' - "manual"模式：支持命名向量（对应anno）、非命名向量（截断多余，补齐不足），缺省则使用自动生成颜色
#' @param gene_order 字符向量，可选，指定基因在图例中的显示顺序；若为NULL则使用基因在数据中出现的顺序
#' @param ribbon_color_scheme 字符，连接带配色方案，可选 "single"、"query" 或 "pident"，默认"pident"
#' @param ribbon_colors 连接带颜色参数：
#' - single: 单一颜色（单值或向量取第一个）
#' - query: 按查询序列映射颜色（命名或非命名向量或单值）
#' - pident: 渐变色阶向量，用于按相似度百分比生成渐变，默认使用蓝到黄的渐变色
#' @param ribbon_alpha 数值，连接带透明度 [0,1]，默认0.35
#' @param ribbon_ctrl_point 向量或列表，贝塞尔控制点，用于调整连接带形状：
#' - 向量：长度为2（单控制点）或4（c1x,c1y,c2x,c2y，双控制点）
#' - 列表：每个元素为包含1-2个控制点的子列表，默认自动计算
#' @param ribbon_gap 数值或向量，序列圆弧与连接带之间的径向距离，默认0.15
#' @param axis_gap 数值或向量，坐标轴与序列圆弧之间径向距离，支持负值，默认0.04
#' @param axis_tick_major_number 整数或向量，每条序列主刻度数，默认5
#' @param axis_tick_major_length 数值或向量，主刻度线长度比例，默认0.02
#' @param axis_tick_minor_number 整数或向量，每两个主刻度之间的次刻度数，默认4
#' @param axis_tick_minor_length 数值或向量，次刻度线长度比例，默认0.01
#' @param axis_label_size 数值或向量，坐标轴刻度文字大小，默认3
#' @param axis_label_offset 数值或向量，坐标轴标签相对于刻度线的偏移比例，默认0
#' @param axis_label_orientation 字符、数值或向量，坐标轴标签方向：
#' - "horizontal"：水平方向
#' - 数值：旋转角度（度）
#' - 向量：长度与序列数一致或命名向量，默认"horizontal"
#' @param rotation 数值，整个图形的旋转角度（度），默认45
#' @param panel_margin 列表，图形边缘留白（t=上, r=右, b=下, l=左），默认list(t=0,r=0,b=0,l=0)
#' @param show_legend 逻辑值，是否显示图例，默认TRUE
#' @param show_axis 逻辑值，是否显示坐标轴及刻度，默认TRUE
#' @param debug 逻辑值，是否输出调试信息，默认FALSE
#' @return ggplot2图形对象，可进一步使用ggplot2函数调整样式


# 加载所需包
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(ggnewscale)
library(grid)


# 处理缺失值
`%||%` <- function (x, y) 
{
  if (is.null(x)) y else x
}

process_panel_margin <- function(arg_list) {
  # 初始化结果列表，设置默认值为0
  result <- list(t = 0, r = 0, b = 0, l = 0)
  
  # 检查输入是否为列表
  if (!is.list(arg_list)) {
    # 处理单值输入的情况
    if (is.numeric(arg_list) && length(arg_list) == 1) {
      value <- arg_list
      result <- list(t = value, r = value, b = value, l = value)
      return(result)
    } else {
      warning("输入不是有效列表或单个数值，将使用默认值")
      return(result)
    }
  }
  
  # 处理空列表
  if (length(arg_list) == 0) {
    return(result)
  }
  
  # 处理命名列表
  if (!is.null(names(arg_list)) && all(names(arg_list) != "")) {
    valid_names <- c("t", "r", "b", "l")
    
    # 遍历输入列表的每个元素
    for (name in names(arg_list)) {
      if (name %in% valid_names) {
        # 检查值是否为数值
        if (is.numeric(arg_list[[name]]) && length(arg_list[[name]]) == 1) {
          result[[name]] <- arg_list[[name]]
        } else {
          warning(paste("参数", name, "不是单个数值，将使用默认值0"))
        }
      } else {
        warning(paste("未知参数", name, "将被忽略"))
      }
    }
  } 
  # 处理非命名列表
  else {
    param_order <- c("t", "r", "b", "l")
    num_args <- length(arg_list)
    
    # 处理单元素非命名列表的情况
    if (num_args == 1 && is.numeric(arg_list[[1]])) {
      value <- arg_list[[1]]
      result <- list(t = value, r = value, b = value, l = value)
      return(result)
    }
    
    # 按新顺序分配值
    for (i in 1:min(num_args, length(param_order))) {
      if (is.numeric(arg_list[[i]]) && length(arg_list[[i]]) == 1) {
        result[[param_order[i]]] <- arg_list[[i]]
      } else {
        warning(paste("位置", i, "的参数不是单个数值，将使用默认值0"))
      }
    }
  }
  
  return(result)
}

# 获取图像四边极点
get_plot_extremes <- function(allRibbon=NULL, seqArcs=NULL, axisLines=NULL, axisTicks=NULL, gene_arrows=NULL, gene_polys = NULL, show_axis = F) {
  # 初始化存储x和y坐标的向量
  x_coords <- numeric(0)
  y_coords <- numeric(0)
  
  # 1. 处理连接带（allRibbon）
  if (!is.null(allRibbon) && nrow(allRibbon) > 0) {
    x_coords <- c(x_coords, allRibbon$x)
    y_coords <- c(y_coords, allRibbon$y)
  }
  # 2. 处理序列弧线（seqArcs，列表转数据框）
  if (!is.null(seqArcs) && length(seqArcs) > 0) {
    seq_df <- do.call(rbind, seqArcs)
    if (nrow(seq_df) > 0) {
      x_coords <- c(x_coords, seq_df$x)
      y_coords <- c(y_coords, seq_df$y)
    }
  }
  # 3. 处理基因标签（gene_arrows）
  if (!is.null(gene_arrows) && nrow(gene_arrows) > 0) {
    x_coords <- c(x_coords, gene_arrows$text_x)
    y_coords <- c(y_coords, gene_arrows$text_y)
  }
  # 4. 处理轴线（axisLines）
  if (show_axis && !is.null(axisLines) && nrow(axisLines) > 0) {
    x_coords <- c(x_coords, axisLines$x)
    y_coords <- c(y_coords, axisLines$y)
  }
  # 5. 处理刻度（axisTicks）
  if (show_axis && !is.null(axisTicks) && nrow(axisTicks) > 0) {
    # x坐标：x0、x1、label_x
    x_coords <- c(x_coords, axisTicks$x0, axisTicks$x1, axisTicks$label_x)
    # y坐标：y0、y1、label_y
    y_coords <- c(y_coords, axisTicks$y0, axisTicks$y1, axisTicks$label_y)
  }
  # 6. 处理基因箭头多边形（gene_polys）：包含x、y坐标（多边形顶点）
  if (!is.null(gene_polys) && nrow(gene_polys) > 0) {
    x_coords <- c(x_coords, gene_polys$x)
    y_coords <- c(y_coords, gene_polys$y)
  }
  
  # 过滤缺失值（NA）
  x_coords <- x_coords[!is.na(x_coords)]
  y_coords <- y_coords[!is.na(y_coords)]
  
  # 计算极值（返回列表，包含上下左右极值）
  list(
    x_min = min(x_coords), # 左极值
    x_max = max(x_coords), # 右极值
    y_min = min(y_coords), # 下极值
    y_max = max(y_coords) # 上极值
  )
}

# 定义自定义的 key 绘制函数
draw_key_gene_arrow <- function(data, params, size) {
  # 五个顶点：矩形左下、矩形左上、矩形右上、箭头尖、矩形右下
  x_pts <- unit(c(0.1, 0.1, 0.6, 0.9, 0.6), "npc")
  y_pts <- unit(c(0.2, 0.8, 0.8, 0.5, 0.2), "npc")
  
  polygonGrob(
    x = x_pts, y = y_pts,
    gp = gpar(
      fill = alpha(data$fill %||% "grey", data$alpha %||% 1),
      col = data$colour %||% "black",
      lwd = data$size %||% 0.5 * .pt
    )
  )
}

# 辅助函数：生成坐标轴主刻度断点
breakPointsFunc <- function(max_value, n = 5, tol = 0.5) {
  if (max_value <= 0) return(c(0, max_value))
  
  # 1. 用 pretty() 生成大致 n 个刻度
  ticks <- pretty(c(0, max_value), n = n)
  ticks <- ticks[ticks >= 0 & ticks <= max_value]
  ticks <- sort(unique(c(0, ticks, max_value)))
  
  # 2. 如果最后一段太小（小于其他分段中位数 * tol），就去掉 penultimate
  if (length(ticks) >= 3) {
    d <- diff(ticks)
    # 其它分段（不含最后一段）的中位数
    med <- median(d[-length(d)])
    # 如果最后一段过小，就剔除 penultimate 点
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

# 处理坐标轴文字方向参数
process_axis_orientation <- function(param, seqs) {
  n <- length(seqs)
  
  # 处理单值"horizontal"
  if (is.character(param) && length(param) == 1 && tolower(param) == "horizontal") {
    return(setNames(rep("horizontal", n), seqs))
  }
  
  # 处理单值数值
  if (is.numeric(param) && length(param) == 1) {
    return(setNames(rep(param, n), seqs))
  }
  
  # 处理向量输入（包括混合类型）
  if (is.vector(param) && length(param) == n) {
    result <- character(n)
    names(result) <- seqs
    
    for (i in seq_along(param)) {
      val <- param[i]
      seq_id <- seqs[i]
      
      if (is.character(val) && tolower(val) == "horizontal") {
        result[seq_id] <- "horizontal"
      } else if (is.numeric(val)) {
        result[seq_id] <- as.character(val)
      } else {
        stop(paste("axis_label_orientation 参数的第", i, "个元素格式不正确，必须是数值或'horizontal'"))
      }
    }
    
    return(result)
  }
  
  # 处理命名向量
  if (!is.null(names(param)) && all(names(param) %in% seqs)) {
    result <- setNames(rep("0", n), seqs)
    for (id in names(param)) {
      val <- param[id]
      if (is.character(val) && tolower(val) == "horizontal") {
        result[id] <- "horizontal"
      } else if (is.numeric(val)) {
        result[id] <- as.character(val)
      } else {
        stop(paste("axis_label_orientation 中", id, "的格式不正确，必须是数值或'horizontal'"))
      }
    }
    return(result)
  }
  
  stop("axis_label_orientation 参数格式不正确，请提供：\n",
       "- 'horizontal'（默认）\n",
       "- 单值数值\n",
       "- 长度与序列数相同的向量（可混合数值和'horizontal'）\n",
       "- 命名向量（名称对应序列ID）")
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


# 修改辅助函数以支持基因标签偏移参数的处理
process_gene_param <- function(param, seqs, param_name, default_value, is_logical = FALSE) {
  n <- length(seqs)
  # 初始化结果列表：每个序列对应一个包含"+"和"-"的向量
  result <- setNames(lapply(seqs, function(id) {
    if (is_logical) {
      c("+" = default_value, "-" = default_value) # 逻辑型参数
    } else {
      c("+" = default_value, "-" = default_value) # 数值型参数
    }
  }), seqs)
  
  if (is.null(param)) {
    return(result)
  }
  
  # 1. 单值情况（所有序列、所有链使用相同值）
  if (length(param) == 1 && !is.list(param)) {
    # 验证参数类型
    if (is_logical && !is.logical(param)) {
      stop(paste(param_name, "为逻辑型参数，必须传入TRUE/FALSE"))
    }
    if (!is_logical && !is.numeric(param)) {
      stop(paste(param_name, "为数值型参数，必须传入数值"))
    }
    val <- param
    # 关键修改：用unname(val)去除名称
    return(setNames(lapply(seqs, function(id) c("+" = unname(val), "-" = unname(val))), seqs))
  }
  
  # 2. 向量情况（非列表）
  if (is.vector(param) && !is.list(param)) {
    # 2.1 命名向量：验证序列名称
    if (!is.null(names(param))) {
      unknown <- setdiff(names(param), seqs)
      if (length(unknown) > 0) {
        warning(paste(param_name, "包含不存在的序列ID:", paste(unknown, collapse = ", ")))
      }
      valid_names <- intersect(names(param), seqs)
      for (id in valid_names) {
        val <- param[id]
        # 验证类型
        if (is_logical && !is.logical(val)) {
          stop(paste(param_name, "的命名向量元素必须为TRUE/FALSE"))
        }
        if (!is_logical && !is.numeric(val)) {
          stop(paste(param_name, "的命名向量元素必须为数值"))
        }
        # 关键修改：用unname(val)去除名称
        result[[id]] <- c("+" = unname(val), "-" = unname(val)) # 向量不区分链，正负链用相同值
      }
      return(result)
    }
    
    # 2.2 未命名向量：验证长度匹配
    if (length(param) != n) {
      stop(paste(param_name, "未命名向量长度必须与序列数相同（当前序列数:", n, "）"))
    }
    for (i in seq_along(seqs)) {
      val <- param[i]
      # 验证类型
      if (is_logical && !is.logical(val)) {
        stop(paste(param_name, "的未命名向量元素必须为TRUE/FALSE"))
      }
      if (!is_logical && !is.numeric(val)) {
        stop(paste(param_name, "的未命名向量元素必须为数值"))
      }
      # 关键修改：用unname(val)去除名称
      result[[seqs[i]]] <- c("+" = unname(val), "-" = unname(val)) # 向量不区分链
    }
    return(result)
  }
  
  # 3. 列表情况（保持不变，因问题不涉及列表）
  if (is.list(param)) {
    # 3.1 单值列表（所有序列的正负链分别设置相同值）
    if (length(param) == 1 && is.null(names(param))) {
      elem <- param[[1]]
      # 检查列表元素是否包含"+"和"-"
      if (!all(c("+", "-") %in% names(elem))) {
        stop(paste(param_name, "单值列表元素必须是包含'+'和'-'的命名向量"))
      }
      # 验证类型
      if (is_logical && (!is.logical(elem["+"]) || !is.logical(elem["-"]))) {
        stop(paste(param_name, "列表元素必须为逻辑值（TRUE/FALSE）"))
      }
      if (!is_logical && (!is.numeric(elem["+"]) || !is.numeric(elem["-"]))) {
        stop(paste(param_name, "列表元素必须为数值"))
      }
      return(setNames(lapply(seqs, function(id) elem[c("+", "-")]), seqs))
    }
    
    # 3.2 命名列表（指定序列的正负链值）
    if (!is.null(names(param))) {
      unknown <- setdiff(names(param), seqs)
      if (length(unknown) > 0) {
        warning(paste(param_name, "列表包含不存在的序列ID:", paste(unknown, collapse = ", ")))
      }
      valid_names <- intersect(names(param), seqs)
      for (id in valid_names) {
        elem <- param[[id]]
        if (!all(c("+", "-") %in% names(elem))) {
          stop(paste(param_name, "列表元素", id, "必须包含'+'和'-'的命名向量"))
        }
        # 验证类型
        if (is_logical && (!is.logical(elem["+"]) || !is.logical(elem["-"]))) {
          stop(paste(param_name, "列表元素", id, "必须为逻辑值（TRUE/FALSE）"))
        }
        if (!is_logical && (!is.numeric(elem["+"]) || !is.numeric(elem["-"]))) {
          stop(paste(param_name, "列表元素", id, "必须为数值"))
        }
        result[[id]] <- elem[c("+", "-")]
      }
      return(result)
    }
    
    # 3.3 未命名列表（按序列顺序指定正负链值）
    if (length(param) != n) {
      stop(paste(param_name, "未命名列表长度必须与序列数相同（当前序列数:", n, "）"))
    }
    for (i in seq_along(seqs)) {
      elem <- param[[i]]
      if (!all(c("+", "-") %in% names(elem))) {
        stop(paste(param_name, "未命名列表第", i, "个元素必须包含'+'和'-'的命名向量"))
      }
      # 验证类型
      if (is_logical && (!is.logical(elem["+"]) || !is.logical(elem["-"]))) {
        stop(paste(param_name, "未命名列表第", i, "个元素必须为逻辑值（TRUE/FALSE）"))
      }
      if (!is_logical && (!is.numeric(elem["+"]) || !is.numeric(elem["-"]))) {
        stop(paste(param_name, "未命名列表第", i, "个元素必须为数值"))
      }
      result[[seqs[i]]] <- elem[c("+", "-")]
    }
    return(result)
  }
  
  # 无效格式
  stop(paste(param_name, "格式错误，支持的传参方式：\n",
             "1. 单值（所有序列/链共用）\n",
             "2. 命名向量（名称为序列ID，所有链共用）\n",
             "3. 未命名向量（长度与序列数相同，所有链共用）\n",
             "4. 单值列表（元素为含'+'/'-'的命名向量，所有序列共用）\n",
             "5. 命名列表（名称为序列ID，元素为含'+'/'-'的命名向量）\n",
             "6. 未命名列表（长度与序列数相同，元素为含'+'/'-'的命名向量）"))
}


# 处理gene_colors参数的辅助函数（strand模式）
process_strand_colors <- function(gene_colors) {
  # 默认值
  default <- c("+" = "#E41A1C", "-" = "#377EB8")
  if (is.null(gene_colors)) {
    return(default)
  }
  # 命名向量处理
  if (!is.null(names(gene_colors))) {
    if (!all(names(gene_colors) %in% c("+", "-"))) {
      stop("'strand'模式下，gene_colors命名向量只能包含'+'和'-'")
    }
    res <- default
    res[names(gene_colors)] <- gene_colors
    return(res)
  } else {
    # 非命名向量处理
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
  # 确定最终的基因顺序
  if (!is.null(gene_order)) {
    # 验证gene_order中的所有元素是否都存在于unique_anno中
    unknown <- setdiff(gene_order, unique_anno)
    if (length(unknown) > 0) {
      stop("'gene_order' 包含未知基因注释: ", paste(unknown, collapse = ","))
    }
    # 确保包含所有unique_anno，按gene_order排序，未在gene_order中的放在最后
    final_order <- c(gene_order, setdiff(unique_anno, gene_order))
  } else {
    final_order <- unique_anno
  }
  n_anno <- length(final_order)
  
  # 生成默认颜色（用于补充）
  if (n_anno <= 9) {
    default_pal <- brewer.pal(n_anno, "Set1")
  } else {
    default_pal <- colorRampPalette(brewer.pal(9, "Set1"))(n_anno)
  }
  names(default_pal) <- final_order
  
  if (is.null(gene_colors)) {
    return(default_pal)
  }
  
  # 命名向量处理（匹配anno）
  if (!is.null(names(gene_colors))) {
    unknown <- setdiff(names(gene_colors), unique_anno)
    if (length(unknown) > 0) {
      warning("'manual'模式下，gene_colors包含未知注释: ", paste(unknown, collapse = ","))
    }
    res <- default_pal
    # 只更新存在于final_order中的颜色
    common_names <- intersect(names(gene_colors), final_order)
    res[common_names] <- gene_colors[common_names]
    return(res)
  } else {
    # 非命名向量处理（先用用户提供的颜色，不足部分用默认颜色补充）
    len <- length(gene_colors)
    res <- character(n_anno)
    
    # 填充用户提供的颜色
    if (len >= 1) {
      res[1:min(len, n_anno)] <- gene_colors[1:min(len, n_anno)]
    }
    
    # 用默认颜色补充剩余部分
    if (len < n_anno) {
      res[(len + 1):n_anno] <- default_pal[(len + 1):n_anno]
    }
    
    names(res) <- final_order
    return(res)
  }
}


# 辅助函数：生成弯曲的序列路径（增加稳定性处理）
generate_curvature_path <- function(start_angle, end_angle, radius, curvature, n_points = 100) {
  # 限制曲率范围
  curvature <- max(-0.99, min(10, curvature))
  
  # 计算中心点和总角度
  center_angle <- (start_angle + end_angle) / 2
  total_angle <- end_angle - start_angle
  
  # 生成均匀分布的角度
  angles <- seq(start_angle, end_angle, length.out = n_points)
  
  # 确保起点和终点的位置固定
  fixed_start_x <- radius * cos(start_angle)
  fixed_start_y <- radius * sin(start_angle)
  fixed_end_x <- radius * cos(end_angle)
  fixed_end_y <- radius * sin(end_angle)
  
  # 根据曲率调整路径
  if (curvature == 0) {
    # 直线
    t <- seq(0, 1, length.out = n_points)
    x <- fixed_start_x + t * (fixed_end_x - fixed_start_x)
    y <- fixed_start_y + t * (fixed_end_y - fixed_start_y)
  } else if (curvature == 1) {
    # 原始圆形路径
    x <- radius * cos(angles)
    y <- radius * sin(angles)
  } else {
    # 使用贝塞尔曲线控制平滑的弯曲，确保起点和终点固定
    # 控制点位置基于曲率调整
    mid_angle <- center_angle
    mid_radius_factor <- 1 + (curvature - 1) * 0.3 # 减小曲率对控制点的影响
    mid_radius <- radius * mid_radius_factor
    
    # 控制点位置
    control_x <- mid_radius * cos(mid_angle)
    control_y <- mid_radius * sin(mid_angle)
    
    # 使用贝塞尔曲线计算路径点
    x <- numeric(n_points)
    y <- numeric(n_points)
    
    for (i in 1:n_points) {
      t <- (i - 1) / (n_points - 1)
      
      # 二次贝塞尔曲线公式
      x[i] <- (1 - t)^2 * fixed_start_x + 2 * (1 - t) * t * control_x + t^2 * fixed_end_x
      y[i] <- (1 - t)^2 * fixed_start_y + 2 * (1 - t) * t * control_y + t^2 * fixed_end_y
    }
  }
  
  # 确保没有NaN或Inf值
  x[is.na(x) | is.infinite(x)] <- 0
  y[is.na(y) | is.infinite(y)] <- 0
  
  data.frame(x = x, y = y)
}

# ggplot2绘图函数
chordPlotFunc <- function(allRibbon,ribbon_alpha,ribbon_color_scheme,ribbon_colors,show_legend,gene_polys,gene_pal,gene_color_scheme,seqArcs,gene_arrows,show_axis,axisLines,axisTicks,axisLabelOrientation,seq_colors,seq_labels,extremes,panel_margin,title) {
  ggplot() +
    # 1. 绘制连接带并设置第一个fill尺度（仅当有有效连接带数据时）
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
    # 连接带颜色尺度（仅当有有效连接带数据时）
    { if (!is.null(allRibbon)) {
      if (ribbon_color_scheme == "pident") {
        scale_fill_stepsn(
          name = "Identity(%)",
          colours = ribbon_colors,
          limits = c(0, 100), 
          breaks = c(0,50,80,90,95,100),
          guide = if (show_legend) guide_colorbar(theme = theme(legend.title.position = "top",legend.key.height = unit(120,"mm")),order = 1,position = "left") else "none"
        )
      } else {
        scale_fill_identity(guide = "none")
      }
    }
    } +
    # 2. 重置fill尺度（核心：使用ggnewscale，仅当有基因数据时）
    { if (nrow(gene_polys) > 0) new_scale_fill() } +
    
    # 3. 绘制基因箭头并设置第二个fill尺度（根据模式映射strand或anno）
    { if (nrow(gene_polys) > 0)
      geom_polygon(
        data = gene_polys,
        aes(x = x, y = y, group = group, 
            fill = if (gene_color_scheme == "strand") strand else anno), # 动态映射填充变量
        color = "black", key_glyph = draw_key_gene_arrow 
      )
    } +
    # 基因箭头颜色尺度（仅当有基因数据时）
    { if (nrow(gene_polys) > 0)
      scale_fill_manual(
        name = if (gene_color_scheme == "strand") "Strand" else "Gene Annotation",
        breaks = if (gene_color_scheme == "strand") c("+","-") else final_gene_order,
        values = gene_pal,
        guide = if (show_legend) guide_legend(order = 3) else "none"
      )
    } +
    # 4. 绘制其他元素
    # 序列弧线
    geom_path(
      data = do.call(rbind, seqArcs),
      aes(x, y, group = seq_id, color = seq_id),
      linewidth = 1.2,
      arrow = arrow(angle = 25, length = unit(2, "mm"), type = "closed")
    ) +
    # 基因标签（仅当需要显示且有基因数据时）
    { if (nrow(gene_arrows) > 0 && gene_label_show)
      geom_text(
        data = gene_arrows,
        aes(x = text_x, y = text_y, label = text, angle = text_angle, 
            hjust = hjust, vjust = vjust), # 使用计算好的对齐参数
        size = gene_label_size,
        color = "black",
        inherit.aes = FALSE
      )
    } +
    # 坐标轴线（仅当show_axis为TRUE时绘制）
    { if (show_axis)
      geom_path(
        data = axisLines,
        aes(x, y, group = seq_id),
        color = "black",
        linewidth = 0.3,
        inherit.aes = FALSE
      )
    } +
    # 刻度线（仅当show_axis为TRUE时绘制）
    { if (show_axis)
      geom_segment(
        data = axisTicks,
        aes(x = x0, y = y0, xend = x1, yend = y1),
        color = "black",
        linewidth = 0.3,
        inherit.aes = FALSE
      )
    } +
    # 刻度标签（仅当show_axis为TRUE时绘制）
    { if (show_axis) {
      # 先创建有标签的刻度数据子集
      label_data <- subset(axisTicks, !is.na(label))
      # 计算每个标签的角度，确保长度匹配
      label_angles <- ifelse(
        axisLabelOrientation[label_data$seq_id] == "horizontal",
        0,
        atan2(label_data$label_y - label_data$y0, label_data$label_x - label_data$x0) * 180 / pi + 90 +
          as.numeric(axisLabelOrientation[label_data$seq_id])
      )
      
      geom_text(
        data = label_data,
        aes(x = label_x, y = label_y, label = label, size = size, group = seq_id),
        inherit.aes = FALSE,
        color = "black",
        angle = label_angles # 使用与数据长度匹配的角度向量
      )
    }
    } +
    scale_size_identity() +
    # 序列颜色
    scale_color_manual(
      name = "Seq ID",
      values = seq_colors,
      labels = seq_labels,
      guide = if (show_legend) guide_legend(order = 2) else "none"
    ) +
    
    # 主题设置
    coord_equal(clip = "off",xlim = c(extremes$x_min - process_panel_margin(panel_margin)$l, extremes$x_max + process_panel_margin(panel_margin)$r), ylim = c(extremes$y_min - process_panel_margin(panel_margin)$b, extremes$y_max + process_panel_margin(panel_margin)$t)) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = .5, size = 20, face = "bold"),
      plot.margin = margin(t = 0,r = 0,b = 0,l = 0),
      legend.background = element_blank(),
      legend.box.spacing = unit(10,"mm"),
      legend.spacing = unit(5,"mm"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, face = "bold"),
      axis.title = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.background = element_blank()
    )
}

# 主函数：绘制多序列弦图
ggchord <- function(
    seq_data,
    ribbon_data = NULL,
    gene_track = NULL,
    title = NULL,
    seq_order = NULL,
    seq_labels = NULL,
    seq_orientation = NULL,
    seq_gap = .03,
    seq_radius = 1.0,
    seq_colors = NULL,
    seq_curvature = 1.0,
    gene_offset = .1,
    gene_width = 0.05,
    gene_label_show = FALSE,
    gene_label_rotation = 0,
    gene_label_size = 2.5,
    gene_label_radial_offset = 0,
    gene_label_circum_offset = 0,
    gene_label_circum_limit = TRUE,
    gene_color_scheme = c("strand", "manual"),
    gene_colors = NULL,
    gene_order = NULL,
    ribbon_color_scheme = c("pident", "query", "single"),
    ribbon_colors = NULL,
    ribbon_alpha = 0.35,
    ribbon_ctrl_point = c(0,0),
    ribbon_gap = .15,
    axis_gap = .04,
    axis_tick_major_number = 5,
    axis_tick_major_length = 0.02,
    axis_tick_minor_number = 4,
    axis_tick_minor_length = 0.01,
    axis_label_size = 3,
    axis_label_offset = 1.5,
    axis_label_orientation = "horizontal",
    rotation = 45,
    panel_margin = 0,
    show_legend = TRUE,
    show_axis = TRUE,
    debug = FALSE
) {
  # 检查必要的包是否安装
  if (!"ggplot2" %in% installed.packages()) stop("需要安装 ggplot2 包")
  if (!"RColorBrewer" %in% installed.packages()) stop("需要安装 RColorBrewer 包")
  if (!"grDevices" %in% installed.packages()) stop("需要安装 grDevices 包")
  if (!"ggnewscale" %in% installed.packages()) stop("需要安装 ggnewscale 包以支持多填充色映射")
  
  ribbon_color_scheme <- match.arg(ribbon_color_scheme)
  gene_color_scheme <- match.arg(gene_color_scheme)
  
  # 1. 验证输入数据格式
  # 验证序列数据（必选）
  required_seq_cols <- c("seq_id", "length")
  if (!all(required_seq_cols %in% colnames(seq_data))) {
    stop("seq_data 必须包含以下列: ", paste(required_seq_cols, collapse = ", "))
  }
  if (any(seq_data$length <= 0)) {
    stop("seq_data 中的 length 必须为正数")
  }
  
  # 验证BLAST数据（可选）
  fb_all <- NULL
  if (!is.null(ribbon_data)) {
    required_ribbon_cols <- c("qaccver", "saccver", "length", "pident", 
                              "qstart", "qend", "sstart", "send")
    if (!all(required_ribbon_cols %in% colnames(ribbon_data))) {
      stop("ribbon_data 必须包含以下列: ", paste(required_ribbon_cols, collapse = ", "))
    }
    fb_all <- ribbon_data # 直接使用传入的已处理比对数据
    if (nrow(fb_all) == 0) warning("ribbon_data 中没有有效比对数据")
    if (debug) cat("使用的比对数据行数:", nrow(fb_all), "\n")
  }
  
  # 验证基因注释数据（可选）
  if (!is.null(gene_track)) {
    required_gene_cols <- c("seq_id", "start", "end", "strand", "anno")
    if (!all(required_gene_cols %in% colnames(gene_track))) {
      stop("gene_track 必须包含以下列: ", paste(required_gene_cols, collapse = ", "))
    }
    if (nrow(gene_track) == 0) warning("gene_track 中没有有效基因注释数据")
    if (debug) cat("使用的基因注释数据行数:", nrow(gene_track), "\n")
    
    # 验证gene_order参数
    if (!is.null(gene_order)) {
      unknown_genes <- setdiff(gene_order, unique(gene_track$anno))
      if (length(unknown_genes) > 0) {
        stop("gene_order 包含未知基因注释: ", paste(unknown_genes, collapse = ", "))
      }
    }
  }
  
  # 2. 处理序列信息
  seqs <- seq_data$seq_id
  lens <- setNames(seq_data$length, seqs)
  if (length(seqs) < 1) stop("seq_data 中至少需要包含1条序列")
  
  # 处理序列顺序
  if (!is.null(seq_order)) {
    if (!all(seq_order %in% seqs)) {
      stop("seq_order 包含未知序列ID: ", paste(setdiff(seq_order, seqs), collapse = ", "))
    }
    seqs <- seq_order
    lens <- lens[seqs] # 按新顺序重新排列长度
  }
  n <- length(seqs) # 序列数量
  
  # 3. 处理序列相关参数
  seq_labels <- process_sequence_param(seq_labels, seqs, "seq_labels", default_value = seqs)
  seqRadius <- process_sequence_param(seq_radius, seqs, "seq_radius", default_value = 1.0)
  ribbonGap <- process_sequence_param(ribbon_gap, seqs, "ribbon_gap", default_value = 0.1)
  axisGap <- process_sequence_param(axis_gap, seqs, "axis_gap", default_value = 0.05)
  axisMaj <- process_sequence_param(axis_tick_major_number, seqs, "axis_tick_major_number", default_value = 5)
  axisMajLen <- process_sequence_param(axis_tick_major_length, seqs, "axis_tick_major_length", default_value = 0.02)
  axisMin <- process_sequence_param(axis_tick_minor_number, seqs, "axis_tick_minor_number", default_value = 4)
  axisMinLen <- process_sequence_param(axis_tick_minor_length, seqs, "axis_tick_minor_length", default_value = 0.01)
  axisLabelOrientation <- process_axis_orientation(axis_label_orientation, seqs)
  labelSize <- process_sequence_param(axis_label_size, seqs, "axis_label_size", default_value = 3)
  labelOffset <- process_sequence_param(axis_label_offset, seqs, "axis_label_offset", default_value = 0)
  orientation <- process_sequence_param(seq_orientation, seqs, "seq_orientation", default_value = 1)
  rot_rad <- rotation * pi / 180 # 转换为弧度
  # 处理基因偏移参数
  geneGap <- process_gene_param(
    param = gene_offset,
    seqs = seqs,
    param_name = "gene_offset",
    default_value = 0.03,
    is_logical = FALSE
  )
  geneWidth <- process_gene_param(
    param = gene_width,
    seqs = seqs,
    param_name = "gene_width",
    default_value = 0.1,
    is_logical = FALSE
  )
  # 处理基因标签偏移参数
  geneLabelRadialOffset <- process_gene_param(
    param = gene_label_radial_offset,
    seqs = seqs,
    param_name = "gene_label_radial_offset",
    default_value = 0,
    is_logical = FALSE
  )
  geneLabelCircumOffset <- process_gene_param(
    param = gene_label_circum_offset,
    seqs = seqs,
    param_name = "gene_label_circum_offset",
    default_value = 0,
    is_logical = FALSE
  )
  geneLabelCircumLimit <- process_gene_param(
    param = gene_label_circum_limit,
    seqs = seqs,
    param_name = "gene_label_circum_limit",
    default_value = TRUE,
    is_logical = TRUE
  )
  # 处理基因标签旋转角度参数
  geneLabelRotation <- process_gene_param(
    param = gene_label_rotation,
    seqs = seqs,
    param_name = "gene_label_rotation",
    default_value = 0,
    is_logical = FALSE
  )
  
  # 处理序列间隙参数
  seq_gap <- process_sequence_param(seq_gap, seqs, "seq_gap")
  if (any(seq_gap < 0 | seq_gap >= 0.5)) {
    stop("seq_gap 必须在 [0, 0.5) 范围内")
  }
  
  # 处理序列弯曲参数
  seq_curvature <- process_sequence_param(seq_curvature, seqs, "seq_curvature", default_value = 1.0)
  
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
  
  # 4. 处理连接带颜色参数（仅当提供ribbon_data时）
  if (!is.null(ribbon_data)) {
    if (is.null(ribbon_colors)) {
      ribbon_colors <- switch(ribbon_color_scheme,
                              single = "steelblue",
                              query = {
                                # 浅色系版本的 seq_colors（与白色按 mix 比例混合）
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
    
    # 验证连接带颜色参数
    if (ribbon_color_scheme == "single") {
      singleCol <- if (length(ribbon_colors) > 1) ribbon_colors[[1]] else ribbon_colors
    } else if (ribbon_color_scheme == "query") {
      queryCols <- process_sequence_param(ribbon_colors, seqs, "ribbon_colors")
    } else if (ribbon_color_scheme == "pident") {
      if (length(ribbon_colors) < 2) {
        stop("pident 模式需指定至少两个颜色用作渐变色阶")
      }
      rampFunc <- colorRampPalette(ribbon_colors)
    }
  }
  
  # 5. 计算基因颜色（基于gene_color_scheme）
  gene_pal <- NULL
  final_gene_order <- NULL
  if (!is.null(gene_track) && nrow(gene_track) > 0) {
    valid_genes <- gene_track[gene_track$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      # 获取唯一的基因注释
      unique_anno <- unique(valid_genes$anno)
      
      # 确定最终的基因顺序
      if (!is.null(gene_order)) {
        final_gene_order <- c(gene_order, setdiff(unique_anno, gene_order))
      } else {
        final_gene_order <- unique_anno
      }
      
      # 根据模式处理颜色
      if (gene_color_scheme == "strand") {
        gene_pal <- process_strand_colors(gene_colors)
      } else if (gene_color_scheme == "manual") {
        gene_pal <- process_manual_colors(gene_colors, unique_anno, gene_order)
      }
    }
  }
  
  # 6. 计算序列弧度和间隙弧度
  total_circ <- 2 * pi # 总圆周弧度
  total_gap_prop <- sum(seq_gap) # 总间隙比例
  
  if (total_gap_prop >= 1) {
    stop("seq_gap 总和不能超过1（无法容纳所有序列）")
  }
  
  seq_total_prop <- 1 - total_gap_prop # 序列总占比
  sum_lens <- sum(lens)
  theta <- (lens / sum_lens) * total_circ * seq_total_prop # 每个序列的弧度
  gap_rads <- total_circ * seq_gap # 每个间隙的弧度（实际角度）
  
  # 7. 计算序列起始和结束角度
  starts <- numeric(n)
  starts[1] <- 0 # 第一个序列从0开始
  
  if (n > 1) {
    for (i in 2:n) {
      starts[i] <- starts[i-1] + theta[i-1] + gap_rads[i-1]
    }
  }
  
  ends <- starts + theta
  names(starts) <- names(ends) <- seqs
  
  # 8. 准备序列外层和内层坐标（修正innerArcs的半径计算）
  nSeg <- 500 # 每个序列的分段数（控制平滑度）
  seqArcs <- lapply(seqs, function(id) {
    # 序列外层圆弧（不变）
    path_data <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id], seq_curvature[id], nSeg
    )
    path_data$seq_id <- id
    if (orientation[id] == -1) {
      path_data <- path_data[nrow(path_data):1, ]
    }
    path_data
  })
  
  innerArcs <- lapply(seqs, function(id) {
    # 连接带所在的内层圆弧（修正：seqRadius + ribbonGap → 正值向外，负值向内）
    # 原理：ribbonGap>0时，半径 = 序列半径 + 间隙（向外）；ribbonGap<0时，半径 = 序列半径 + 负数 = 变小（向内）
    path_data <- generate_curvature_path(
      starts[id], ends[id], seqRadius[id] + ribbonGap[id], seq_curvature[id], nSeg
    )
    if (orientation[id] == -1) {
      path_data <- path_data[nrow(path_data):1, ]
    }
    path_data
  })
  names(innerArcs) <- seqs
  
  
  # —— 新增：对每条序列生成高分辨率参考路径 —— 
  seq_refs <- lapply(seqs, function(id) {
    ref_n <- 2000
    # 基准半径取最外层：序列半径
    r0 <- seqRadius[id]
    path <- generate_curvature_path(starts[id], ends[id], r0, seq_curvature[id], n = ref_n)
    angles <- seq(starts[id], ends[id], length.out = ref_n)
    list(path = path, angles = angles, r0 = r0)
  })
  names(seq_refs) <- seqs
  
  # —— 新增：通用映射函数 —— 
  map_to_curve <- function(angle, radius, ref) {
    # 找到最接近角度的参考点
    idx <- which.min(abs(ref$angles - angle))
    base <- ref$path[idx, ]
    # 计算切线和法线
    if (idx < nrow(ref$path)) {
      dx <- ref$path$x[idx+1] - base$x
      dy <- ref$path$y[idx+1] - base$y
    } else {
      dx <- base$x - ref$path$x[idx-1]
      dy <- base$y - ref$path$y[idx-1]
    }
    norm <- c(-dy, dx)
    norm <- norm / sqrt(sum(norm^2))
    # 沿法线偏移
    offset <- radius - ref$r0
    c(x = base$x + norm[1] * offset,
      y = base$y + norm[2] * offset)
  }
  
  
  # —— 新增：利用 map_to_curve 生成刻度线和标签坐标 —— 
  axisTicks <- do.call(rbind, lapply(seqs, function(id) {
    ref <- seq_refs[[id]]
    # 基准半径：序列半径 + axisGap
    r0 <- ref$r0 - axisGap[id]
    
    # 主次刻度位置
    majors <- breakPointsFunc(lens[id], axisMaj[id])
    minors <- unlist(lapply(seq_len(length(majors)-1), function(i) {
      seq(majors[i], majors[i+1], length.out = axisMin[id] + 2)[-c(1, axisMin[id] + 2)]
    }))
    pts <- data.frame(pos = c(majors, minors),
                      is_major = c(rep(TRUE, length(majors)), rep(FALSE, length(minors))))
    
    # 每个刻度点都做映射
    do.call(rbind, lapply(seq_len(nrow(pts)), function(j) {
      p <- pts[j,]
      frac <- if (orientation[id] == 1) p$pos / lens[id] else 1 - p$pos / lens[id]
      angle <- starts[id] + frac * (ends[id] - starts[id])
      
      # 基准点、刻度头点、标签点
      base <- map_to_curve(angle, r0, ref)
      dir <- if (axisGap[id] >= 0) -1 else 1
      len <- if (p$is_major) axisMajLen[id] else axisMinLen[id]
      tip <- map_to_curve(angle, r0 + len * dir, ref)
      lbl <- map_to_curve(angle, r0 + len * (1.5 + labelOffset[id]) * dir, ref)
      
      data.frame(
        x0 = base[1], y0 = base[2],
        x1 = tip[1], y1 = tip[2],
        label = if (p$is_major) as.character(p$pos) else NA,
        label_x = lbl[1], label_y = lbl[2],
        size = labelSize[id],
        seq_id = id
      )
    }))
  }))
  
  
  # —— 新增：利用 map_to_curve 生成坐标轴整条线 —— 
  axisLines <- do.call(rbind, lapply(seqs, function(id) {
    ref <- seq_refs[[id]]
    # 基准半径：序列半径 + axisGap
    r0 <- ref$r0 - axisGap[id]
    # 把整个角度区间分成 nSeg 段
    angles <- seq(starts[id], ends[id], length.out = nSeg)
    # 对每个角度都做映射
    pts <- t(sapply(angles, function(angle) {
      map_to_curve(angle, r0, ref)
    }))
    data.frame(x = pts[,1], y = pts[,2], seq_id = id)
  }))
  
  # 10. 生成连接带数据（仅当提供ribbon_data时）
  allRibbon <- NULL
  if (!is.null(ribbon_data) && !is.null(fb_all) && nrow(fb_all) > 0) {
    ribbons <- list()
    cntValid <- cntInvalid <- 0
    for (i in seq_len(nrow(fb_all))) {
      row <- fb_all[i,]
      q <- row$qaccver
      s <- row$saccver
      if (q == s || !q %in% seqs || !s %in% seqs) { 
        cntInvalid <- cntInvalid + 1
        next
      }
      
      # 修正查询序列连接带起点/终点坐标（使用 seqRadius + ribbonGap 确保方向正确）
      q_ref <- seq_refs[[q]]
      q_frac_start <- if (orientation[q] == 1) (row$qstart - 1)/lens[q] else 1 - (row$qstart - 1)/lens[q]
      q_angle_start <- starts[q] + q_frac_start * (ends[q] - starts[q])
      # 连接带半径 = 序列半径 + ribbonGap（正值向外，负值向内）
      q_start_coord <- map_to_curve(q_angle_start, seqRadius[q] + ribbonGap[q], q_ref)
      
      q_frac_end <- if (orientation[q] == 1) (row$qend - 1)/lens[q] else 1 - (row$qend - 1)/lens[q]
      q_angle_end <- starts[q] + q_frac_end * (ends[q] - starts[q])
      q_end_coord <- map_to_curve(q_angle_end, seqRadius[q] + ribbonGap[q], q_ref)
      
      # 生成查询序列的比对片段（弯曲路径）
      q_angles <- seq(q_angle_start, q_angle_end, length.out = 50)
      q_coords <- do.call(rbind, lapply(q_angles, function(angle) {
        map_to_curve(angle, seqRadius[q] + ribbonGap[q], q_ref)
      }))
      segQ <- data.frame(x = q_coords[,1], y = q_coords[,2])
      
      # 修正目标序列连接带起点/终点坐标
      s_ref <- seq_refs[[s]]
      s_frac_start <- if (orientation[s] == 1) (row$sstart - 1)/lens[s] else 1 - (row$sstart - 1)/lens[s]
      s_angle_start <- starts[s] + s_frac_start * (ends[s] - starts[s])
      s_start_coord <- map_to_curve(s_angle_start, seqRadius[s] + ribbonGap[s], s_ref)
      
      s_frac_end <- if (orientation[s] == 1) (row$send - 1)/lens[s] else 1 - (row$send - 1)/lens[s]
      s_angle_end <- starts[s] + s_frac_end * (ends[s] - starts[s])
      s_end_coord <- map_to_curve(s_angle_end, seqRadius[s] + ribbonGap[s], s_ref)
      
      # 生成目标序列的比对片段（弯曲路径）
      s_angles <- seq(s_angle_start, s_angle_end, length.out = 50)
      s_coords <- do.call(rbind, lapply(s_angles, function(angle) {
        map_to_curve(angle, seqRadius[s] + ribbonGap[s], s_ref)
      }))
      segS <- data.frame(x = s_coords[,1], y = s_coords[,2])
      
      # 确定贝塞尔曲线控制点（修复ribbon_ctrl_point参数处理）
      if (!is.null(ribbon_ctrl_point)) {
        # 处理列表形式：每个元素对应一个连接带的两个控制点(c1, c2)
        if (is.list(ribbon_ctrl_point)) {
          # 确保列表索引不越界
          cp_idx <- ifelse(i > length(ribbon_ctrl_point), length(ribbon_ctrl_point), i)
          cp <- ribbon_ctrl_point[[cp_idx]]
          # 确保每个控制点元素包含两个点
          if (length(cp) >= 2) {
            c1 <- cp[[1]] # 起点曲线控制点
            c2 <- cp[[2]] # 终点曲线控制点
          } else {
            # 如果元素不足，用第一个点作为默认
            c1 <- c2 <- if (length(cp) == 1) cp[[1]] else c(0, 0)
          }
        } else {
          # 处理向量形式：长度为2时表示单个点，长度为4时表示两个点(c1x,c1y,c2x,c2y)
          if (length(ribbon_ctrl_point) == 2) {
            # 单个控制点应用于所有曲线
            c1 <- c2 <- ribbon_ctrl_point
          } else if (length(ribbon_ctrl_point) == 4) {
            # 分别指定起点和终点控制点
            c1 <- ribbon_ctrl_point[1:2]
            c2 <- ribbon_ctrl_point[3:4]
          } else {
            # 无效长度时使用默认（圆心）
            warning("ribbon_ctrl_point向量长度必须为2或4，已使用默认值")
            c1 <- c2 <- c(0, 0)
          }
        }
      } else {
        # 自动计算适合曲率的控制点（原始逻辑保留）
        mid_angle_q <- (q_angle_start + q_angle_end) / 2
        mid_angle_s <- (s_angle_start + s_angle_end) / 2
        mid_point_q <- map_to_curve(mid_angle_q, seqRadius[q] + ribbonGap[q] * 0.5, q_ref)
        mid_point_s <- map_to_curve(mid_angle_s, seqRadius[s] + ribbonGap[s] * 0.5, s_ref)
        c1 <- colMeans(rbind(mid_point_q, mid_point_s))
        c2 <- c1
      }
      
      # 生成贝塞尔曲线（使用修复后的c1和c2）
      b1 <- bezier_pts(
        as.numeric(segQ[1, ]), # 起点曲线起点（查询序列起点）
        as.numeric(segS[1, ]), # 起点曲线终点（目标序列起点）
        c1, c1, # 使用起点控制点c1
        n = 50
      ) 
      b2 <- bezier_pts(
        as.numeric(segQ[nrow(segQ), ]), # 终点曲线起点（查询序列终点）
        as.numeric(segS[nrow(segS), ]), # 终点曲线终点（目标序列终点）
        c2, c2, # 使用终点控制点c2
        n = 50
      )
      
      # 构建闭合多边形
      poly <- rbind(
        segQ, # 1. 查询序列片段：查询起点 → 查询终点（正向）
        b2, # 2. 终点贝塞尔曲线：查询终点 → 目标终点（正向）
        segS[nrow(segS):1, ], # 3. 目标序列片段：目标终点 → 目标起点（反向，与b2终点衔接）
        b1[nrow(b1):1, ] # 4. 起点贝塞尔曲线：目标起点 → 查询起点（反向，闭合回到起点）
      )
      
      # 设置填充色
      if (ribbon_color_scheme == "pident") {
        poly$pident <- row$pident
      } else {
        poly$fill <- switch(ribbon_color_scheme,
                            single = singleCol,
                            query = queryCols[q])
      }
      poly$group <- i # 分组标识
      ribbons[[length(ribbons) + 1]] <- poly
      cntValid <- cntValid + 1
    }
    
    if (debug) {
      cat("有效连接带:", cntValid, "无效连接带:", cntInvalid, "\n")
    }
    if (cntValid > 0) {
      allRibbon <- do.call(rbind, ribbons)
    } else {
      allRibbon <- NULL
      warning("没有有效连接带可绘制")
    }
  }
  
  
  # 11. 处理基因注释箭头 —— 关键修改基因偏移方向的部分
  gene_polys <- data.frame()
  if (!is.null(gene_track) && nrow(gene_track) > 0) {
    valid_genes <- gene_track[gene_track$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      for (i in seq_len(nrow(valid_genes))) {
        gene <- valid_genes[i, ]
        sid <- gene$seq_id
        strand <- gene$strand
        anno <- gene$anno
        
        # 校验宽度
        width <- geneWidth[[sid]][ strand ]
        if (!is.numeric(width) || width <= 0) width <- 0.1
        
        # 计算起/终角度（保持不变）
        seq_len <- lens[sid]
        sp <- min(gene$start, gene$end); ep <- max(gene$start, gene$end)
        if (ep <= sp) next
        frac_sp <- if (orientation[sid]==1) sp/seq_len else 1-sp/seq_len
        frac_ep <- if (orientation[sid]==1) ep/seq_len else 1-ep/seq_len
        a_start <- starts[sid] + frac_sp*(ends[sid]-starts[sid])
        a_end <- starts[sid] + frac_ep*(ends[sid]-starts[sid])
        if (strand == "-") { tmp <- a_start; a_start <- a_end; a_end <- tmp }
        
        # 生成角度向量和宽度因子（保持不变）
        n_body <- 30; n_head <- 15
        body_ang <- seq(a_start, a_start + 0.6*(a_end - a_start), length.out = n_body)
        head_ang <- seq(tail(body_ang,1), a_end, length.out = n_head)
        angs <- c(body_ang, head_ang)
        total_pt <- length(angs)
        widths <- c(rep(1, n_body), seq(1, 0, length.out = n_head))
        
        # 基准半径计算 —— 核心修改：根据链方向反转gene_offset的作用方向
        # 正向链(+): gene_offset>0时向外偏移（seqRadius - offset）
        # 负向链(-): gene_offset>0时向内偏移（seqRadius + offset）
        if (strand == "+") {
          r0 <- seqRadius[sid] - geneGap[[sid]][strand] # 正向链：减号保持向外
        } else {
          r0 <- seqRadius[sid] + geneGap[[sid]][strand] # 负向链：加号实现向内
        }
        
        # 计算内外半径（保持不变）
        outer_r <- r0 + (width/2)*widths
        inner_r <- r0 - (width/2)*widths
        
        # 原始极坐标表（保持不变）
        orig_ang <- c(angs, rev(angs))
        orig_rad <- c(outer_r, rev(inner_r))
        
        # 映射每一点（保持不变）
        ref <- seq_refs[[sid]]
        pts_list <- mapply(function(ang, rad) {
          map_to_curve(ang, rad, ref)
        }, orig_ang, orig_rad, SIMPLIFY = FALSE)
        
        mapped <- do.call(rbind, pts_list)
        
        # 构建多边形（保持不变）
        gene_poly <- data.frame(
          x = mapped[,1],
          y = mapped[,2],
          group = i,
          anno = anno,
          strand = strand,
          ord = seq_len(2 * total_pt)
        )
        
        gene_polys <- rbind(gene_polys, gene_poly)
      }
    }
  }
  
  
  # 12. 处理基因标签，确保与箭头偏移同步
  gene_arrows <- data.frame()
  if (!is.null(gene_track) && nrow(gene_track) > 0 && gene_label_show) {
    valid_genes <- gene_track[gene_track$seq_id %in% seqs, ]
    if (nrow(valid_genes) > 0) {
      gene_arrows <- do.call(rbind, lapply(seq_len(nrow(valid_genes)), function(i) {
        gene <- valid_genes[i, ]
        sid <- gene$seq_id
        strand <- gene$strand
        seq_len <- lens[sid]
        ref <- seq_refs[[sid]]
        orient <- orientation[sid] # 1=正向，-1=反向
        
        # 1. 计算基因中点相对位置（不变）
        sp <- min(gene$start, gene$end); ep <- max(gene$start, gene$end)
        frac_mid <- (sp + ep) / (2 * seq_len)
        
        # 应用周向偏移（不变）
        circum_ratio <- geneLabelCircumOffset[[sid]][strand]
        if (geneLabelCircumLimit[[sid]][strand]) {
          gene_length_ratio <- (ep - sp) / seq_len
          max_offset_ratio <- gene_length_ratio * 0.5
          circum_ratio <- pmin(max_offset_ratio, pmax(-max_offset_ratio, circum_ratio))
        }
        frac_mid <- frac_mid + circum_ratio
        frac_mid <- pmin(1, pmax(0, frac_mid))
        
        # 反向序列调整中点位置（不变）
        if (orient != 1) frac_mid <- 1 - frac_mid
        
        # 2. 找参考路径索引及切线向量(dx,dy)（不变）
        ref_n <- length(ref$angles)
        idx <- round(frac_mid * (ref_n - 1)) + 1
        idx <- pmin(ref_n, pmax(1, idx))
        if (idx < ref_n) {
          dx <- ref$path$x[idx + 1] - ref$path$x[idx]
          dy <- ref$path$y[idx + 1] - ref$path$y[idx]
        } else {
          dx <- ref$path$x[idx] - ref$path$x[idx - 1]
          dy <- ref$path$y[idx] - ref$path$y[idx - 1]
        }
        dx <- dx * orient # 根据序列方向修正切线方向
        dy <- dy * orient
        
        # 3. 计算箭头中心点（核心修改：与箭头使用相同的基准半径逻辑）
        width <- geneWidth[[sid]][ strand ]
        
        # 与箭头保持一致的基准半径计算
        if (strand == "+") {
          r0 <- seqRadius[sid] - geneGap[[sid]][strand] # 正向链：向外偏移
        } else {
          r0 <- seqRadius[sid] + geneGap[[sid]][strand] # 负向链：向内偏移
        }
        
        # 计算箭头中心点半径（基于箭头的内外半径中间）
        center_r <- r0 # 箭头中心点就是基准半径
        
        # 计算中心点坐标
        center_pt <- map_to_curve(angle = ref$angles[idx], radius = center_r, ref = ref)
        
        # 4. 计算法线向量（确保与箭头偏移方向一致）
        normal_x <- -dy # 基础法线方向（垂直切线）
        normal_y <- dx
        nl <- sqrt(normal_x^2 + normal_y^2)
        if (nl > 0) {
          normal_x <- normal_x / nl # 归一化
          normal_y <- normal_y / nl
        }
        
        # 根据链方向和序列方向调整法线方向
        direction_factor <- ifelse(strand == "+", 1, -1) * orient
        normal_x <- normal_x * direction_factor
        normal_y <- normal_y * direction_factor
        
        # 5. 计算标签位置（与箭头偏移同步）
        text_x <- center_pt[1] - normal_x * geneLabelRadialOffset[[sid]][strand]
        text_y <- center_pt[2] - normal_y * geneLabelRadialOffset[[sid]][strand]
        
        # 6. 计算文本角度及对齐方式（不变）
        base_angle <- atan2(dy, dx) * 180 / pi
        text_angle <- base_angle + 90 + geneLabelRotation[[sid]][strand]
        
        # 自动对齐逻辑（根据方向调整左/右对齐）
        if (strand == "+" && orient == 1) {
          hjust <- 1
        } else if (strand == "+" && orient != 1) {
          hjust <- 0
        } else if (strand == "-" && orient == 1) {
          hjust <- 0
        } else {
          hjust <- 1
        }
        
        # 确保文本方向正确（不颠倒）
        text_angle <- (text_angle + 360) %% 360
        if (text_angle > 90 && text_angle < 270) {
          text_angle <- text_angle + 180
          hjust <- 1 - hjust
        }
        text_angle <- text_angle %% 360
        
        data.frame(
          text = gene$anno,
          text_x = text_x,
          text_y = text_y,
          text_angle = text_angle,
          hjust = hjust,
          vjust = 0.5,
          seq_id = sid,
          group = i,
          stringsAsFactors = FALSE
        )
      }))
    }
  }
  
  
  # 定义旋转函数
  rotate_df <- function(df) {
    if (all(c("x","y") %in% names(df))) {
      x0 <- df$x; y0 <- df$y
      df$x <- x0 * cos(rot_rad) - y0 * sin(rot_rad)
      df$y <- x0 * sin(rot_rad) + y0 * cos(rot_rad)
    }
    if (all(c("x0","y0","x1","y1") %in% names(df))) {
      X <- df$x0; Y <- df$y0
      df$x0 <- X * cos(rot_rad) - Y * sin(rot_rad)
      df$y0 <- X * sin(rot_rad) + Y * cos(rot_rad)
      X1 <- df$x1; Y1 <- df$y1
      df$x1 <- X1 * cos(rot_rad) - Y1 * sin(rot_rad)
      df$y1 <- X1 * sin(rot_rad) + Y1 * cos(rot_rad)
    }
    if (all(c("label_x","label_y") %in% names(df))) {
      LX <- df$label_x; LY <- df$label_y
      df$label_x <- LX * cos(rot_rad) - LY * sin(rot_rad)
      df$label_y <- LX * sin(rot_rad) + LY * cos(rot_rad)
    }
    if (all(c("x_start","y_start","x_end","y_end") %in% names(df))) {
      Xs <- df$x_start; Ys <- df$y_start
      df$x_start <- Xs * cos(rot_rad) - Ys * sin(rot_rad)
      df$y_start <- Xs * sin(rot_rad) + Ys * cos(rot_rad)
      Xe <- df$x_end; Ye <- df$y_end
      df$x_end <- Xe * cos(rot_rad) - Ye * sin(rot_rad)
      df$y_end <- Xe * sin(rot_rad) + Ye * cos(rot_rad)
    }
    if (all(c("text_x","text_y") %in% names(df))) {
      TX <- df$text_x; TY <- df$text_y
      df$text_x <- TX * cos(rot_rad) - TY * sin(rot_rad)
      df$text_y <- TX * sin(rot_rad) + TY * cos(rot_rad)
      df$text_angle <- df$text_angle + rotation
    }
    df
  }
  
  # 应用旋转
  seqArcs <- lapply(seqArcs, rotate_df)
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
  
  extremes <- get_plot_extremes(
    allRibbon = allRibbon,
    seqArcs = seqArcs,
    axisLines = axisLines,
    axisTicks = axisTicks,
    gene_polys = gene_polys,
    gene_arrows = NULL,
    show_axis = show_axis
  )
  
  # 13. 绘制图形（使用ggnewscale实现双fill映射）
  chordPlot <- chordPlotFunc(allRibbon,ribbon_alpha,ribbon_color_scheme,ribbon_colors,show_legend,gene_polys,gene_pal,gene_color_scheme,seqArcs,gene_arrows,show_axis,axisLines,axisTicks,axisLabelOrientation,seq_colors,seq_labels,extremes,panel_margin,title)
  
  return(chordPlot)
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
  #gene_track = gene_track,
  #title = "Multi-sequence Chord Diagram with Gene Annotations"
  
)

ggview::ggview(p_final, width = 12, height = 12)


