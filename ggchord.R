#' ggchord: 多序列弦图（带基因注释箭头）
#'
#' 绘制多序列间的BLAST比对弦图，并在序列上用箭头显示基因注释，箭头方向由strand决定
#'
#' @param seq_data data.frame/tibble，包含序列信息，必须包含列：
#'   - seq_id: 序列唯一标识
#'   - length: 序列长度
#' @param ribbon_data data.frame/tibble，包含BLAST比对结果，必须包含列：
#'   - qaccver: 查询序列ID
#'   - saccver: 目标序列ID
#'   - length: 比对长度
#'   - pident: 序列相似度百分比
#'   - qstart: 查询序列起始位置
#'   - qend: 查询序列结束位置
#'   - sstart: 目标序列起始位置
#'   - send: 目标序列结束位置
#' @param gene_track data.frame/tibble，包含基因注释信息，必须包含列：
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
#' @param ribbon_color_scheme 字符，连接带配色方案，可选 "single"、"query" 或 "pident"，默认"single"
#' @param ribbon_colors 连接带颜色参数：
#'   - single: 单一颜色（单值或向量取第一个）
#'   - query: 按查询序列映射颜色（命名或非命名向量或单值）
#'   - pident: 渐变色阶向量，用于按相似度百分比生成渐变
#' @param ribbon_alpha 数值，连接带透明度 [0,1]，默认0.6
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
#' @param debug 逻辑值，是否打印调试信息，默认FALSE
#' @return 返回一个 ggplot2 对象
#' @import ggplot2
#' @import RColorBrewer
#' @import grDevices
#' @export

# 加载所需包
library(ggplot2)
library(RColorBrewer)
library(grDevices)

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
  # 初始化结果列表，默认值为default的正负链
  result <- setNames(lapply(seqs, function(id) c("+" = default, "-" = default)), seqs)
  
  if (is.null(gene_offset)) {
    return(result)
  }
  
  # 处理单值情况
  if (length(gene_offset) == 1 && !is.list(gene_offset) && is.numeric(gene_offset)) {
    val <- gene_offset
    return(setNames(lapply(seqs, function(id) c("+" = val, "-" = val)), seqs))
  }
  
  # 处理向量情况（长度与序列数相同）
  if (is.vector(gene_offset) && !is.list(gene_offset) && length(gene_offset) == n) {
    if (!is.null(names(gene_offset))) {
      # 命名向量，按名称匹配
      for (id in seqs) {
        if (id %in% names(gene_offset)) {
          val <- gene_offset[id]
          result[[id]] <- c("+" = val, "-" = val)
        }
      }
    } else {
      # 非命名向量，按顺序匹配
      for (i in seq_along(seqs)) {
        val <- gene_offset[i]
        result[[seqs[i]]] <- c("+" = val, "-" = val)
      }
    }
    return(result)
  }
  
  # 处理列表情况
  if (is.list(gene_offset)) {
    # 检查列表中的名称
    list_names <- names(gene_offset)
    if (is.null(list_names)) {
      # 非命名列表，按顺序匹配序列
      if (length(gene_offset) != n) {
        stop("gene_offset 列表长度与序列数不匹配，且未命名")
      }
      for (i in seq_along(seqs)) {
        elem <- gene_offset[[i]]
        # 处理每个元素
        if (length(elem) == 1 && is.numeric(elem)) {
          # 单值，正负链相同
          result[[seqs[i]]] <- c("+" = elem, "-" = elem)
        } else if (is.vector(elem) && all(names(elem) %in% c("+", "-"))) {
          # 命名向量，包含+/-
          current <- result[[seqs[i]]] # 初始为默认
          if ("+" %in% names(elem)) current["+"] <- elem["+"]
          if ("-" %in% names(elem)) current["-"] <- elem["-"]
          result[[seqs[i]]] <- current
        } else {
          stop(paste("gene_offset 列表中第", i, "个元素格式错误，应为单值或包含'+', '-'的命名向量"))
        }
      }
    } else {
      # 命名列表，按名称匹配
      for (id in list_names) {
        if (!id %in% seqs) {
          stop(paste("gene_offset 列表包含未知序列ID:", id))
        }
        elem <- gene_offset[[id]]
        current <- result[[id]] # 初始为默认
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
  
  # 无法处理的格式
  stop("gene_offset 格式错误，请提供单值、向量或列表")
}

# 主函数：绘制多序列弦图
ggchord <- function(
    seq_data,
    ribbon_data,
    gene_track,
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
    ribbon_color_scheme    = c("single","query","pident"),
    ribbon_colors          = NULL,
    ribbon_alpha           = 0.6,
    ribbon_ctrl_point      = NULL,
    ribbon_gap             = 0.1,
    axis_gap               = 0.05,
    axis_tick_major_number        = 5,
    axis_tick_major_length = 0.02,
    axis_tick_minor_number        = 4,
    axis_tick_minor_length = 0.01,
    axis_label_size        = 3,
    axis_label_offset      = 0.1,
    rotation               = 45,
    show_legend            = TRUE,
    debug                  = FALSE
) {
  # 检查必要的包是否安装
  if (!"ggplot2" %in% installed.packages()) stop("需要安装 ggplot2 包")
  if (!"RColorBrewer" %in% installed.packages()) stop("需要安装 RColorBrewer 包")
  if (!"grDevices" %in% installed.packages()) stop("需要安装 grDevices 包")
  
  ribbon_color_scheme <- match.arg(ribbon_color_scheme)
  
  # 1. 验证输入数据格式
  # 验证序列数据
  required_seq_cols <- c("seq_id", "length")
  if (!all(required_seq_cols %in% colnames(seq_data))) {
    stop("seq_data 必须包含以下列: ", paste(required_seq_cols, collapse = ", "))
  }
  if (any(seq_data$length <= 0)) {
    stop("seq_data 中的 length 必须为正数")
  }
  
  # 验证BLAST数据
  required_ribbon_cols <- c("qaccver", "saccver", "length", "pident", 
                            "qstart", "qend", "sstart", "send")
  if (!all(required_ribbon_cols %in% colnames(ribbon_data))) {
    stop("ribbon_data 必须包含以下列: ", paste(required_ribbon_cols, collapse = ", "))
  }
  fb_all <- ribbon_data  # 直接使用传入的已处理比对数据
  if (nrow(fb_all) == 0) stop("ribbon_data 中没有有效比对数据")
  if (debug) cat("使用的比对数据行数:", nrow(fb_all), "\n")
  
  # 验证基因注释数据
  required_gene_cols <- c("seq_id", "start", "end", "strand", "anno")
  if (!all(required_gene_cols %in% colnames(gene_track))) {
    stop("gene_track 必须包含以下列: ", paste(required_gene_cols, collapse = ", "))
  }
  if (nrow(gene_track) == 0) warning("gene_track 中没有有效基因注释数据")
  if (debug) cat("使用的基因注释数据行数:", nrow(gene_track), "\n")
  
  # 2. 处理序列信息
  seqs <- seq_data$seq_id
  lens <- setNames(seq_data$length, seqs)
  if (length(seqs) < 2) stop("seq_data 中至少需要包含2条序列")
  
  # 处理序列顺序
  if (!is.null(seq_order)) {
    if (!all(seq_order %in% seqs)) {
      stop("seq_order 包含未知序列ID: ", paste(setdiff(seq_order, seqs), collapse = ", "))
    }
    seqs <- seq_order
    lens <- lens[seqs]  # 按新顺序重新排列长度
  }
  n <- length(seqs)  # 序列数量
  
  # 3. 处理序列相关参数
  seq_labels      <- process_sequence_param(seq_labels, seqs, "seq_labels", default_value = seqs)
  seqRadius       <- process_sequence_param(seq_radius, seqs, "seq_radius", default_value = 1.0)
  ribbonGap       <- process_sequence_param(ribbon_gap, seqs, "ribbon_gap", default_value = 0.1)
  axisGap         <- process_sequence_param(axis_gap, seqs, "axis_gap", default_value = 0.05)
  # 处理基因偏移参数（新逻辑）
  geneGap         <- process_gene_offset(gene_offset, seqs, default = 0.03)
  geneWidth       <- process_sequence_param(gene_width, seqs, "gene_width", default_value = 0.1)
  axisMaj         <- process_sequence_param(axis_tick_major_number, seqs, "axis_tick_major_number", default_value = 5)
  axisMajLen      <- process_sequence_param(axis_tick_major_length, seqs, "axis_tick_major_length", default_value = 0.02)
  axisMin         <- process_sequence_param(axis_tick_minor_number, seqs, "axis_tick_minor_number", default_value = 4)
  axisMinLen      <- process_sequence_param(axis_tick_minor_length, seqs, "axis_tick_minor_length", default_value = 0.01)
  labelSize       <- process_sequence_param(axis_label_size, seqs, "axis_label_size", default_value = 3)
  labelOffset     <- process_sequence_param(axis_label_offset, seqs, "axis_label_offset", default_value = 0)
  orientation     <- process_sequence_param(seq_orientation, seqs, "seq_orientation", default_value = 1)
  rot_rad <- rotation * pi / 180  # 转换为弧度
  
  # 处理序列间隙参数
  seq_gap <- process_sequence_param(seq_gap, seqs, "seq_gap")
  if (any(seq_gap < 0 | seq_gap >= 0.5)) {
    stop("seq_gap 必须在 [0, 0.5) 范围内")
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
  
  # 4. 处理连接带颜色参数
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
  
  # 5. 计算序列弧度和间隙弧度
  total_circ <- 2 * pi  # 总圆周弧度
  total_gap_prop <- sum(seq_gap)  # 总间隙比例
  
  if (total_gap_prop >= 1) {
    stop("seq_gap 总和不能超过1（无法容纳所有序列）")
  }
  
  seq_total_prop <- 1 - total_gap_prop  # 序列总占比
  sum_lens <- sum(lens)
  theta <- (lens / sum_lens) * total_circ * seq_total_prop  # 每个序列的弧度
  gap_rads <- total_circ * seq_gap  # 每个间隙的弧度（实际角度）
  
  # 6. 计算序列起始和结束角度
  starts <- numeric(n)
  starts[1] <- 0  # 第一个序列从0开始
  
  if (n > 1) {
    for (i in 2:n) {
      starts[i] <- starts[i-1] + theta[i-1] + gap_rads[i-1]
    }
  }
  
  ends <- starts + theta
  names(starts) <- names(ends) <- seqs
  
  # 7. 准备序列外层和内层坐标
  nSeg <- 500  # 每个序列的分段数（控制平滑度）
  seqArcs <- lapply(seqs, function(id) {
    angs <- seq(starts[id], ends[id], length.out = nSeg)
    if (orientation[id] == -1) angs <- rev(angs)  # 反向序列
    data.frame(x = seqRadius[id] * cos(angs), 
               y = seqRadius[id] * sin(angs), 
               seq_id = id)
  })
  
  innerArcs <- lapply(seqs, function(id) {
    angs <- seq(starts[id], ends[id], length.out = nSeg)
    if (orientation[id] == -1) angs <- rev(angs)
    r <- seqRadius[id] - ribbonGap[id]  # 内层半径（连接带起始位置）
    data.frame(x = r * cos(angs), y = r * sin(angs))
  })
  names(innerArcs) <- seqs
  
  # 8. 生成坐标轴刻度和线
  axisTicks <- do.call(rbind, lapply(seqs, function(id) {
    majors <- breakPointsFunc(lens[id], axisMaj[id])
    minors <- unlist(lapply(seq_len(length(majors)-1), function(i) {
      seq(majors[i], majors[i+1], length.out = axisMin[id] + 2)[-c(1, axisMin[id] + 2)]
    }))
    pts <- c(majors, minors)
    do.call(rbind, lapply(pts, function(v) {
      frac <- v / lens[id]  # 相对位置比例
      # 计算角度（考虑方向）
      ang <- if (orientation[id] == 1) {
        starts[id] + frac * (ends[id] - starts[id])
      } else {
        ends[id] - frac * (ends[id] - starts[id])
      }
      # 基准半径（序列半径 + 坐标轴间隙）
      base_r <- seqRadius[id] + axisGap[id]
      # 刻度长度（根据主/次刻度区分）
      tick_len <- if (v %in% majors) axisMajLen[id] else axisMinLen[id]
      # 方向（axis_gap >=0 向外，否则向内）
      dir <- if (axisGap[id] >= 0) 1 else -1
      end_r <- base_r + tick_len * dir
      # 标签位置
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
  
  # 坐标轴连接线
  axisLines <- do.call(rbind, lapply(seqs, function(id) {
    angs <- seq(starts[id], ends[id], length.out = nSeg)
    if (orientation[id] == -1) angs <- rev(angs)
    r <- seqRadius[id] + axisGap[id]
    data.frame(x = r * cos(angs), y = r * sin(angs), seq_id = id)
  }))
  
  # 9. 生成连接带数据
  ribbons <- list()
  cntValid <- cntInvalid <- 0
  for (i in seq_len(nrow(fb_all))) {
    row <- fb_all[i,]
    q <- row$qaccver
    s <- row$saccver
    if (q == s) {  # 跳过自身比对
      cntInvalid <- cntInvalid + 1
      next
    }
    if (!q %in% seqs || !s %in% seqs) {  # 跳过不在序列表中的比对
      cntInvalid <- cntInvalid + 1
      next
    }
    
    # 获取序列的内层弧线数据
    qr <- innerArcs[[q]]
    sr <- innerArcs[[s]]
    
    # 计算比对位置在弧线上的索引
    qi0 <- map_idx((row$qstart - 1)/lens[q], seq(starts[q], ends[q], length.out = nSeg))
    qi1 <- map_idx((row$qend - 1)/lens[q], seq(starts[q], ends[q], length.out = nSeg))
    si0 <- map_idx((row$sstart - 1)/lens[s], seq(starts[s], ends[s], length.out = nSeg))
    si1 <- map_idx((row$send - 1)/lens[s], seq(starts[s], ends[s], length.out = nSeg))
    
    # 过滤过短的比对
    if (abs(qi1 - qi0) < 1 || abs(si1 - si0) < 1) {
      cntInvalid <- cntInvalid + 1
      next
    }
    
    # 提取比对片段
    segQ <- qr[qi0:qi1, , drop = FALSE]
    segS <- sr[si0:si1, , drop = FALSE]
    
    # 确定贝塞尔曲线控制点
    if (!is.null(ribbon_ctrl_point)) {
      if (is.list(ribbon_ctrl_point)) {
        cp <- ribbon_ctrl_point[[i]]
        c1 <- cp[[1]]
        c2 <- cp[[2]]
      } else {
        c1 <- ribbon_ctrl_point
        c2 <- ribbon_ctrl_point
      }
    } else {
      c1 <- c2 <- c(0, 0)  # 默认控制点为圆心
    }
    
    # 生成贝塞尔曲线
    b1 <- bezier_pts(as.numeric(segQ[1, ]), as.numeric(segS[1, ]), c1, c1)  # 起点曲线（查询起点→目标起点）
    b2 <- bezier_pts(as.numeric(segQ[nrow(segQ), ]), as.numeric(segS[nrow(segS), ]), c2, c2)  # 终点曲线（查询终点→目标终点）
    
    # 构建闭合多边形
    poly <- rbind(
      segQ,                   # 1. 查询序列片段：查询起点 → 查询终点（正向）
      b2,                     # 2. 终点贝塞尔曲线：查询终点 → 目标终点（正向）
      segS[nrow(segS):1, ],   # 3. 目标序列片段：目标终点 → 目标起点（反向，与b2终点衔接）
      b1[nrow(b1):1, ]        # 4. 起点贝塞尔曲线：目标起点 → 查询起点（反向，闭合回到起点）
    )
    
    # 设置填充色
    if (ribbon_color_scheme == "pident") {
      poly$pident <- row$pident
    } else {
      poly$fill <- switch(ribbon_color_scheme,
                          single = singleCol,
                          query  = queryCols[q])
    }
    poly$group <- i  # 分组标识
    ribbons[[length(ribbons) + 1]] <- poly
    cntValid <- cntValid + 1
  }
  
  if (debug) {
    cat("有效连接带:", cntValid, "无效连接带:", cntInvalid, "\n")
  }
  if (cntValid == 0) {
    stop("没有有效连接带可绘制，请检查输入数据或参数")
  }
  allRibbon <- do.call(rbind, ribbons)
  
  # 10. 处理基因注释箭头（曲边三角头箭头）
  gene_polys <- data.frame()
  if (nrow(gene_track) > 0) {
    valid_genes <- gene_track[gene_track$seq_id %in% seqs, ]
    for (i in seq_len(nrow(valid_genes))) {
      gene   <- valid_genes[i, ]
      sid    <- gene$seq_id
      strand <- gene$strand
      
      # 1. 计算 ang_s, ang_e（已兼顾 seq_orientation 和 strand）
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
      
      # 2. 半径 & 参数（根据链方向获取偏移值）
      r0        <- seqRadius[sid] + geneGap[[sid]][strand]
      width     <- geneWidth[sid]
      r_out     <- r0 + width/2
      r_in      <- r0 - width/2
      head_frac <- 0.3
      a_break   <- a_s + (1-head_frac)*(a_e - a_s)
      
      # 3. 构造 body 部分
      n_body <- 40
      ang_body <- seq(a_s, a_break, length.out = n_body)
      body_out <- data.frame(
        x = r_out * cos(ang_body),
        y = r_out * sin(ang_body),
        group = i, ord = seq_len(n_body)
      )
      body_in <- data.frame(
        x = r_in * cos(rev(ang_body)),
        y = r_in * sin(rev(ang_body)),
        group = i, ord = seq_len(n_body) + 3*n_body
      )
      
      # 4. 构造 head 边：参数 t=[0,1]
      n_head <- 40
      t_seq  <- seq(0, 1, length.out = n_head)
      # 外侧头边
      ang_out <- a_break + t_seq*(a_e - a_break)
      r_out_edge <- r_out + t_seq*(r0 - r_out)
      edge_out <- data.frame(
        x = r_out_edge * cos(ang_out),
        y = r_out_edge * sin(ang_out),
        group = i, ord = seq_len(n_head) + n_body
      )
      # 内侧头边（逆序角度，收敛到同一尖点）
      ang_in <- a_break + (1 - t_seq)*(a_e - a_break)
      r_in_edge <- r_in + (1 - t_seq)*(r0 - r_in)
      edge_in <- data.frame(
        x = r_in_edge * cos(ang_in),
        y = r_in_edge * sin(ang_in),
        group = i, ord = seq_len(n_head) + 2*n_body
      )
      
      gene_polys <- rbind(gene_polys, body_out, edge_out, edge_in, body_in)
    }
  }
  
  # 10. 处理基因注释箭头（旧版，保留但优先使用曲边箭头）
  gene_arrows <- data.frame()
  if (nrow(gene_track) > 0 && gene_label_show) {
    # 过滤不在当前序列中的基因注释
    valid_genes <- gene_track[gene_track$seq_id %in% seqs, ]
    
    if (nrow(valid_genes) > 0) {
      gene_arrows <- do.call(rbind, lapply(1:nrow(valid_genes), function(i) {
        gene <- valid_genes[i, ]
        seq_id <- gene$seq_id
        start <- gene$start
        end <- gene$end
        strand <- gene$strand
        anno <- gene$anno
        
        # 确保start < end，方便后续计算
        if (start > end) {
          temp <- start
          start <- end
          end <- temp
        }
        
        # 计算基因在序列上的相对位置
        start_frac <- start / lens[seq_id]
        end_frac <- end / lens[seq_id]
        
        # 计算基因起点和终点的角度（考虑序列方向）
        if (orientation[seq_id] == 1) {
          start_ang <- starts[seq_id] + start_frac * (ends[seq_id] - starts[seq_id])
          end_ang <- starts[seq_id] + end_frac * (ends[seq_id] - starts[seq_id])
        } else {
          start_ang <- ends[seq_id] - start_frac * (ends[seq_id] - starts[seq_id])
          end_ang <- ends[seq_id] - end_frac * (ends[seq_id] - starts[seq_id])
        }
        
        # 确定箭头的半径（在序列线外侧或内侧，根据链方向获取偏移）
        arrow_radius <- seqRadius[seq_id] + geneGap[[seq_id]][strand]
        
        # 计算箭头起点和终点坐标
        if (strand == "+") {
          # 与序列方向相同
          x_start <- arrow_radius * cos(start_ang)
          y_start <- arrow_radius * sin(start_ang)
          x_end <- arrow_radius * cos(end_ang)
          y_end <- arrow_radius * sin(end_ang)
        } else {
          # 与序列方向相反
          x_start <- arrow_radius * cos(end_ang)
          y_start <- arrow_radius * sin(end_ang)
          x_end <- arrow_radius * cos(start_ang)
          y_end <- arrow_radius * sin(start_ang)
        }
        
        # 计算注释文本位置（箭头中间稍微向外一点）
        mid_ang <- (start_ang + end_ang) / 2
        text_radius <- arrow_radius + 0.15  # 文本位置比箭头稍向外
        text_x <- text_radius * cos(mid_ang)
        text_y <- text_radius * sin(mid_ang)
        
        # 确定文本角度
        text_angle <- mid_ang * (180/pi) - 0  # 转换为度并调整方向
        if (text_angle > 90) text_angle <- text_angle - 180
        if (text_angle < -90) text_angle <- text_angle + 180
        
        data.frame(
          x_start = x_start,
          y_start = y_start,
          x_end = x_end,
          y_end = y_end,
          text = anno,
          text_x = text_x,
          text_y = text_y,
          text_angle = text_angle,
          seq_id = seq_id,
          group = i
        )
      }))
    } else {
      warning("gene_track 中没有与 seq_data 匹配的序列ID")
    }
  }
  
  # 定义旋转函数
  rotate_df <- function(df) {
    # 对 (x,y)
    if (all(c("x","y") %in% names(df))) {
      x0 <- df$x; y0 <- df$y
      df$x <-  x0 * cos(rot_rad) - y0 * sin(rot_rad)
      df$y <-  x0 * sin(rot_rad) + y0 * cos(rot_rad)
    }
    # 对 (x0,y0) 和 (x1,y1)
    if (all(c("x0","y0") %in% names(df))) {
      X <- df$x0; Y <- df$y0
      df$x0 <-  X * cos(rot_rad) - Y * sin(rot_rad)
      df$y0 <-  X * sin(rot_rad) + Y * cos(rot_rad)
      X1 <- df$x1; Y1 <- df$y1
      df$x1 <-  X1 * cos(rot_rad) - Y1 * sin(rot_rad)
      df$y1 <-  X1 * sin(rot_rad) + Y1 * cos(rot_rad)
    }
    # 对标签位置 (label_x,label_y)
    if (all(c("label_x","label_y") %in% names(df))) {
      LX <- df$label_x; LY <- df$label_y
      df$label_x <-  LX * cos(rot_rad) - LY * sin(rot_rad)
      df$label_y <-  LX * sin(rot_rad) + LY * cos(rot_rad)
    }
    # 对基因箭头 (x_start,y_start) 和 (x_end,y_end)
    if (all(c("x_start","y_start","x_end","y_end") %in% names(df))) {
      Xs <- df$x_start; Ys <- df$y_start
      df$x_start <-  Xs * cos(rot_rad) - Ys * sin(rot_rad)
      df$y_start <-  Xs * sin(rot_rad) + Ys * cos(rot_rad)
      Xe <- df$x_end; Ye <- df$y_end
      df$x_end <-  Xe * cos(rot_rad) - Ye * sin(rot_rad)
      df$y_end <-  Xe * sin(rot_rad) + Ye * cos(rot_rad)
    }
    # 对基因文本位置 (text_x,text_y)
    if (all(c("text_x","text_y") %in% names(df))) {
      TX <- df$text_x; TY <- df$text_y
      df$text_x <-  TX * cos(rot_rad) - TY * sin(rot_rad)
      df$text_y <-  TX * sin(rot_rad) + TY * cos(rot_rad)
      # 调整文本角度
      df$text_angle <- df$text_angle + rotation
    }
    df
  }
  
  # 应用旋转
  seqArcs   <- lapply(seqArcs,   rotate_df)
  innerArcs <- lapply(innerArcs, rotate_df)
  axisLines <- rotate_df(axisLines)
  axisTicks <- rotate_df(axisTicks)
  allRibbon <- rotate_df(allRibbon)
  if (nrow(gene_arrows) > 0) {
    gene_arrows <- rotate_df(gene_arrows)
  }
  if (nrow(gene_polys) > 0) {
    gene_polys <- rotate_df(gene_polys)
    gene_polys <- gene_polys[with(gene_polys, order(group, ord)), ]
  }
  
  # 11. 绘制图形
  p <- ggplot() +
    # 绘制连接带
    geom_polygon(
      data = allRibbon,
      aes(
        x, y, group = group,
        fill = if (ribbon_color_scheme == "pident") pident else fill
      ),
      alpha = ribbon_alpha, 
      color = "grey", 
      linewidth = 0.2
    ) +
    # 绘制序列弧线
    geom_path(
      data = do.call(rbind, seqArcs),
      aes(x, y, group = seq_id, color = seq_id),
      linewidth = 1.2,
      arrow = arrow(angle = 25, length = unit(2, "mm"), type = "closed")
    ) +
    # 基因注释曲边箭头
    { if (nrow(gene_polys) > 0)
      geom_polygon(
        data  = gene_polys,
        aes(x = x, y = y, group = group),
        fill  = "blue",
        color = NA
      )
    } +
    # 绘制基因标签
    { if (nrow(gene_arrows) > 0 && gene_label_show)
      geom_text(
        data = gene_arrows,
        aes(x = text_x, y = text_y, label = text, angle = text_angle),
        size = gene_label_size,
        color = "black",
        inherit.aes = FALSE
      )
    } +
    # 绘制坐标轴线
    geom_path(
      data = axisLines,
      aes(x, y, group = seq_id),
      color = "black",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    # 绘制刻度线
    geom_segment(
      data = axisTicks,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      color = "black",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    # 绘制刻度标签
    geom_text(
      data = subset(axisTicks, !is.na(label)),
      aes(x = label_x, y = label_y, label = label, size = size),
      inherit.aes = FALSE,
      color = "black"
    ) +
    # 统一设置大小尺度
    scale_size_identity() +
    # 序列颜色
    scale_color_manual(
      name = "Seq ID",
      values = seq_colors,
      labels = seq_labels,
      guide = guide_legend(order = 1)
    ) +
    # 连接带颜色设置
    {
      if (ribbon_color_scheme == "pident") {
        scale_fill_stepsn(
          name = "Identity(%)",
          colours = ribbon_colors,
          limits = c(0, 100), 
          breaks = c(0,50,80,90,95,100),
          guide = guide_colorbar(theme = theme(legend.title.position = "top",
                                               legend.key.height = unit(90,"mm")),
                                 order = 3)
        )
      } else {
        scale_fill_identity(guide = "none")
      }
    } +
    
    # 主题设置
    theme_void() +
    coord_equal(clip = "off") +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.margin = margin(t = 0,r = 0,b = 0,l = 0),
      legend.box.spacing = unit(50,"mm"),
      legend.spacing = unit(20,"mm"),
      legend.position = if (show_legend) "right" else "none",
      legend.text = element_text(face = "bold", size = 10),
      legend.title = element_text(size = 14,face = "bold")
    )
  
  return(p)
}


# 示例使用方法
# 1. 读取数据
#读取序列长度数据
seq_data <- read.delim("seq_track.tsv", sep = "\t", stringsAsFactors = FALSE)

# 读取基因注释数据
gene_track <- read.delim("gene_track.tsv", sep = "\t", stringsAsFactors = FALSE)

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
  seq_gap = .03,
  seq_radius = c(4,3,2,1),
  seq_orientation = c(-1, 1, -1, -1),
  gene_offset = list(.3, 
                     .25, 
                     .4, 
                     c("+"=.3,"-"=-.1)),
  gene_width  = .08,
  gene_label_show = F,
  ribbon_gap = 0.2,
  ribbon_color_scheme = "pident",
  axis_gap = .1,
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

