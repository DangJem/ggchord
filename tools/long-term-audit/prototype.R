# 仅供设计验证的隔离原型，不导出、不挂接包内构建流水线。
# A/B/C/D 与混合注释均为明确标注的合成测试数据，不冒充科研样本。
library(dplyr)
out <- "tools/long-term-audit"
write_csv <- function(x, name) utils::write.csv(x, file.path(out, name), row.names = FALSE)

register_instances <- function(seqs, seq_order) {
  stopifnot(!anyDuplicated(seqs$accver), all(seq_order %in% seqs$accver))
  ids <- names(seq_order)
  if (is.null(ids)) ids <- paste0("instance-", seq_along(seq_order))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop("invalid instance IDs")
  data.frame(instance_id = ids, accver = unname(seq_order),
    declared_order = seq_along(seq_order), stringsAsFactors = FALSE)
}

# 路由发生在 Coord 之前；不以 pident、行顺序或屏幕距离猜副本。
route_hsps <- function(hits, instances, between = NULL) {
  if (!is.null(between)) {
    if (length(between) != 2L || anyNA(between)) stop("between needs two instances")
    if (is.numeric(between)) {
      if (any(between != floor(between)) || any(!between %in% instances$declared_order)) stop("unknown between index")
      between <- instances$instance_id[match(between, instances$declared_order)]
    }
    if (any(!between %in% instances$instance_id) || anyDuplicated(between)) stop("invalid between")
  }
  routes <- vector("list", nrow(hits))
  for (i in seq_len(nrow(hits))) {
    row <- hits[i, ]
    q <- instances$instance_id[instances$accver == row$qaccver]
    s <- instances$instance_id[instances$accver == row$saccver]
    candidates <- expand.grid(q_instance = q, s_instance = s, stringsAsFactors = FALSE)
    if (!is.null(between)) candidates <- candidates[
      candidates$q_instance != candidates$s_instance &
      candidates$q_instance %in% between & candidates$s_instance %in% between, , drop = FALSE]
    if (nrow(candidates) != 1L) stop("HSP ", row$source_hsp_id, ": ambiguous or absent route; specify between")
    routes[[i]] <- cbind(row, candidates, displayed_relation_id = paste0("display-", row$source_hsp_id))
  }
  if (!length(routes)) return(data.frame())
  bind_rows(routes)
}

expect_error <- function(expr) {
  result <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  stopifnot(result)
}

seqs <- data.frame(accver = c("A", "B", "C", "D"), length = c(1000L, 500L, 200L, 800L))
hits <- data.frame(source_hsp_id = paste0("h", 1:4),
  qaccver = c("A", "B", "C", "D"), saccver = c("B", "C", "D", "A"),
  qstart = c(100L, 20L, 40L, 300L), qend = c(300L, 100L, 90L, 500L),
  sstart = c(30L, 180L, 200L, 700L), send = c(130L, 80L, 400L, 900L),
  pident = c(95, 85, 90, 97))
source_copy <- hits
aba_hits <- bind_rows(hits[1, ], mutate(hits[1, ], source_hsp_id = "h5", pident = 80))
instances <- register_instances(seqs, c(A_high = "A", B = "B", A_low = "A"))
routes <- bind_rows(
  route_hsps(aba_hits |> filter(pident >= 90), instances, c("A_high", "B")),
  route_hsps(aba_hits |> filter(pident < 90), instances, c("B", "A_low")))
stopifnot(nrow(routes) == nrow(aba_hits), !anyDuplicated(routes$source_hsp_id),
  routes$q_instance == c("A_high", "A_low"), all(routes$s_instance == "B"),
  identical(hits, source_copy))
expect_error(route_hsps(aba_hits, instances))
expect_error(register_instances(seqs, c(x = "A", x = "B")))
expect_error(route_hsps(aba_hits, instances, c(1.5, 2)))
numeric_routes <- bind_rows(route_hsps(aba_hits[1, ], instances, c(1, 2)),
  route_hsps(aba_hits[2, ], instances, c(2, 3)))
stopifnot(identical(routes, numeric_routes))
write_csv(instances, "prototype-instances.csv")
write_csv(routes, "prototype-routes.csv")

# 实体默认 + 实例覆盖，向量位置不构成优先级。
resolve_orientation <- function(values, instances, default = 1) {
  if (anyDuplicated(names(values)) || any(!names(values) %in% c(instances$instance_id, instances$accver))) stop("unknown or duplicate key")
  vapply(seq_len(nrow(instances)), function(i) {
    for (key in c(instances$instance_id[i], instances$accver[i]))
      if (key %in% names(values)) return(unname(values[[key]]))
    default
  }, numeric(1))
}
stopifnot(identical(resolve_orientation(c(A_low = -1, A = 1), instances), c(1, 1, -1)),
  identical(resolve_orientation(c(A = 1, A_low = -1), instances), c(1, 1, -1)))

# 辐条环形邻接只筛选真实输入；N=2 不产生两份关系。
radial_select <- function(hits, order) {
  if (length(order) < 2L || anyDuplicated(order)) stop("need distinct spokes")
  index <- match(hits$qaccver, order); target <- match(hits$saccver, order)
  delta <- abs(index - target)
  hits[!is.na(delta) & (delta == 1L | delta == length(order) - 1L), , drop = FALSE]
}
closed <- radial_select(hits, seqs$accver)
open <- radial_select(filter(hits, source_hsp_id != "h4"), seqs$accver)
stopifnot(nrow(closed) == 4L, nrow(open) == 3L,
  !"h4" %in% open$source_hsp_id,
  nrow(radial_select(hits[1, ], c("A", "B"))) == 1L)
write_csv(closed, "prototype-radial-with-closure.csv")
write_csv(open, "prototype-radial-without-closure.csv")

# 两轴 HSP 端点直接映射，负向段保留负斜率，不推断内部 mismatch。
dot <- hits |>
  transmute(source_hsp_id, qaccver, saccver, x = qstart, xend = qend,
    y = sstart, yend = send, slope = (send - sstart) / (qend - qstart))
stopifnot(dot$slope[1] > 0, dot$slope[2] < 0)
write_csv(dot, "prototype-dotplot.csv")

# equal 采用边界投影 [start-1,end]；当前包惯例需桥接，不改现有映射。
boundaries <- seqs |>
  mutate(source_length = length, display_length = 1,
    first_boundary = 0, last_boundary = 1, first_center = .5 / length,
    last_center = (length - .5) / length, factor = 1 / length)
stopifnot(all(boundaries$last_boundary == 1), identical(seqs$length, c(1000L, 500L, 200L, 800L)))
write_csv(boundaries, "prototype-equal-boundaries.csv")

# Feature 类型策略原型：显式拆分图层，不猜未知类型，不自动隐式排除其他层。
features <- data.frame(source_feature_id = paste0("f", 1:5), accver = "A",
  type = c("CDS", "tRNA", "rRNA", "promoter", "repeat"),
  start = seq(10, 90, by = 20), end = seq(19, 99, by = 20), strand = c("+", "+", "-", ".", "."))
genes <- features |> filter(type == "CDS")
other <- features |> filter(is.na(type) | type != "CDS")
stopifnot(!length(intersect(genes$source_feature_id, other$source_feature_id)),
  setequal(c(genes$source_feature_id, other$source_feature_id), features$source_feature_id))
write_csv(features, "prototype-features.csv")

# Count 源坐标分组与 ggplot2 结果比较；不将 count 偷换成分箱。
events <- data.frame(accver = "A", pos = c(10, 10, 20), weight = c(2, 3, 1))
counted <- events |> group_by(accver, pos) |> summarise(count = sum(weight), .groups = "drop")
native <- ggplot2::ggplot_build(ggplot2::ggplot(events, ggplot2::aes(x = pos, weight = weight)) + ggplot2::geom_bar())$data[[1]]
stopifnot(identical(counted$count, native$count), identical(counted$pos, as.numeric(native$x)))
write_csv(counted, "prototype-count.csv")

# 仅一个 PNG 概念图，不是新 Coord 已实现，也不是发布设备验收。
grDevices::png(file.path(out, "prototype-geometry.png"), width = 1200, height = 900, res = 120)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(c(0, 1000), c(0, 1000), type = "n", xlab = "query source position", ylab = "subject source position", main = "Synthetic HSP endpoints (not a public Coord)")
segments(dot$x, dot$y, dot$xend, dot$yend, col = c("#3B7A8C", "#B64D4D", "#3B7A8C", "#3B7A8C"), lwd = 3)
draw_radial <- function(data, title) {
  theta <- pi / 2 - 2 * pi * (seq_len(nrow(seqs)) - 1) / nrow(seqs)
  point <- function(ids, pos) {
    i <- match(ids, seqs$accver)
    r <- .25 + pos / max(seqs$length)
    cbind(r * cos(theta[i]), r * sin(theta[i]))
  }
  plot(c(-1.4, 1.4), c(-1.4, 1.4), type = "n", asp = 1, xlab = "layout x", ylab = "layout y", main = title)
  a <- point(seqs$accver, 0); b <- point(seqs$accver, seqs$length)
  segments(a[, 1], a[, 2], b[, 1], b[, 2], lwd = 3)
  text(b[, 1] * 1.12, b[, 2] * 1.12, labels = seqs$accver)
  q <- point(data$qaccver, (data$qstart + data$qend) / 2)
  s <- point(data$saccver, (data$sstart + data$send) / 2)
  segments(q[, 1], q[, 2], s[, 1], s[, 2], col = "#3B7A8C", lwd = 2)
}
draw_radial(closed, "Spokes: supplied edges AB/BC/CD/DA")
draw_radial(open, "Spokes: no DA input, no closing edge")
plot(c(0, 1000), c(.5, 3.5), type = "n", xlab = "source position (proportional)", ylab = "display instance", yaxt = "n", main = "A-B-A: route before projection")
axis(2, 1:3, rev(instances$instance_id), las = 1)
segments(0, 3:1, c(1000, 500, 1000), 3:1, lwd = 3)
for (i in seq_len(nrow(routes))) {
  q <- 4 - match(routes$q_instance[i], instances$instance_id)
  s <- 4 - match(routes$s_instance[i], instances$instance_id)
  segments(mean(c(routes$qstart[i], routes$qend[i])), q,
    mean(c(routes$sstart[i], routes$send[i])), s, lwd = 2, col = c("#3B7A8C", "#B64D4D")[i])
}
grDevices::dev.off()
jsonlite::write_json(list(status = "isolated design prototype only", input = "synthetic",
  assertions = c("unique source routes", "ambiguous routes rejected", "named/index routes agree",
  "instance override wins", "no phantom DA edge", "signed dotplot slope", "equal boundary projection",
  "disjoint feature partition", "weighted source count equals ggplot2"),
  public_api_implemented = FALSE), file.path(out, "prototype-results.json"), auto_unbox = TRUE, pretty = TRUE)
cat("Isolated model and geometry assertions passed\n")
