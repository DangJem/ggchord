# API 变化与迁移记录

**状态：v0.16.0 开发中，首批修复已实现并通过局部直接验证；尚未正式发布。**

本文件只列实际变化与有明确标记的待决策事项。没有新增、删除或重命名公开函数，
也没有改变参数默认值；完整发布验收仍须用户明确授权后执行。

## 已实现变化

| 问题 | 旧行为 | 当前行为 | 用户影响与替代写法 | 删除版本 |
| --- | --- | --- | --- | --- |
| F01 | clean 的非 NULL 错误类型可选表被静默丢为 NULL | 明确报出 ribbon_data/gene_data 类型错误 | 传 data.frame/tibble；有意不提供时用 NULL；支持的区间对象先显式适配 | 不适用，未删除参数 |
| F02 | 反向 ribbon 排序报告的 original_value/new_value 相同 | 分别记录交换前后值 | 调用不变；读取 report 时得到真实变更；数据排序和 direction 语义不变 | 不适用 |
| F11 关联修复 | plasmid wrapper 的 shape 单参被误判为双形状参数 | 单参别名正常转发 | 既有 shape 调用恢复；显式两参仍报原冲突 | 不适用 |
| F31 | ribbon/gene 的缺失或空 ID 在 unknown_id 策略前报错 | 遵循已公开的 drop/error/keep，并保留原值/原行报告 | 默认 drop 会移除这些行；需立即停止时显式 error；keep 仅保留供审计，不保证可绘图 | 不适用，未删除参数 |
| F32 | feature 的 factor 起止标签被当成等级码，例如 100/300 变 1/2 | 按数值标签解释 factor 顶层起止列 | 既有调用恢复真实坐标；无需按等级重编码；非数值标签不是有效坐标 | 不适用 |

### F01：明确输入，而不是丢弃数据

以下示例的 seqs/ribbons 是完整教学输入；不代表实际科研数据。

```r
library(ggchord)
seqs <- data.frame(accver = c("q", "s"), length = c(20, 20))
ribbons <- data.frame(qaccver = "q", saccver = "s",
  qstart = 2, qend = 8, sstart = 15, send = 3)

# 旧错误调用曾静默得到 ribbon_data = NULL；当前明确报错。
# clean_ggchord_data(seqs, ribbon_data = "not a table")

# 显式传入数据；可选表有意缺席时才使用 NULL。
cleaned <- clean_ggchord_data(seqs, ribbon_data = ribbons, gene_data = NULL)
```

`GRanges`、`GInteractions` 等不是 clean 的新增输入类型。需要时先使用对应的
`as_ggchord_feature()` / `as_ggchord_link()` 转成表，再调用 clean；适配器自身的
坐标、链和可选依赖要求仍适用。不要为了绕过报错把本应保留的数据改成 NULL。

### F02：调用不变，修正报告

```r
sorted <- cleaned$report[cleaned$report$action == "sort", ]
stopifnot(nrow(sorted) == 1L,
  sorted$original_value == "qstart=2, qend=8, sstart=15, send=3",
  sorted$new_value == "qstart=2, qend=8, sstart=3, send=15",
  cleaned$ribbon_data$direction == "reverse",
  ribbons$sstart == 15, ribbons$send == 3)
```

修复限定于 ribbon 排序报告；gene 分支原本已能区分交换前后值。默认排序政策与
prepare 的方向保留差异没有在本轮调整，详见 F03 待决策事项。

### F11：恢复已有 shape 简写

```r
features <- data.frame(accver = "q", start = 2, end = 8,
  strand = "+", type = "CDS")
p <- ggchord(seqs[1, ]) + geom_seq() +
  scale_feature_fill_manual(values = c(CDS = "black")) +
  geom_feature_plasmid(data = features, shape = "block") + coord_circular()
built <- ggplot2::ggplot_build(p)
stopifnot(all(built$data[[2]]$feature_fill == "black"))
```

同样可以继续使用 `feature_shape = "block"`；本轮没有决定哪个名称长期保留。
`geom_feature_plasmid(feature_shape = "arrow", shape = "block")` 仍明确报错，
没有采用静默优先规则。默认几何、布局与用户 scale 优先行为保留。

### F31：缺失、空和未知 ID 遵循同一策略

用户于 2026-10-30 选择按既有文档处理，修复提交为
`2c90371ce2a27ee4083c07a396150e78c46040b0`。默认值仍是 `"drop"`，但缺失/空 ID
由旧前置报错改为按该策略移除；需要保留原先停止效果时显式使用 `"error"`。

```r
seq_ref <- data.frame(accver = c("q", "s"), length = c(20, 20))
bad_links <- data.frame(qaccver = NA_character_, saccver = "",
  qstart = 2, qend = 8, sstart = 15, send = 3)
bad_genes <- data.frame(accver = NA_character_, start = -5, end = -1,
  strand = "+", anno = NA_character_)

dropped_ids <- clean_ggchord_data(seq_ref, bad_links, bad_genes)
stopifnot(nrow(dropped_ids$ribbon_data) == 0L,
  nrow(dropped_ids$gene_data) == 0L,
  nrow(dropped_ids$report) == 3L,
  all(dropped_ids$report$row == 1L),
  sum(is.na(dropped_ids$report$original_value)) == 2L,
  sum(dropped_ids$report$original_value == "", na.rm = TRUE) == 1L)

kept_ids <- clean_ggchord_data(seq_ref, bad_links, bad_genes,
  unknown_id = "keep", empty_annotation = "replace")
stopifnot(identical(kept_ids$ribbon_data, bad_links),
  identical(kept_ids$gene_data, bad_genes),
  all(kept_ids$report$action == "keep"))

# 如需遇到缺失/空/未知 ID 就停止，而不是移行：
# clean_ggchord_data(seq_ref, bad_links, unknown_id = "error")
```

连接两端分别报告问题；一行同时存在两个错误 ID 时会有两条报告，源行号不变。
`original_value` 是字符型表示：NA 保留缺失，空串保留空串；原始 NaN 的文字表示为
`"NaN"`，可通过源行追溯原值。原始缺失掩码参与匹配，即使参考中有合法字面 ID
`"NaN"`，数值 NaN 也不会被当成该序列；keep 返回的原始 ID 值不被改写。

`error` 会停止并指出表、列、源行和问题值，不返回部分清洗结果。`keep` 对无参考
序列的行不清坐标/注释，也不保证通过后续验证。seq_data 的 ID 仍须非缺失、非空且唯一，
数值有限性和 gene strand 检查仍适用。prepare 的适配器仍在清洗步骤前拒绝缺失 ID；
因此需先独立 clean/drop 或显式修正 ID，再调用 prepare，而不是依赖其 steps$clean。

### F32：feature 因子列保留坐标标签

修复前，下面的标签可能被静默绘成 factor 等级 1/2；现在直接使用 100/300，
无需更改调用或修改原始数据。数值列及既有数值文本转换路径保持原义。

```r
factor_features <- data.frame(accver = "A",
  start = factor("100", levels = c("100", "300")),
  end = factor("300", levels = c("100", "300")), strand = "+", type = "CDS")
factor_plot <- ggchord(data.frame(accver = "A", length = 1000)) + geom_seq() +
  geom_feature(data = factor_features, feature_shape = "block") + coord_collinear()
factor_export <- export_ggchord_layout(factor_plot, include = "feature")
stopifnot(identical(unique(factor_export$feature$start), 100),
  identical(unique(factor_export$feature$end), 300),
  is.factor(factor_features$start), is.factor(factor_features$end))
```

这修正顶层起止列，不把任意类别标签猜成坐标，也不放宽已有有限数值要求。
复合 segment 子表的类型仍按其原有合同提供。

## 尚未实施的迁移

F03 默认差异、F07 getter/export、F09 gene/feature、F10 既有偏移弃用、F11 唯一
形状入口、标签及连接命名、颜色别名和重导出均尚未新增删除决定。后续若改变，
必须补入旧/新调用、替代方案、影响、错误提示、用户确认状态及确切删除版本；
不能将待确认的删除日期写成承诺。详细方案见 [决策台账](API-AUDIT-DECISIONS.md)。
