# API 变化与迁移记录

## v0.16.3 本地开发变化

本版保留 getter/export/preview/diagnose 四个入口、原签名和全部已描述字段，
无公开入口或参数删除版本。推荐消费表格用 export，完整布局用 getter；不是弃用。

| 问题 | 原行为 | 当前行为与迁移 |
| --- | --- | --- |
| F48 构建后空导出 | `ggplot_build(p); export_ggchord_layout(p)` 可返回空组件，复用 layer 还可能共 ID | 原调用修复，无需删除 build 或先 export；每次出现的 layer 均独立注册，原图层不改写 |
| F26 缓存 | 普通标签可能沿用另一设备的估计 | 所有布局检查设备/context；`get_chord_layout(p, build=FALSE)` 在上下文变更后明确失败，需要 `get_chord_layout(p)` 重建 |
| F49 预览 | `view_ggchord(p)` 可覆盖 p 的布局缓存 | 原调用返回相同类型的临时图片路径，改用渲染副本；正式保存仍用 ggsave |

```r
library(ggchord)
s <- data.frame(accver="A", length=1000)
p <- ggchord(s) + geom_seq()
invisible(ggplot2::ggplot_build(p))
# 原先此顺序可能为空；现在直接保留序列几何。
out <- export_ggchord_layout(p, include="seq")
stopifnot(nrow(out$seq) > 0)
```

字段、单位、来源和有限验证见 [引擎契约](API-ENGINE-CONTRACTS-v0.16.3.md)；
没有将 registry 字段降为内部后删除。

## v0.16.2 本地开发变化

本版不删除、重命名或新增公开函数参数，删除版本均不适用。保留 file/files 双入口、
clean 的 sort 和 prepare 的 keep 默认；完整职责见
[数据契约](API-DATA-CONTRACTS-v0.16.2.md)。以下调用均为明确标注的教学输入。

| 问题 | 原行为 | 当前行为与迁移 |
| --- | --- | --- |
| F42 互惠重叠 | 用短分母会删除被包含的小区块；降序端点漏检 | 调用不变，交集必须覆盖每个区间的阈值；降序有效，opposite 关系保留。若研究意图确是删除任意包含关系，请自行显式筛选，不再借 reciprocal overlap 实现 |
| F43 merge 来源与关系 | 未合并行的 alignment length 被包络覆盖；自定义分组可能跨 pair，顺序交叉可误合并 | 调用不变；原行统计保留，pair 始终分组，默认保留交叉区块；合并结果 `.length_estimated=TRUE`，report 新增 source_length/source_pident 列表。不把合并 length 当新比对统计 |
| F44 线性甲基化上下文 | linear 也把末端接到开头，可能虚构阻断上下文 | `circular=FALSE` 原调用恢复线性含义；确为环形才显式 TRUE，不能为了恢复旧误报而改 DNA |
| F45 引物标签 factor | 标签位置和 primer_start/end 显示等级码 | 调用不变，按 factor 的真实数值标签；原表不变 |
| F46 读取与来源 | 混用 file/files 时 file glob 不展开；数字式 ID 可被自动转换；缺物理源行 | 两参数均统一解析，规范化路径去重且保持首次顺序；ID 原样，新增 `.source_line`。有意重复文件时请读取一次后显式 rbind 并标注复制来源 |
| F46 输入失败 | 畸形行可填充/错位，编码未明确，错误阈值可能 recycling | 不齐字段、无效 UTF-8/路径和非法单值阈值给明确错误；需先修复输入或按已知源编码转成 UTF-8。不是通过清洗猜值 |
| F47 单bp方向 | 零端点顺序会否定已有reverse，或被聚成假同向包络；exact可丢显式不同方向 | 原调用恢复方向保留；merge/bundle保留单bp端点独立，dedup/density使用已有方向证据，否则unknown；不要求伪造非零端点 |

```r
library(ggchord)
seqs <- data.frame(accver = c("A", "B"), length = 1000)
blocks <- data.frame(qaccver = "A", saccver = "B", qstart = c(1, 101),
  qend = c(1000, 200), sstart = c(1, 101), send = c(1000, 200),
  length = c(1100, 100), pident = 95)
# 旧调用不改；过去小区块可能被移除，现在两条均保留。
distinct <- deduplicate_ggchord_ribbons(blocks, by = "overlap")
stopifnot(nrow(distinct$data) == 2L)
unchanged <- merge_ggchord_ribbons(blocks[1, ])
stopifnot(unchanged$data$length == 1100, !unchanged$data$.length_estimated)
adjacent <- transform(blocks, qend = c(100, 200), send = c(100, 200),
  length = c(120, 130))
merged <- merge_ggchord_ribbons(adjacent)
stopifnot(merged$data$.length_estimated,
  identical(merged$report$source_length[[1]], c(120, 130)))

f <- tempfile(fileext = ".blast")
writeLines(c("# source comment",
  "001\tNA\t95\t10\t0\t0\t1\t10\t10\t1\t0\t1"), f)
# file/files 指向同一文件时只读一次；无需迁移到单一名字。
imported <- read_blast(file = f, files = f, source_file = TRUE)
stopifnot(nrow(imported) == 1L, imported$qaccver == "001",
  imported$saccver == "NA", imported$.source_line == 2L)
unlink(f)

patterns <- data.frame(enzyme = "custom", motif = "ACGT", blocked_by = "Dcm",
  methylation_motif = "CCWGG")
linear <- find_restriction_sites(c(A = "AGGACGTCC"), database = patterns,
  circular = FALSE, methylation = "dam_dcm")
stopifnot(all(linear$methylation_status == "unblocked"))

primers <- data.frame(accver = "A", start = factor("100", levels = c("100", "300")),
  end = factor("300", levels = c("100", "300")), strand = "+", name = "primer")
labels <- export_ggchord_layout(ggchord(seqs[1, ]) + geom_seq() +
  geom_primer_label_repel(data = primers) + coord_circular(), include = "labels")$labels
stopifnot(all(stats::na.omit(labels$primer_start) == 100),
  all(stats::na.omit(labels$primer_end) == 300))
```

`.source_line` 是物理文件行，`as_*()` 的 `.source_row` 是本次对象输入行；
GenBank 使用 `.source_record`/`.source_feature`。不要把这几类来源当成同一个索引。
如果旧代码按固定列数或数字列序号取输入，应改为按字段名取列。
BLAST `na.strings`/`colClasses` 仍可通过 dots 显式提供；调用者自定义转换会影响保真。

非 UTF-8 文件须明确其真实编码后转换；例如已知 Latin-1 时可用
`writeLines(enc2utf8(readLines(old_file, encoding = "latin1")), utf8_file, useBytes = TRUE)`，
再调用 reader。不猜测未知编码。解析失败先检查路径、源行或选定解释器的错误输出。

bundle 绝对端点分组与标准遗传码候选未构成新增算法承诺。bundle 保留相对方向聚合，
需保留绝对来源方向时使用实际存在的 qstrand/sstrand 等 group_by 列。common-feature
protein 模式明确只作标准码六框架搜索；source geneticCode 元数据不选择翻译表，
非标准码用 `mode="dna"` 或外部 code-aware 分析，不能把无命中解释为生物学不存在。

## v0.16.1 本地开发变化

以下为本版实际迁移；所有删除均发生在 **v0.16.1**，旧调用明确报错，不再执行兼容
分支。用户已明确取消此前“等到 v1.0.0 才删除”的期限。保留 colour/color 拼写别名。

| 旧调用中的参数/取值 | 当前写法 | 影响 |
| --- | --- | --- |
| `geom_feature(feature_shape = "block")`（含 plasmid） | `shape = "block"` | 固定值只有一个名称；`aes(feature_shape = type)` 仍为映射 |
| `gene_width`、`feature_width` | `width` | gene/feature/预设共用宽度名称和既有单位 |
| gene 标签 `gene_label_*` | 对应的 `label_*` | orientation/side/overlap/rotation/offset/wrap/layout/fit/max_lines/min_size/segment 参数统一；`max_overlaps` 和固定 `size` 不变 |
| ribbon 的 `ribbon_gap`、`ribbon_ctrl_point` | `link_gap`、`link_ctrl_point` | 与 point link 共用名称，保留各自默认与几何 |
| `gene_offset = .1`、`feature_offset = .1` | `position = position_strand(offset = .1)` | 放置统一由 Position 负责 |
| `leader = "trunk"` | `leader = "radial"` | 删除无独立路由语义的兼容取值 |

旧逐链 offset 将负链值作为幅度再翻转；Position 的显式逐链值已经是有符号距离。
例如旧 `gene_offset = c("+" = .1, "-" = .2)` 应改为
`position_strand(offset = c("+" = .1, "-" = -.2))`。逐序列嵌套链设置同样转换负链值；
单个正标量直接迁移。不能把旧正幅度原样作为负链的正距离。

theme 的 `text`、`plot.*` 和通用 `legend.*` 移出显式签名，原来的命名调用仍可通过
`...` 使用，如 `theme_ggchord(text = element_text(size = 12), legend.position = "bottom")`。
IDE 是否展开上游 dots 由 IDE 决定；本包能保证自身 formals 仅列专属设置。
这些签名整理后的调用应使用参数名，不依赖旧的位置匹配顺序。

以下为完整教学输入，不代表真实科研数据：

```r
library(ggchord)
seqs <- data.frame(accver = c("A", "B"), length = c(1000, 1000))
genes <- data.frame(accver = "A", start = 100, end = 300,
  strand = "+", anno = "gene")
features <- transform(genes, type = "CDS", label = "feature")
links <- data.frame(qaccver = "A", saccver = "B", qstart = 100,
  qend = 300, sstart = 400, send = 600, pident = 90)
ggchord(seqs, links, genes) + geom_seq() +
  geom_gene(width = .05, position = position_strand(.1)) +
  geom_feature(data = features, shape = "block", width = .08) +
  geom_gene_label_repel(label_layout = "radial", label_fit = "wrap") +
  geom_link_ribbon(link_gap = .03, link_ctrl_point = c(0, 0))
```

- 固定样式优先：旧调用 `geom_link_line(aes(link_colour = group), colour = "red")`
  以前显示映射颜色，现在显示 red；要继续映射，写
  `geom_link_line(aes(link_colour = group))`。ribbon 的 alpha/colour/linetype、feature
  形状采用同一规则；省略固定参数不会以默认值压住映射。未删除参数，无删除版本。
- `geom_restriction_site(leader = "straight")` 现在使用直接连线，`leader = "elbow"`
  使用一个折点；这两种旧调用曾使用 radial 路由。需要保留旧图形时改用
  `leader = "radial"`；默认 radial 不变。直线可能穿越骨架或其他连线，
  不自动改为折线；没有承诺任意密度下无相交。
- `coord_chord(rotation = 45)` 调用与几何保持不变，帮助修正为从正 x 轴逆时针旋转；
  原来的 clockwise 文案与源码旋转矩阵不符。无需取反参数迁移。
- `scale_seq_position_continuous(breaks = ...)` 配置序列轴；
  `geom_seq_scale(length = 200)` 为 collinear 添加独立 bp 参照条。二者不互为别名，
  没有重命名或删除版本。
- gene 是共享 feature 构造内核的定向箭头快捷入口，仍保留 gene 独立美学和严格
  正负链要求。未知链继续使用 feature；未删除便利入口或22个第三方重导出。

## v0.16.0 历史修复记录

**状态：v0.16.0 本地开发完成，尚未正式发布。以下行为表描述该版当时事实；
v0.16.1 的参数迁移以前节为准。**

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
| F33 | 原始 label 等列仅因名字在保护名单中而被丢弃，gene_fill = label 报找不到列 | 几何尚无同名受保护字段时保留原列供映射 | 原有映射调用恢复；几何已经计算的受保护字段仍优先，未解决所有同名列冲突 | 不适用 |
| F34 | NULL 控制点帮助未区分坐标，错误暗示 chord 也从端点推导 | 帮助明确 chord 用 c(0, 0)，collinear 从端点推导 | 调用与图形不变；chord 想显式表达既有默认可写 ribbon_ctrl_point = c(0, 0) | 不适用，文档修正 |
| F35 | chord 变换后的轴 break 被当成原 bp，刻度标注与位置不一致 | 主/次 break 反变换后按原 bp 放置 | 原调用自动恢复正确位置，无需重编码 breaks；序列几何不变 | 不适用 |
| F36 | NULL 主/次刻度可能被自动刻度回退覆盖 | 显式空值保持空 tick 集合，轴线保留 | 继续用 breaks=NULL/minor_breaks=NULL 关闭相应刻度 | 不适用 |
| F37 | 默认连续ribbon alpha scale的NULL range导致palette调用错误 | NULL使用ggplot2原生默认range | 直接使用scale_ribbon_alpha_continuous()；已有显式range调用不变 | 不适用 |
| F38 | region过滤后可能附回错误源行的ID/坐标，factor坐标按等级码绘制 | 保留过滤前来源行，factor按数值标签解释 | 既有调用恢复真实坐标与来源，不需要重编码；被过滤的行不会替换保留几何的元数据 | 不适用 |
| F41 | 自定义引物同ID不同DNA时，后者被丢弃且name错配 | 完整aliases检测ID/sequence冲突并报原错误 | 不同oligo显式提供不同ID；同DNA别名仍共ID | 不适用 |

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

### F33：恢复缺失的输入元数据映射

修复 `aa02ad42827270a5bbbedff7e9b2559197407834` 前，下例的 `gene_fill = label`
会因展开后的几何缺少原始 label 列而报错；当前无需改名即可使用。数据角色映射
仍在输入表上求值，视觉映射随后在展开几何上求值。

```r
metadata_genes <- data.frame(id = "A", lo = 100, hi = 260, direction = "+",
  label = "alpha", annotation_tag = "shown", x = NA_real_)
metadata_original <- metadata_genes
metadata_plot <- ggchord(data.frame(accver = "A", length = 1000)) + geom_seq() +
  geom_gene(aes(accver = id, start = lo, end = hi, strand = direction,
    anno = label, gene_fill = label), data = metadata_genes) +
  geom_gene_label(aes(accver = id, start = lo, end = hi, strand = direction,
    anno = annotation_tag), data = metadata_genes) +
  scale_gene_fill_manual(values = c(alpha = "#113355")) + coord_collinear()
metadata_build <- ggplot2::ggplot_build(metadata_plot)
stopifnot(nrow(metadata_build$data[[2]]) > 0L,
  all(metadata_build$data[[2]]$gene_fill == "#113355"),
  all(is.finite(metadata_build$data[[2]]$x)),
  nrow(metadata_build$data[[3]]) > 0L,
  all(metadata_build$data[[3]]$label == "shown"),
  identical(metadata_genes, metadata_original))
```

输入 x 不替换已计算坐标，输入 label 不替换标签图层已经计算的 shown 文本。
保护只针对几何实际已具备的受保护列；本次修复不承诺同名原始字段与计算字段
总能同时用于视觉映射，也不改变固定样式与映射的优先规则。

### F35/F36：变换轴和关闭刻度使用标准写法

```r
axis_sequences <- data.frame(accver = "A", length = 1000)
axis_base <- ggchord(axis_sequences) + geom_seq()
axis_logged <- axis_base + scale_seq_position_continuous(
  breaks = c(10, 100), minor_breaks = NULL, transform = "log10")
axis_empty <- axis_base + scale_seq_position_continuous(
  breaks = NULL, minor_breaks = NULL)
axis_ticks <- get_chord_layout(axis_logged)$axis_ticks
axis_no_ticks <- get_chord_layout(axis_empty)
stopifnot(nrow(axis_ticks) == 2L,
  identical(as.character(axis_ticks$label), c("10", "100")),
  nrow(axis_no_ticks$axis_ticks) == 0L,
  nrow(axis_no_ticks$axis_lines) > 0L)
```

变换控制轴的 break 选择和标签；几何位置仍对应原始 bp，不能自行先把 breaks
改成 log 值。上述图形只用于验证轴设置，未改动输入序列长度或生物学注释。

F37 恢复无参 `scale_ribbon_alpha_continuous()` 的正常构造与映射；若需要确定的
输出范围，仍可显式写 `scale_ribbon_alpha_continuous(range = c(.2, .9))`。

F38 的 `geom_seq_region(data = regions)` 调用不变。factor 起止列应具有真实数值
标签；未知序列或无效坐标仍按既有布局规则过滤，保留行的 source_row 指向过滤前
输入位置，原表不被修改。没有将 region 的过滤策略推广成所有 geom 的共同政策。

### F41：自定义引物先明确ID与DNA的一一对应

旧输入若把不同DNA都标为同一个primer_id，现在会得到原有说明
`each primer_id must identify one distinct sequence`，不再静默丢失序列或错配name。
确为不同oligo时为各行提供不同ID；同一oligo的别名可共用ID和相同DNA。

```r
custom_primers <- data.frame(primer_id = c("oligo1", "oligo2"),
  name = c("one", "two"), sequence = c("ACGTACGT", "TTTTGGGG"))
primer_matches <- find_primer_bindings(c(target = "ACGTACGTAAAATTTTGGGG"),
  primers = custom_primers, set = "all", circular = FALSE,
  min_annealed_bases = 8)
stopifnot(all(c("oligo1", "oligo2") %in% primer_matches$primer_id))
```

没有自动生成替代ID、改原始DNA或设定新的删除版本；这是恢复已有检查的修复。

## 尚未实施的迁移


F34 已按用户选择保留历史默认并澄清帮助。F03 在 v0.16.2 保留任务相关默认；
F09–F15 的 v0.16.1 取舍见上文，不再列为待决定。F07 getter/export 仍归 v0.16.3。
其他后续契约若改变，
必须补入旧/新调用、替代方案、影响、错误提示、用户确认状态及确切删除版本；
不能将待确认的删除日期写成承诺。详细方案见 [决策台账](API-AUDIT-DECISIONS.md)。
