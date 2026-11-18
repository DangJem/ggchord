## 迁移到 0.50.0

下表左列入口均在 0.50.0 移除。后续章节保留各历史版本的用法。

| 原调用或参数 | 0.50.0 写法 |
| --- | --- |
| `geom_*()` 曾忽略的拼写，如 `geom_seq(linewditth = 2)` | v0.50.0 起会报错。改写为 `geom_seq(linewidth = 2)`。统计专用参数随所选 Stat 传入，例如区间数据使用 `geom_track_line(stat = "coverage", binwidth = 4)`。 |
| 缩写形参，如 `position_offset(offse = .1)` 或 `guide_seq_axis(sequence = "A")` | 改用 `position_offset(offset = .1)` 和 `guide_seq_axis(sequences = "A")`。本次重构的入口要求写全具名参数，经 `...` 转发也适用。按位置传值仍可写 `position_offset(.1)`。 |
| 显式传入一个 `width` 同时设置外部 Stat 的邻域半径和 Geom 的柱宽 | 将 Stat 的半径参数改名为 `radius`，分别设置两项，例如 `geom_track_col(stat = MyStat, radius = 4, width = 2)`。 |
| `ggchord(seq_data = ref, gene_data = f, ribbon_data = h)` | 从 `ggchord(ref)` 开始，分别加 `geom_feature(data = f, preset = "gene")` 和 `geom_link_ribbon(data = h)`。也可写 `ggchord(prepare_ggchord_data(ref, ribbon_data = h, gene_data = f))`。有多张可用的准备表时在图层中指定数据。 |
| 用 `position_track("depth", width = w)` 指定物理带 | 加一次 `track_band("depth", width = w)`，在同带各层写 `track = "depth"`。默认带可省略轨道 ID。堆叠、并排和微移交给各层的 `position`。 |
| `geom_seq()` 的 `seq_order`、`seq_orientation`、`seq_gap`、`length_mode` | 添加 `layout_sequences(order = ids, orientation = directions, gap = gaps, length_mode = "equal")`。骨架线条与固定样式继续放在 `geom_seq()`。 |
| `position_plasmid(offset)` 和 `position_feature_stack(side = ...)` | 偏移改用 `position_offset(offset)`。原 `side = "outside"` 用 `base_position = position_offset(-0.1), direction = "inward"`。原 `"inside"` 用 `position_offset(0.1)` 及 `"outward"`。原 `"strand"` 对应默认分道。 |
| `geom_gene()`、`geom_feature_plasmid()` | 使用 `geom_feature(preset = "gene")` 或 `geom_feature(preset = "plasmid")`，传入相同图层数据。 |
| `geom_feature_label_repel()` | 使用 `geom_feature_label(placement = "repel")`。原先的放置策略交给 `label_layout`，原 `label_fit` 改为 `overflow`。多个实体层均可匹配时用 `anchor_layer = "features"` 指定目标。 |
| `scale_seq_position_continuous(sequences = ids)` | 改用 `scale_seq_position_continuous(guide = guide_seq_axis(sequences = ids))`。 |

使用阶段导出的轨道统计前，先检查 `layer_stage_status$status`。统计成功并返回零行时为
`available`。Stat 失败时，统计及几何阶段均标为 `failed`，`reason` 保留原因。
部分面板失败时，结果仍可包含成功面板的行。

运行[迁移示例](examples/migrate-v0.50.0.R)可绘制双样本比较，并在结果表中检查每个窗口的来源行。两面板共用 10 bp 参考和数值域。第二张图先按完整特征表分道，再对选中的特征加标签。

## 0.29.0 数据与映射

| 工作流 | 当前写法 |
| --- | --- |
| 独立清洗与反向连接 | 默认排序保留相对方向；需要方向端点时对原表使用 `reversed_interval = "keep"`。 |
| 自定义 BLAST 列 | 未指定类型的自定义元数据保留文本。可用 `colClasses = c(score_custom = "numeric")` 读入数值测量，也可在读入后用 `mutate(score_custom = as.numeric(score_custom))` 转换。其他自定义标识保留原拼写。 |
| 轨道身份 | `aes(accver = id)` 的结果须为单值或逐行等长。原先较短的向量会循环补齐，现在报错。 |
| 科学示例 | 真实数据使用 `strand = source_strand`，GFF3 可用 `unstranded = "keep"` 保留未知链。 |
| 注释方向 | 同时提供 `strand` 和 `directionality`（或别名 `direction`）时，按报错行核对来源，或只把可信字段交给图层。 |
| 函数形式的图层数据 | 每次重新绘制时，每个图层读取一次。方向推断、几何及随布局导出的图层输入使用同一批记录。 |
| 准备结果 | 将原本逐项传递的默认表简写为 `ggchord(prepared)`。命名集合由图层选取，如 `geom_feature(data = prepared$features$cds)`。在 `feature_data` 上重命名列用 `mapping = list(feature = c(strand = "source_strand"))`。 |
| 定量轨道 | 折线和瓦片均可用 `track = position_track("gc")` 指定物理轨道。普通折线支持 `aes(x = pos, y = gc)`。闭合圆形折线另设 `closed = TRUE`，区间瓦片继续使用起止位置。观察计数用 `geom_track_bar()`，已有数值高度用 `geom_track_col()`。 |
| 自动标签 | 在 `geom_feature_label_repel()` 中用 `placement = "callout", overflow = "ellipsis"` 允许外置省略，或用 `placement = "feature", overflow = "shrink"` 按特征空间缩小字号。已有的 `label_layout`、`label_fit` 组合仍按原规则解释。一次调用只选一组。 |
| 追溯图层数据 | `export_ggchord_layout(p, data_stage = c("source", "mapped", "statistic", "geometry"))$layer_data` 按图层 ID 返回所选阶段。source 是适配后的输入。原生轨道的映射在 ggplot2 内执行，其 mapped 为空，可查 statistic。数据函数在本次导出重新求值。 |
| 主题与图例 | 全图字体用 `theme(text = element_text(family = "sans"))`，单个连接带图例标题用 `labs(ribbon_fill = "Identity")`，序列标签字号用 `theme(ggchord.seq.label = element_text(size = 10))`。 |
| 调整连接带 | 先用 `geom_link_ribbon()` 查看默认几何。端点留白用 `link_gap`，局部避让用 `link_avoid = "smooth"`，需要改变中段曲线时再指定 `link_ctrl_point`。 |
| 起始方向 | 依次对齐 12、3、6、9 点钟方向时，`coord_circular(rotation = ...)` 分别取 0、−90、180、90，`coord_radial_synteny(start_angle = ...)` 分别取 90、0、270、180 度。前者的 `direction` 控制序列坐标递增方向，后者的 `clockwise` 控制后续辐条次序。 |
| 连接带统计 | `geom_link_ribbon(stat = "ribbon_bundle", bins = 20, weight = "count")` 等价于原有 `stat_ribbon_bundle(bins = 20, weight = "count")`。保留单条连接并按局部密度调整透明度可用 `geom_link_ribbon(aes(ribbon_alpha = after_stat(density)), stat = "ribbon_density", bins = 20, weight = "count")`。只有 bundle 接受 `min_bundle`。 |
| 特征与轨道位置 | 默认分道可写 `geom_feature(position = "feature_stack")`。自定义间距用 `position_feature_stack(spacing = 0.15)`。定量层用 `track = position_track("depth")` 声明物理轨道，另用 `position = ggplot2::position_dodge(width = 10)` 调整源坐标。传错对象的调用现在会报错。 |
| 轨道统计与来源 | 例如 `geom_track(aes(x = pos), data = events, geom = "bar", stat = ggplot2::StatBin, binwidth = 10, track = position_track("events"))`。导出的 `source_rows` 在移动几何之前确定，每个汇总值和加权计数只列出相应的有效行。 |
| 区间覆盖 | 用 `geom_track_bar(data = intervals, stat = "coverage", binwidth = 100, track = position_track("coverage"))` 绘制各窗口的并集覆盖比例。改用 `aes(y = after_stat(mean_depth))` 可让重叠记录分别贡献深度。 |
| 分道方向 | `base_position` 决定第一道，`direction` 决定后续道。旧 `side = "outside"` 等价于 `base_position = position_plasmid(-0.1), direction = "inward"`，旧 `side = "inside"` 等价于 `base_position = position_plasmid(0.1), direction = "outward"`。默认按链放置，旧调用继续有效。 |
| 分阶段配色 | 在定量轨道上，`aes(fill = stage(start = condition, after_scale = scales::alpha(fill, 0.5)))` 先按 condition 分组并训练轨道填充标度，再调整已得到的颜色。原始列名可含 after_stat 或 after_scale。 |

## 0.19.0：静态跨面板组合

既有单图调用及默认行为保持不变。先创建独立图形，再用
`ggchord_fusion(center) + fusion_panel(detail, "detail") +
fusion_ribbon(links, qpanel="center", spanel="detail")`组合。
继续使用标准q/s字段；只有组合路由需要面板或实例选择。
跨面板层使用标准`aes(fill=...)`与`scale_fill_*()`，角色标度仍属于各自子图。
通过`ggsave()`保存，用`export_ggchord_fusion()`核对物理几何，不调用`ggplot_build()`。
要求区间完整可见且单位相同；这是增量接口，没有重命名或移除旧API。

[English](MIGRATION.md) · [现代汉语](MIGRATION-Hans.md)

# 迁移到 ggchord 候选 API

## 0.18.0 的径向共线性

在支持的序列、Feature和Track组合中使用`coord_radial_synteny()`，默认保持源长度
比例；只有显式`geom_seq(length_mode = "equal")`才归一显示长度。环形邻接包含
已有的末条/首条关系，不补造闭合比对。新视图的完整区间前缘使用按输入方向排序的
`[min-1,max]`边界，已有坐标约定不变，位置轴从1开始。径向反查须提供所选骨架上
的x/y表。不支持共享直线比例尺、branch、多窗口拼行及圆图专有注释，会明确报错。
双语view-expansion教程提供真实来源、物理Track、实例路由和完整支持矩阵。


## 0.18.0 的 HSP Dotplot

使用相同的序列表和 Ribbon 表，组合 `geom_alignment_segment()` 与
`coord_dotplot()`，无需 `geom_seq()`。点连接继续用 `geom_link_line()`；
Alignment 样式复用 `link_*` 映射及 `scale_link_*()`。端点必须已规范化方向，
`read_paf()` 已完成此步骤；不能直接把升序 target 区间加负链当作正向线段。
未知方向会报告并省略，没有精细比对轨迹时不生成逐碱基路径。
位置反查指定 `map_ggchord_position(..., axis = "x"/"y")`；聚焦后的原始来源
坐标用 `source = TRUE`。窗口轴默认仍显示窗口内相对坐标。已有弦图、单圆图和
共线性代码无需因此改写。


## 0.17.0 的统一注释输入

推荐用 `ggchord(seq_data = seqs, feature_data = annotations)` 提供基因和其他特征。
`geom_feature()` 默认使用全表；没有显式图层数据时，基因层依次从 `feature_class`、
`type`、`role` 中选择 `gene`/`CDS`，tRNA/rRNA 使用通用特征层。缺失链按未知处理，
通用层也接受双向注释。

`gene_data` 继续兼容，包括第三位置参数；每次构造至多提示一次迁移，渲染不重复提示，
尚未安排删除版本。无类型的旧表保留原基因语义。迁移时先给方向已知的旧基因表添加
`type = "gene"`，再传给 `feature_data`；也可直接使用 `geom_gene(data = genes)`。
两个构造输入不能同时提供。图层数据优先于全局表；组合基因与其他特征时显式筛选
互斥子集，避免把同一全表绘制两遍。

```r
library(ggchord)
library(dplyr)
data(seq_data_example)
data(gene_data_example)
annotations <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand, type = "CDS")
p <- ggchord(seq_data_example, feature_data = annotations) +
  geom_seq() +
  geom_gene()
```

`export_ggchord_layout()` 的 `layer_sources` 提供注释来源身份；原几何 `source_row`
仍指向解析后的图层输入行。默认注释层共享构造输入的来源 ID；显式图层表与函数采用
图层局部身份，不按数值内容猜测来源，也不会去重删除相同的真实源记录。

## 0.17.0 的轨道数值与序列坐标轴

`geom_track_bar()` 现在按原生 bar 语义统计观测数；`geom_track_col()` 以
`stat = "identity"` 显示已给定的数值高度。两个函数均未删除。迁移旧数值柱时，
将函数名换成 `geom_track_col()`，保留数据、数值映射及物理 band。

| 旧用法 | 替代写法 | 变化 |
| --- | --- | --- |
| 用 `geom_track_bar(aes(track_value = gc), data = bins, ...)` 显示给定值 | `geom_track_col(aes(track_value = gc), data = bins, ...)` | 数值高度改用 col，bar 统计行数。 |
| 用 `geom_track(..., geom = "bar")` 显示给定值 | `geom_track(..., geom = "col")` | 统计观测数时用 bar 并映射 `aes(x = position)`。 |
| 示例中显式从0开始的序列刻度向量 | `scale_seq_position_continuous(n.breaks = 3)` | 自动序列轴从1开始，优先保留真实首末位置。 |

`track = position_track(...)` 声明物理 band，`position` 调整数据；此前
`position = position_track(...)` 的物理 band 写法继续兼容。源坐标窗口保留原始
基因组区间，仍可显式设置 breaks/labels；闭圆接缝合并1/长度。Identity 标尺仍从0开始，
采用0/50/80/90/95/100；其他自动托管的数值pident映射也使用同一域与刻度，
显式用户标度继续优先。

## 0.16.7 跨行连接带边框

`geom_link_ribbon()` 对共线性非相邻行采用 `nonadjacent_linewidth = 0.4`，
独立于普通固定 `linewidth`。设为 `nonadjacent_linewidth = 0` 时，
跨行带沿用图层的普通线宽。显式线宽映射由映射和标度控制。
该参数只改变线宽，不改变默认的虚线线型。

## 0.16.7 圆图起点

`coord_circular()` 保持默认 `rotation = 0`，将起点放在 12 点钟。旋转以此为
基准逆时针计算：90 对应九点钟，-90 对应三点钟。要保留 v0.16.6 或此前
v0.16.7 开发快照中显式指定角度的效果，将原角度减去 90 度；省略 rotation
的调用保持当前开发版的默认绘图效果。坐标默认仍顺时针增加；弦图旋转契约不变。

## 范围与时间 {#scope}
下列写法是当前候选替代，不是 v1.0.0 长期兼容承诺。0.16.0–0.16.4 阶段成果已随 v0.16.5 交付。
当前安装标签为 v0.16.6。RC 审计保留
灵活参数、拼写别名、重导出及 getter/export 入口，没有新增弃用。
表中的旧调用用于解释迁移，不应作为当前例子执行。
历史证据只提供删除上界时，不猜测精确首删 tag。

## 替代表 {#table}
| 旧写法（勿执行） | 当前写法 | 删除 / 守卫 | 含义 |
| --- | --- | --- | --- |
| `geom_ribbon()` | `geom_link_ribbon()` | ≤ 0.16.1 | 连接图层改名，不影响 ggplot2::geom_ribbon()。 |
| `ggchord(..., rotation = 20)` | `ggchord(...) + coord_chord(rotation = 20)` | ≤ 0.16.1 | 全局坐标使用 + 组合。 |
| `geom_feature(feature_shape = "block", feature_width = .05)` | `geom_feature(shape = "block", width = .05)` | 0.16.1 | 映射角色 feature_shape 保留。 |
| `geom_gene(gene_offset = c("+" = .1, "-" = .2))` | `geom_gene(position = position_strand(c("+" = .1, "-" = -.2)))` | 0.16.1 | 显式 Position 值带符号。 |
| `geom_feature(feature_offset = .1)` | `geom_feature(position = position_strand(.1))` | 0.16.1 | 使用 Position 对象。 |
| `geom_gene(gene_width = .05)` | `geom_gene(width = .05)` | 0.16.1 | 固定几何统一使用 width。 |
| `geom_link_ribbon(ribbon_gap = .01)` | `geom_link_ribbon(link_gap = .01)` | 0.16.1 | 使用 link_gap，控制点同样使用 link_ctrl_point。 |
| `geom_gene_label_repel(gene_label_layout = "aligned")` | `geom_gene_label_repel(label_layout = "auto")` | 0.16.1 | 同时迁移参数名和取值。 |
| `geom_gene_label_repel(label_rotation = 90)` | `geom_gene_label(label_rotation = 90)` | 0.16.4 guard | 原先忽略的手动控制现报错，也可删除它使用自动布局。 |
| `geom_feature_label_repel(label_orientation = "horizontal")` | `geom_feature_label(label_orientation = "horizontal")` | 0.16.4 guard | 未知链 feature 使用固定 feature 标签层。 |
| `geom_restriction_site(leader = "trunk")` | `geom_restriction_site(leader = "radial")` | 0.16.1 | radial 默认值保持。 |
| `theme_ggchord(base_size = 10)` | `theme_ggchord() + theme(text = element_text(size = 10))` | ≤ 0.16.1 | 文字外观属于 theme 元素。 |

## 运行当前替代 {#replacement}
在新会话中运行下面的准备代码，再执行当前替代写法：


```r
library(dplyr)
seqs <- seq_data_example |>
  select(accver, length)
links <- ribbon_data_example |>
  select(qaccver, saccver, qstart, qend, sstart, send)
# 随包示例用 source_strand 单独保存原始注释链。
genes <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand)
```


```r
library(dplyr)
data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19 |>
  select(accver, length)
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source |>
  select(accver, start, end, strand, type, anno)
```



```r
p_migrated <- ggchord(plasmid_seq) +
  geom_seq() +
  geom_feature(data = features, shape = "block", width = .05) +
  geom_feature_label(data = features, label_rotation = 0) +
  coord_circular()
p_offsets <- ggchord(seqs, gene_data = genes) +
  geom_seq() +
  geom_gene(position = position_strand(c("+" = .1, "-" = -.2)))
```


旧负链 offset 幅度在显式 Position 中需要负号。映射仍使用
`aes(feature_shape = type)` 配合 shape 标度，固定参数则为 `shape`。
自动标签使用 `label_layout`，手动几何控制使用相应固定标签层。
`label_layout = "aligned"` 改为 `"auto"`。已删除名称明确报错，不保留静默执行路径。

## 保留的支持 {#retained}
五个 colour/color 拼写别名、全部上游重导出、gene/注释/plasmid 快捷入口及 getter/export
两类几何接口均保留。各坐标默认差异和已记录结果字段保持。NULL/NA/缺参语义应查阅
接口说明 (`help("ggchord_api_contract", package = "ggchord")`)和当前参考；NULL 不统一表示关闭。

## v0.16.7 的弦图默认旋转

弦图默认旋转改为 0 度。需要 45 度布局时，显式添加 `coord_chord(rotation = 45)`。原生蛋白图声明 `sequence_unit = "aa"`，并提供完整蛋白长度。

## v0.17.0 基础架构迁移

`geom_track_bar()` 现在按源位置统计观测数（`stat="count"`）。旧调用若已提供数值高度，
请改为 `geom_track_col()`，原映射与 `position=position_track("signal")` 可以保留；
也可以使用 `track=position_track("signal")` 并单独指定数据 Position。
通用 `geom_track()` 支持 point/line/path/step/col/bar/tile/rect/area/ribbon/segment/
linerange/errorbar/text 和 identity/count/bin/summary；bin 必须显式给出 binwidth 或
breaks。任意第三方 Geom/Stat 扩展不属于本版。

prepare/clean/focus/segment 工作流接受 `feature_data` 并返回同名结果；旧 gene 调用继续
返回 `gene_data`。两者不能同时提供。源记录、未知链和元数据保留。

`geom_seq(seq_order=c(A_high="A", B="B", A_low="A"))` 声明显示副本。
各连接层先筛选真实 HSP，再指定 `between=c("A_high","B")` 或 `between=c(2,3)`，
路由必须唯一。提供稳定的 `source_hsp_id` 可检查跨层重复显示；确实需要重复时显式
设置 `allow_repeated=TRUE`。未提供时来源身份只在层内有效，不根据相同值推断身份。
同一图对象可以切换 Chord/Collinear。Sankey 自己声明 query 左列、subject 右列副本，
不将 between 解释为多阶段流量。

`coord_collinear(seq_axis="vertical")` 使用纵向骨架，反查坐标传入最终 y。
`coord_sankey()` 默认过滤同源实体自身连接，保留输入的互反记录；首批 MVP 不支持
共享分支。坐标映射时可从导出的 `display_instances` 选择无歧义实例。

`geom_seq(length_mode="equal")` 只改变显示分配，真实长度、坐标和局部轴仍保留；
禁止具有误导性的共享绝对比例尺。默认刻度显示真实首位 1 与末位，闭圆接缝合并两者；
显式 breaks 始终优先。

Ribbon 默认恢复此前的 Identity 填充并保留完整区间端面。
显式 `gradient = TRUE` 才按曲线从 query 序列色渐变至 subject 序列色。
显式 `tip="subject"` 才收尖：箭肩保留平移后的目标端点跨度，尖点仅是装饰标记，
不代表唯一匹配坐标。
`tip="none"` 保留完整端面；`aes(ribbon_fill=pident)` 选择 Identity 模式，固定域为
0–100，0/50/80/90/95/100 区间按数值比例显示。固定 fill、用户映射和用户标度优先。

渐变填充沿真实曲线使用设备分辨率网格，轮廓保持向量路径，避免透明分片接缝，并保留弯曲的来源端面。

自动 Identity 图例不再把通用图例色块尺寸当成整条数值标尺的长度，而是按边界文字测量和区间比例预留长度；独立设置 Ribbon 图例尺寸或显式 guide 可覆盖它。圆图 parallel 坐标轴文字逐字符沿弧排布；水平或数值角度仍保持直线文字。
