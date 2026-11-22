# API 版本选择

需要使用旧 API 的代码请选用 v0.29.0，从 v0.50.0 起，接口与旧版不兼容。

<details>
<summary>维护者内部变动记录</summary>

以下档案按版本保存接口移除项、替代调用和核验过的行为，供维护者追溯兼容性决策和复核历史审计。

## 迁移到 0.50.0

下表记录 0.50.0 移除的调用和发生变化的行为。后续章节注明历史示例适用的版本。

| 原调用或参数 | 0.50.0 写法 |
| --- | --- |
| 原先两个未命名图层共享颜色或填色标度 | 在两层分别写 `id = NULL` 或 `id = NA` 即可共享默认标度。添加普通显式标度也会同时作用于两层。省略 ID 可分别训练标度，通过具名图层和 `scale_for()` 可单独指定标度。 |
| 原先省略 `anchor_layer`，直接使用图级数据放置独立标签 | 设置 `anchor_layer = NULL` 可绘制独立标签。省略时标签跟随唯一兼容目标；有多个候选时，为目标图层命名并在 `anchor_layer` 中填写该名称。 |
| 特征用 `aes(accver = ref, fill = accver)` 时，填色曾跟随重映射后的序列标识 | 各个起始映射都读取图层源列。保留此式可按原始序列标识着色。若要按绘制所用的参考着色，写 `aes(accver = ref, fill = ref)`。 |
| 点轨道按类别使用不同形状，或只需隐藏填色图例 | 用 `aes(shape = kind)` 映射类别，再用 `scale_shape_manual(values = c(a = 21, b = 24))` 指定符号。在轨道层设 `show.legend = c(fill = FALSE)` 可隐藏该层对填色图例的贡献。其他图层仍可提供该图例。 |
| 在同一高度轴混用覆盖率和平均深度 | 两种度量分带绘制。自定义高度可声明 `measure = "read_depth"`。跨图比较时将 `shared_track_limits(p1, p2)` 得出的范围设给各自的轨道标度。 |
| 分样本观察各面板的局部数值变化 | 使用 `scale_track_value_continuous("depth", facet_domain = "panel")` 并保留局部轴。默认在各样本间使用共同的高度范围。 |
| 按覆盖率给整窗着色 | 使用 `geom_track_tile(aes(fill = after_stat(coverage)), stat = "coverage", binwidth = 4, data = intervals)`。末窗按实际宽度绘制。 |
| 将双端记录分派到面板 | 原始连接带或点连接可选择源列，例如 `facet_wrap(vars(saccver))`。Dotplot 写明两端，例如 `facet_grid(rows = vars(qaccver), cols = vars(saccver))`。 |
| 替换骨架表中有重复键，或把两个不同的聚焦区间都命名为 W | 每个参考 ID 保留一行，不同窗口给不同的 `locus_id`。重复显示同一序列用 `layout_sequences(order = c(left = "A", right = "A"))`。 |
| 计数轨道的样本列也叫 `count` | 把样本键复制到 `sample`，再用 `facet_wrap(vars(sample))`；导出的 `count` 保留计数含义。 |
| 查询源阶段时也读取默认几何 | 需要轨道几何时同时给 `include = "track"` 和 `data_stage = "source"`，也可请求源和几何两个阶段。 |
| 按显示边缘提供区域高亮 | 改用一基闭端点。单位置高亮占一个序列单位。限制位点需使用已知参考上的 bp 坐标。 |
| 叠加的 `geom_seq()` 层为同一序列声明了不同的长度或单位 | 各层使用同一张参考表，分别设置样式。需要比较不相容的参考时分别构图。 |
| 轨道、区域高亮或点连接声明 bp，骨架为 aa，或反过来 | 提供骨架单位下的坐标。核酸区间要画到蛋白质参考上，可通过对应的 CDS 投影，也可在各自单位的参考上分别绘图。 |
| 轨道表达式依赖重映射后的序列 ID，或筛选后的空表缺少映射列 | 先在源表中生成新 ID 列，让各个映射直接引用它。空图层可由列齐全的输入表筛为零行，以保留列名和类型。 |
| `geom_track_col(aes(fill = kind), fill = "red", stat = "coverage", binwidth = 10)` | 要合并区间，去掉填色映射。要分组统计并统一画红色，使用 `aes(group = kind)` 和 `fill = "red"`。 |
| 依赖两张独立特征表的行内容相同来绑定标签 | 向实体和标签传同一张表，或提供匹配的 `feature_id` 及来源区间。筛选、重排后的子集按键匹配。 |
| `geom_*()` 曾忽略的拼写，如 `geom_seq(linewditth = 2)` | v0.50.0 起会报错。改写为 `geom_seq(linewidth = 2)`。统计专用参数随所选 Stat 传入，例如区间数据使用 `geom_track_line(stat = "coverage", binwidth = 4)`。 |
| 缩写形参，如 `position_offset(offse = .1)` 或 `guide_seq_axis(sequence = "A")` | 改用 `position_offset(offset = .1)` 和 `guide_seq_axis(sequences = "A")`。本次重构的入口要求写全具名参数，经 `...` 转发也适用。按位置传值仍可写 `position_offset(.1)`。 |
| 显式传入一个 `width` 同时设置外部 Stat 的邻域半径和 Geom 的柱宽 | 将 Stat 的半径参数改名为 `radius`，分别设置两项，例如 `geom_track_col(stat = MyStat, radius = 4, width = 2)`。 |
| `ggchord(seq_data = ref, gene_data = f, ribbon_data = h)` | 从 `ggchord(ref)` 开始，分别加 `geom_feature(data = f, preset = "gene")` 和 `geom_link_ribbon(data = h)`。已有来源分类的主表可写 `ggchord(prepare_ggchord_data(ref, ribbon_data = h, feature_data = f))`。有多张可用的准备表时在图层中指定数据。 |
| 准备、校验、清洗、聚焦、分段函数的 `gene_data`，以及准备函数的 `mapping$gene` | 0.50.0 的第三个参数为 `feature_data`。自定义列可写 `prepare_ggchord_data(ref, feature_data = f, mapping = list(feature = c(strand = "source_strand")))`。主注释输出、转换记录、校验行数分别读取 `$feature_data`、`$report$conversions$feature` 和 `$data_summary$n_feature`。附加表仍按名称独立保存，例如 `features = list(cds = cds, sites = sites)`。 |
| 已保存且含非空 `$gene_data` 的准备对象，或将无分类表按第三位置传入并继承基因预设 | 用原输入重建。基因预设继承时，依次查找 `feature_class`、`type`、`role`，按第一个存在的字段选 gene/CDS。原始资料确认所有记录均为基因时，无分类表可写 `feature_data = dplyr::mutate(f, type = "gene")`。要直接绘制有方向且无分类列的基因表，可写 `ggchord(ref) + geom_seq() + geom_feature(data = f, preset = "gene")`。继承记录的来源标识前缀为 `feature:`，显式图层使用该层自己的行标识。 |
| 用 `position_track("depth", width = w)` 指定物理带 | 加一次 `track_band("depth", width = w)`，在同带各层写 `track = "depth"`。默认带可省略轨道 ID。堆叠、并排和微移交给各层的 `position`。 |
| `geom_seq()` 的 `seq_order`、`seq_orientation`、`seq_gap`、`length_mode` | 添加 `layout_sequences(order = ids, orientation = directions, gap = gaps, length_mode = "equal")`。骨架线条与固定样式继续放在 `geom_seq()`。 |
| `position_plasmid(offset)` 和 `position_feature_stack(side = ...)` | 偏移改用 `position_offset(offset)`。原 `side = "outside"` 用 `base_position = position_offset(-0.1), direction = "inward"`。原 `"inside"` 用 `position_offset(0.1)` 及 `"outward"`。原 `"strand"` 对应默认分道。 |
| `geom_gene()`、`geom_feature_plasmid()` | 使用 `geom_feature(preset = "gene")` 或 `geom_feature(preset = "plasmid")`，传入相同图层数据。 |
| `geom_feature_label_repel()` | 使用 `geom_feature_label(placement = "repel")`。原先的放置策略交给 `label_layout`，原 `label_fit` 改为 `overflow`。多个实体层均可匹配时用 `anchor_layer = "features"` 指定目标。 |
| `scale_seq_position_continuous(sequences = ids)` | 改用 `scale_seq_position_continuous(guide = guide_seq_axis(sequences = ids))`。 |

使用阶段导出的轨道统计前，先检查 `layer_stage_status$status`。统计成功并返回零行时为
`available`。Stat 失败时，统计及几何阶段均标为 `failed`，`reason` 保留原因。
部分面板失败时，结果仍可包含成功面板的行。
普通构图遇到轨道 Stat 错误会停止。可用 `get_chord_layout(p, build = FALSE)` 读取上次成功的布局，修正输入后再构图。
零行输入也保留内置统计的已知计算列，可照常选取 `count`、`prop` 等字段。

[审计示例](examples/migrate-v0.50.0.R)记录双样本比较及各窗口的来源行。两面板共用 10 bp 参考和数值域。第二张图先按完整特征表分道，再对选中的特征加标签。

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

本节的构造参数和示例适用于 v0.17.0。当前版本的等价调用见 [v0.50.0 特征与基因示例](#replacement)。

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

## v0.16.5–v0.16.6 的迁移范围 {#scope}
0.16.0–0.16.4 阶段成果已随 v0.16.5 交付。v0.16.6 发布候选审计保留了灵活参数格式、
拼写别名、重导出及布局查询和导出入口。下表记录这两个版本的替代写法。历史证据仅能
确定删除上界时，删除栏只给出该上界。本地版本请参照篇首的 v0.50.0 迁移表和下方示例。

## 历史替代表 {#table}
| 旧写法 | v0.16.5–v0.16.6 的替代写法 | 删除 / 守卫 | 含义 |
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

## 在 v0.50.0 运行等效示例 {#replacement}
在新 R 会话中使用 ggchord v0.50.0 和 dplyr 运行下列示例。

```r
library(ggchord)
library(dplyr)

seqs <- seq_data_example |>
  select(accver, length)
# 随包示例用 source_strand 单独保存原始注释链。
genes <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand)

plasmid_seq <- plasmid_example_pUC19 |>
  select(accver, length)
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source |>
  select(accver, start, end, strand, type, anno)

p_migrated <- ggchord(plasmid_seq) +
  geom_seq() +
  geom_feature(data = features, shape = "block", width = .05) +
  geom_feature_label(data = features, label_rotation = 0) +
  coord_circular()

p_offsets <- ggchord(seqs) +
  geom_seq() +
  geom_feature(
    data = genes,
    preset = "gene",
    position = position_strand(c("+" = .1, "-" = -.2))
  )

p_migrated
p_offsets
```

旧负链 offset 幅度在显式 Position 中需要负号。按类别控制形状时，使用
`aes(feature_shape = type)` 配合形状标度。固定几何使用 `shape`。
自动标签用 `label_layout` 排布，手动几何控制用相应固定标签层。
`label_layout = "aligned"` 改为 `"auto"`。

## v0.16.6 审计保留的接口 {#retained}
v0.16.6 审计保留了五个 colour/color 拼写别名、上游重导出、gene/feature/plasmid
快捷入口，以及布局查询和导出接口。v0.50.0 用 `geom_feature(preset = "gene")`
和 `geom_feature(preset = "plasmid")` 承接原快捷入口。各入口对 `NULL`、`NA`
及省略参数的处理见 `help("ggchord_api_contract", package = "ggchord")` 和相应函数帮助。

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

</details>
