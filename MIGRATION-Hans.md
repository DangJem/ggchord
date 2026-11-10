[English](MIGRATION.md) · [现代汉语](MIGRATION-Hans.md)

# 迁移到 ggchord 候选 API

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

Ribbon 默认按曲线从 query 序列色渐变至 subject 序列色，默认保留完整区间端面。
显式 `tip="subject"` 才收尖：箭肩保留平移后的目标端点跨度，尖点仅是装饰标记，
不代表唯一匹配坐标。
`tip="none"` 保留完整端面；`aes(ribbon_fill=pident)` 选择 Identity 模式，固定域为
0–100，0/50/80/90/95/100 区间按数值比例显示。固定 fill、用户映射和用户标度优先。

渐变填充沿真实曲线使用设备分辨率网格，轮廓保持向量路径，避免透明分片接缝，并保留弯曲的来源端面。
