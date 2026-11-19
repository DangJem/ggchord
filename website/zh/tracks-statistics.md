{{code:courseblockbg}}


<a id="tracks"></a>

## 统计观测数或显示给定数值 {#tracks}

`bar` 按源位置统计观测数，`col` 显示已提供的数值。各 accession/window 在视图投影或
显示复制前独立计算 Stat 和数据 Position。省略 `track` 时，各层使用同一默认带。用 `track_band("depth")` 声明具名带，各层通过 `track = "depth"` 引用；声明可放在引用前后。`position` 选择 stack、dodge 等源数据调整。

从随包比对中选取两个参考，按 query 起点统计观测数。同一位置的多条记录会累加。

```r
ids <- c("MT108731.1", "MT118296.1")
seqs <- subset(seq_data_example, accver %in% ids)
hits <- subset(ribbon_data_example,
  qaccver %in% ids & saccver %in% ids & qaccver != saccver
)
ggchord(seqs) +
  geom_seq() +
  track_band("query_starts") +
  geom_track_bar(aes(x = qstart, accver = qaccver), data = hits, track = "query_starts") +
  coord_collinear()
```

{{figure:foundationcount}}

这里的 `geom_track_bar()` 也可以写作 `geom_track(geom = "bar")`，两种写法使用相同的起点和计数规则。


| 几何 | 默认 Stat | 输入含义 |
| --- | --- | --- |
| bar | count | x/pos 观测，可选 weight |
| col | identity | 已给定 y/track_value，区间或位置 |
| point、line、path、step、text | identity | 源 x 与局部 y，text 还需 label |
| tile、rect、area、ribbon、segment、linerange、errorbar | identity | 对应 ggplot2 原生范围字段 |

可用统计包括 identity/count/bin/summary、区间覆盖，以及[已声明的估计组合](tracks-statistics.md#estimates)。bin 需显式给出 binwidth 或 breaks。自定义 ggplot2 Stat 可从 source_row 读取本图层的源行号，并用 source_rows 列表返回各项结果的成员；省略成员时，导出将其标为未知。

源位置使用 x orientation，原生文字字号目前使用 mm，区间几何只支持完整轮廓。其他方向、字号单位或局部轮廓会明确报错。原生文字样式和路径箭头继续转发。

轨道使用 ggchord 坐标。固定量尺的样本分面按源表的分面列分配观测，每个样本内分别统计，共享具名带的数值域。用 `coord_collinear(seq_axis = "vertical")` 选择纵向投影。统计结果保留 dodge 等位置调整前的来源。

<a id="interval-coverage"></a>

## 按区间计算每个窗口的覆盖 {#interval-coverage}

`stat_coverage(data = intervals, binwidth = 4)` 默认画出每个窗口内至少被一条区间覆盖的位置比例。改用 `aes(y = after_stat(mean_depth))` 可查看平均深度。端点从 1 开始，两端都计入。跨窗区间分别贡献各窗的交集，平均深度对每条输入记录的交集单独计长。

以长度 10 的合成序列为例，第 1、2、3 行区间分别为 3–6、4–5、10–10。位置 4 和 5 各被两条记录覆盖。末窗只有两个位置，分母是 2。

| 窗口 | 不同覆盖位置数 | 覆盖率 | 平均深度 | 来源行 |
| --- | --- | --- | --- | --- |
| 1–4 | 2 | 2/4 | 3/4 | 1, 2 |
| 5–8 | 2 | 2/4 | 3/4 | 1, 2 |
| 9–10 | 1 | 1/2 | 1/2 | 3 |

有观测组的空窗返回零覆盖、零深度和空来源列表；没有观测的序列或样本不产生统计行。每条序列或已聚焦窗口从局部位置 1 开始分箱。也可用 `windows` 提供显式窗口，各窗分别统计交集，可互相重叠。`binwidth` 与显式窗口只选一项。区间和窗口的单位应与参考一致。默认拒绝越界区间，`bounds = "clip"` 先将区间裁至参考范围。多段和跨起点记录先拆成单段。

<a id="independent-bands"></a>

## 两条独立轨道与显示顺序 {#independent-bands}

同一 GC 既可用柱高也可用颜色编码。这里明确声明 gc_height 和 gc_colour 两个 ID；高度值域 0–1 与颜色值域 0–1 分别设置，物理带宽不同。第二幅只把颜色轨道 order 从 2 改为 0，使其靠近骨架。order 是有限数值优先级，不是测量值，编号不必从 1 开始；同优先级按 ID 排序。同一 ID 的各层沿用带的 side/order/width/gap，跨独立窗口不会隐式连接折线。

{{code:comparemultitrackdata}}
{{code:comparetracktwo}}
{{figure:comparetracktwo}}
{{code:comparetrackorder}}
{{figure:comparetrackorder}}

跨图比较同一测量值时，先用 [shared_track_limits()](../../man/shared_track_limits.Rd) 从两个图取得共同范围，再把该范围传给两边的轨道标度。

<a id="window-policy"></a>

## DNA 窗口与折线边界 {#window-policy}

此处 GC 为 `mean(bases %in% c("G", "C"))`：模糊碱基仍计入分母、对分子贡献为零；末窗口使用实际碱基数。当前 pUC19 不含模糊碱基，换输入时应明确并记录策略，不能假定自动忽略。局部聚焦裁剪区间并保留原汇总，不会重新计算缩短窗口的 GC。折线按输入样本及分组连接，`closed` 显式控制圆周闭合，`gap` 属于轨道的物理间隔。不能用闭合线暗示被省略源窗口之间存在连续观测。

<a id="estimates"></a>

## 原生源空间估计 {#estimates}

`geom_track()` 还会在**投影前**委托原生 ggplot2 执行 `smooth`、`density` 和 `ecdf`。它们是给定观测的估计，不是新增观测，按 accession/窗口及显式 group 独立计算。估计的 `source_rows` 记录原生缺失值处理后提交给估计器的组，不表示每点影响力、阳性 read support，也不是自定义模型内部进一步过滤的报告。

| 统计 | 已声明几何 | 输入和显式选择 | 数据 Position |
| --- | --- | --- | --- |
| `smooth` | `line`、`path`、`point`、`ribbon` | 源 `x` 与原数值 `y`；选择 method/formula，ribbon 需要 SE 边界 | `"identity"` 或原生 PositionIdentity |
| `density` | `line`、`path`、`point`、`area` | `x` 中的位置及可选权重；选择带宽，`after_stat(density)` 单位是源位置单位的倒数 | `"identity"` 或原生 PositionIdentity |
| `ecdf` | `line`、`path`、`point`、`step` | `x` 中的位置；`after_stat(ecdf)` 是累计比例；显式 `pad = FALSE` 避免无限源锚点 | `"identity"` 或原生 PositionIdentity |

density area 需指定 `position = "identity"`，因为 area 默认堆叠。`bounds` 用于控制源序列边界的裁剪。可投影的几何见本章支持表。

复用上文准备的真实 `bins` 和 `track_seq`。下面显式线性拟合作为 GC 比例的显示摘要，原生置信 ribbon 与均值线共享一个 Track：

{{code:courseblockbi}}

拟合使用给定的原始源 x/y，之后 Track 数值标度的变换只映射估计，不会暗中重拟合 log 空间模型；这与原生 ggplot 在 Stat 之前放置变换有别。如需在变换空间拟合，应事先显式准备、命名变换变量并声明单位。事件位置的核密度不是测序深度，ECDF 是累计分布，不是逐位点比对。

真实 pUC19 GC 窗口及 G/C 位置已与原生 ggplot2 数值和组来源逐项对照，并核对曲线/竖向投影及显式缺失值控制。`check-track-estimates.R` 保存数值报告和实际 CI/均值图。有限源边界及物理 band 约束仍适用；自动摆放不改变模型或制造观测。

<a id="track-focus"></a>

## 轨道与窗口使用相同源坐标 {#track-focus}

这里聚焦 pUC19 的 500–1500 bp，并按 interval 裁剪同一批实际 GC 窗口。focus_ggchord_track 不重新计算窗口均值，不在边界创造观测，也不会用一条线跨过省略的窗口。需要重新计算更细粒度 GC 时，应回到原 DNA 明确新窗口与分母。

{{code:comparetrackfocus}}
{{console:comparetrackfocus}}

## 对照同一组位置的计数与已汇总数值 {#count-and-col}

前例的 `bar` 统计每个 query 起点上的 HSP。下面先按位置计数，再用 `col` 画出相同的总数。`track` 命名物理轨道，`position` 决定数值在轨道内的排列。

```r
start_totals <- foundation_events |>
  dplyr::count(accver, pos, name = "hsp_count")
p_start_totals <- ggchord(data = foundation_seqs) +
  track_band("query_starts") +
  geom_seq() +
  geom_track(
    aes(x = pos, y = hsp_count), data = start_totals,
    geom = "col", track = "query_starts",
    position = "identity"
  ) +
  coord_collinear()
```
