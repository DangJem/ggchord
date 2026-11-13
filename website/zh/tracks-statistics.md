{{code:courseblockbg}}

{{prepare:trackdata,tracks,foundationdata}}

## 统计观测数或显示给定数值 {#tracks}

`bar` 按源位置统计观测数，`col` 显示已提供的数值。各 accession/window 在视图投影或
显示复制前独立计算 Stat 和数据 Position。`track=position_track(...)` 声明物理 band，
`position` 独立选择 stack、dodge 等原生数据调整。旧的
`position=position_track(...)` 物理 band 写法继续接受。

{{code:foundationcount}}

{{figure:foundationcount}}

下面的快捷入口在相同物理轨道中统计同一组观测。

{{code:courseblockbh}}


| 几何 | 默认 Stat | 输入含义 |
| --- | --- | --- |
| bar | count | x/pos 观测，可选 weight |
| col | identity | 已给定 y/track_value，区间或位置 |
| point、line、path、step、text | identity | 源 x 与局部 y，text 还需 label |
| tile、rect、area、ribbon、segment、linerange、errorbar | identity | 对应 ggplot2 原生范围字段 |

可用统计包括 identity/count/bin/summary，以及已声明的 [smooth/density/ecdf 估计组合](tracks-statistics.html#estimates)。bin 必须显式给出 binwidth 或 breaks。
支持原生 group 与 after_stat 映射；数值轴保留统计量单位。任意第三方 Geom/Stat 插件不在支持范围内，不支持的组合明确拒绝。既有 track line/tile 快捷层保留其原有行为。

源位置使用 x orientation，原生文字字号目前使用 mm，区间几何只支持完整轮廓。其他方向、字号单位或局部轮廓会明确报错。原生文字样式和路径箭头继续转发。

首批原生 Track 要求单面板和 ggchord 坐标系。分面、`coord_flip()`、raster、boxplot 和 violin 尚不支持，显式请求会报错。纵向投影使用 `coord_collinear(seq_axis = "vertical")`。统计结果在 dodge 等数据位置调整前保留源成员。矩形和 tile 使用输入边界，不虚构零基线。

## 两条独立轨道与显示顺序 {#independent-bands}

同一 GC 既可用柱高也可用颜色编码。这里明确声明 gc_height 和 gc_colour 两个 ID；高度值域 0–1 与颜色值域 0–1 分别设置，物理带宽不同。第二幅只把颜色轨道 order 从 2 改为 0，使其靠近骨架。order 是有限数值优先级，不是测量值，编号不必从 1 开始；同优先级按 ID 排序。相同 ID 的多层需一致的 side/order/width/gap，跨独立窗口不会隐式连接折线。

{{code:comparemultitrackdata}}
{{code:comparetracktwo}}
{{figure:comparetracktwo}}
{{code:comparetrackorder}}
{{figure:comparetrackorder}}

## DNA 窗口与折线边界 {#window-policy}

此处 GC 为 `mean(bases %in% c("G", "C"))`：模糊碱基仍计入分母、对分子贡献为零；末窗口使用实际碱基数。当前 pUC19 不含模糊碱基，换输入时应明确并记录策略，不能假定自动忽略。局部聚焦裁剪区间并保留原汇总，不会重新计算缩短窗口的 GC。折线按输入样本及分组连接，`closed` 显式控制圆周闭合，`gap` 属于轨道的物理间隔。不能用闭合线暗示被省略源窗口之间存在连续观测。

## 原生源空间估计 {#estimates}

`geom_track()` 还会在**投影前**委托原生 ggplot2 执行 `smooth`、`density` 和 `ecdf`。它们是给定观测的估计，不是新增观测，按 accession/窗口及显式 group 独立计算。估计的 `source_rows` 记录原生缺失值处理后提交给估计器的组，不表示每点影响力、阳性 read support，也不是自定义模型内部进一步过滤的报告。

| 统计 | 已声明几何 | 输入和显式选择 | 数据 Position |
| --- | --- | --- | --- |
| `smooth` | `line`、`path`、`point`、`ribbon` | 源 `x` 与原数值 `y`；选择 method/formula，ribbon 需要 SE 边界 | `"identity"` 或原生 PositionIdentity |
| `density` | `line`、`path`、`point`、`area` | `x` 中的位置及可选权重；选择带宽，`after_stat(density)` 单位是源位置单位的倒数 | `"identity"` 或原生 PositionIdentity |
| `ecdf` | `line`、`path`、`point`、`step` | `x` 中的位置；`after_stat(ecdf)` 是累计比例；显式 `pad = FALSE` 避免无限源锚点 | `"identity"` 或原生 PositionIdentity |

这些统计的其他组合明确报错。density area 必须显式 `position = "identity"`，因为 area 否则默认堆叠。`bounds` 仍是 Track 源边界策略，不转发原生 density 的边界校正参数。facet、boxplot、violin、raster 及任意第三方 Geom/Stat 对象仍不支持；存在一个数值名称不等于具备投影契约。

复用上文准备的真实 `bins` 和 `track_seq`。下面显式线性拟合作为 GC 比例的显示摘要，原生置信 ribbon 与均值线共享一个 Track：

{{code:courseblockbi}}

拟合使用给定的原始源 x/y，之后 Track 数值标度的变换只映射估计，不会暗中重拟合 log 空间模型；这与原生 ggplot 在 Stat 之前放置变换有别。如需在变换空间拟合，应事先显式准备、命名变换变量并声明单位。事件位置的核密度不是测序深度，ECDF 是累计分布，不是逐位点比对。

真实 pUC19 GC 窗口及 G/C 位置已与原生 ggplot2 数值和组来源逐项对照，并核对曲线/竖向投影及显式缺失值控制。`check-track-estimates.R` 保存数值报告和实际 CI/均值图。有限源边界及物理 band 约束仍适用；自动摆放不改变模型或制造观测。

## 轨道与窗口使用相同源坐标 {#track-focus}

这里聚焦 pUC19 的 500–1500 bp，并按 interval 裁剪同一批实际 GC 窗口。focus_ggchord_track 不重新计算窗口均值，不在边界创造观测，也不会用一条线跨过省略的窗口。需要重新计算更细粒度 GC 时，应回到原 DNA 明确新窗口与分母。

{{code:comparetrackfocus}}
{{console:comparetrackfocus}}
