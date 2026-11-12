<div class="task-intro"><p>从真实 DNA 计算 GC 比例，再逐步叠加柱与滑动均值。读者应能说清每个轨道的单位、窗口与基线。</p></div>

<span id="prepare"></span>

## 从实际 DNA 计算数据轨道 {#domain}
本例读取随包 pUC19，每 100 bp 计算
一次 GC 比例；末尾不足 100 bp 的窗口按实际碱基数计算。

{{code:trackdata,frontgctable}}

{{output:frontgctable}}

`bins` 的每一行是一段实际序列窗口；先检查起止位置及 `gc`，再把数值映射到轨道高度与颜色。

### 绘制 GC 比例 {#draw-gc}

{{code:tracks}}
{{figure:tracks}}

值域为 0–1，单位是比例。换成其他测量值时，要声明对应单位和有效范围。
`position_track()` 决定放置与物理宽度，`scale_track_value_continuous()` 决定数值映射，
二者不能互相替代。

## 同一 ID 共享同一条轨道 {#sharing}
同名轨道上的各层共享数值域，bar 的基线也参与范围计算；它们必须声明一致的宽度和
放置方式。颜色标度按轨道 ID 独立设置。tile 表示区间类别，line 连接观测位置，
bar 表示相对基线的区间值。不要把相同数值域误认为必须位于同一条轨道。

## 检查裁剪和数值范围 {#diagnostics}
`export_ggchord_layout()` 的 `track_registry` 与 `track_diagnostics` 记录值域、
物理尺寸、原生位移及裁剪情况。被裁剪或省略的行不等于测量值为零。
局部视图使用 `focus_ggchord_track()` 并共享相同窗口。

## 叠加柱与滑动均值 {#multiple}
沿用上面的 GC 窗口，加一条三窗口滑动均值线。

{{code:multitrack}}
{{figure:multitrack}}

柱和线都来自同一组实际 GC 值，线使用三个窗口的滑动均值；两端没有完整窗口时省略，
不填造观测值。折线位置是整数窗口中点；柱与线共用同一个轨道 ID。

## 用颜色绘制 tile 轨道 {#tile}

沿用本页从真实 DNA 得到的 `track_seq` 与 `bins`。tile 要求 `accver, start, end`，并不要求数值高度；这里把同一个 `gc` 映射到 `track_fill`。

{{code:tracktile}}
{{figure:tracktile}}

与柱图使用相同 9 × 8 英寸尺寸和 100 bp 窗口。柱高表示数值；tile 的固定带宽不随 GC 改变，颜色才表示 0–1 的 GC 比例。它不是另一组实验观测，也不会自动平均重叠区间。此图比较编码方式，**不是单参数几何对照**。

| 参数或映射 | 默认值 | 本例选择与限制 |
| --- | --- | --- |
| `position` | 必须显式提供 | `position_track("gc_colour")` 声明独立轨道 |
| `width` | `position_track()` 为 4 mm | 本例 10 mm；是物理带宽 |
| `side` / `order` | outside / 1 | 多轨道按 order 排列，越小越靠近骨架 |
| `gap` | 1 mm | 到前一带的物理间隔，不是序列缺口 |
| `track_fill` | 不自动推断 | 在 `aes()` 映射 gc；固定颜色写 `fill` |
| `bounds` | `"error"` | 越界默认报错；clip 会报告裁剪 |
| `alpha` / `linewidth` | 1 / 0.25 | 固定透明度和边框线宽；默认边框 colour 为 NA |

连续颜色用 `scale_track_fill_gradientn()`；类别用 `scale_track_fill_manual()`。二者都要填写对应轨道 ID。不要把比例误标为百分数。若要叠加本页柱图，把该 tile 图层加到 `p_tracks`，保留不同 ID：这是两条独立轨道。相同 ID 必须使用一致的 side/order/width/gap。

试把 tile 的 width 从 10 改为 4 mm，其他设置不动；数值颜色不变，只有物理带宽变化。自备数据必须声明窗口、单位、NA 和重叠区间的处理。源坐标越界先查输入，不能用 clip 掩盖错误。更多见 [`geom_track_tile()`](reference/geom_track_tile.html) 和 [`position_track()`](reference/position_track.html)。

## 只改变 tile 的物理宽度 {#tile-width}

回看上面的 10 mm 色带，保持 9 × 8 英寸、输入窗口、颜色和数值域完全相同，仅把 `position_track(width)` 改为 4 mm：

{{code:comparetrackthin}}
{{figure:comparetrackthin}}

颜色仍表示同一批 GC 比例，变窄的是物理带宽。`geom_track_tile()` 负责区间色块，Position 负责带宽；它们不改变窗口里的碱基组成。把 4 改为 6 mm 练习一次，保留同一导出尺寸。柱图的数值高度另由 `track_value` 和 value scale 控制，不能用缩窄带宽代替数值归一化。

## 完成后检查 {#check-result}
bins 的 gc 在 0–1 之间，每行对应一个实际窗口；最后窗口可短于 100 bp。柱高与折线共用 0–1 的比例尺度。

<div class="next-step">继续：<a href="output.html">预览和保存图形 →</a></div>

## 保存本页结果 {#save-result}
在工作目录中保存刚刚创建的图；预览时使用相同宽高。

```r
ggsave("pUC19-tracks.png", p_multitrack, width = 9, height = 8, units = "in", dpi = 144)
```

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

```r
estimate_bins <- bins |> dplyr::mutate(pos = (start + end) / 2)
estimate_band <- position_track("gc_fit", width = unit(12, "mm"))
p_estimate <- ggchord(track_seq) +
  geom_seq() +
  geom_track(
    aes(x = pos, y = gc), data = estimate_bins, geom = "ribbon", stat = "smooth",
    method = "lm", formula = y ~ x, se = TRUE, position = "identity",
    fill = "#387D88", alpha = 0.2, track = estimate_band
  ) +
  geom_track(
    aes(x = pos, y = gc), data = estimate_bins, geom = "line", stat = "smooth",
    method = "lm", formula = y ~ x, se = FALSE, colour = "#285B65", track = estimate_band
  ) +
  scale_track_value_continuous("gc_fit", name = "GC fit", unit = "fraction") +
  coord_collinear()
p_estimate
```

拟合使用给定的原始源 x/y，之后 Track 数值标度的变换只映射估计，不会暗中重拟合 log 空间模型；这与原生 ggplot 在 Stat 之前放置变换有别。如需在变换空间拟合，应事先显式准备、命名变换变量并声明单位。事件位置的核密度不是测序深度，ECDF 是累计分布，不是逐位点比对。

真实 pUC19 GC 窗口及 G/C 位置已与原生 ggplot2 数值和组来源逐项对照，并核对曲线/竖向投影及显式缺失值控制。`check-track-estimates.R` 保存数值报告和实际 CI/均值图。有限源边界及物理 band 约束仍适用；自动摆放不改变模型或制造观测。
