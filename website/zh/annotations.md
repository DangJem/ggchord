<div class="task-intro"><p>把真实注释画为独立图层，选择合适形状，并让重叠区间和文字共享堆叠位置。</p></div>

<span id="prepare"></span>

## 保留注释语义 {#stack}
用 pUC19 识别结果的前六条区间绘制圆图。
这里保留 pUC19 中实际重叠的注释和原本未知的链方向。
两层使用同一个堆叠位置，文字跟随对应注释。

{{code:featurestack}}
{{figure:featurestack}}

`position_feature_stack()` 解决注释轨道重叠。标签层使用同一 Position，让文字锚点
跟随轨道。标签子集或重排时优先提供唯一 `feature_id`。严格的 gene 快捷层要求
已知 `+/-`，通用注释保留未知/双向语义。

## 按含义选择形状 {#shapes}
| 几何 | 常见用途 | 重要区别 |
| --- | --- | --- |
| `arrow`, `compact_arrow` | 有方向的编码/注释区间 | 已知方向来自证据，不由显示猜测 |
| `block`, `capped_line` | 区域或不填充区间 | 选择形状不推断生物学功能 |
| `promoter_arrow`, `primer_arrow`, `primer_arc` | 已声明的启动子/引物式注释 | 计算匹配与绘图分开 |
| `marker`, `chevron`, `lollipop` | 自定义区间/类点强调 | 保留真实区间元数据，显示形状属于映射 |

`geom_gene()`、`geom_feature_plasmid()` 是有用预设，默认值并不相同。
`position_strand()` 提供带符号位移，`position_plasmid()` 提供圆图放置预设，
后者是局部 Position，不新增 Coord。

## 高亮区间而不改写观测 {#decoration}
{{code:decorations}}
{{figure:decorations}}

图中高亮 pUC19 复制起点的实际范围。区域叠层有自己的数据和 `region_fill` 标度。`geom_seq_label()` 命名序列，
`geom_seq_center_label()` 在单圆图中心提供名称/长度，`geom_seq_break()` 标记裁剪或
分段末端。`geom_seq_scale()` 是仅用于共线性图的独立序列位置参照条，与控制序列基因组轴
的 `scale_seq_position_continuous()` 不同。普通 ggplot2 注释应使用
[位置映射](views.html#mapping)，不能把 bp 直接当作笛卡尔 x/y。

## 完成后检查 {#check-result}
先检查区间与链，再看形状与标签。未知链仍保留在通用注释中；堆叠改变显示位置，不改变原始基因组坐标。

<div class="next-step">继续：<a href="labels.html">放置和调整标签 →</a></div>

比较真实 CDS 的 shouldered/flush/triangle 三种箭头，见[箭头几何对照](feature-geometry.html#arrow-heads)。

## 区间高亮的侧向、宽度与偏移 {#region-parameters}

从第一条真实 HSP 的 query 端保留原始 ID 与区间。默认宽度 0.08、偏移 0、侧向 inside；后图分别只改变 outside、宽度或偏移。它们是局部几何量，不是 bp。固定填色写在 geom 外观参数，数据类别映射使用 region_fill 和其专用标度。不能把高亮当作新的功能注释。

{{code:comparefinishchord}}
{{code:compareregioninside}}
{{figure:compareregioninside}}
{{code:compareregionoutside}}
{{figure:compareregionoutside}}
{{code:compareregionwidth}}
{{figure:compareregionwidth}}
{{code:compareregionoffset}}
{{figure:compareregionoffset}}
