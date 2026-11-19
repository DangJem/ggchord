四条组装序列长度各异，部分区间互相匹配。图中把直线和弧形骨架排在不同位置，以便同时读取 CDS 标签和跨序列连接。

<a id="map"></a>

## 绘制完整图 {#map}

该图使用 12 × 10 英寸画布、144 dpi。脚本返回绘图对象，可传给 `ggplot2::ggsave()` 按相同尺寸保存。

{{example:case-chord}}

<a id="reading"></a>

## 读图 {#reading}

从序列标签开始，按骨架刻度读取 bp 位置。每条带的两端包围已观测的比对区间，填色表示一致性百分比。骨架两侧的箭头标出具有来源链记录的 CDS。

<a id="inputs"></a>

## 来源记录与展示选择 {#inputs}

`flagship_comparison_example$chord` 保存四条全长、12 条 HSP 和 20 条入图 CDS。先保留长度至少 900 bp 且一致性至少 75% 的比对，再从每对序列取最长的至多三条。`layout_sequences()` 和 `geom_seq()` 的具名设置控制各序列的显示方向、半径与曲率。

[从一条比对带开始](chord.md) · [在径向辐条上查看同一输入](case-radial.md)
