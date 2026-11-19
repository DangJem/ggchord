四条来源序列都在左右两列出现。即使两个 accession 在行序中相距较远，也能从左列 query 端追踪到右列 subject 端。

<a id="map"></a>

## 绘制完整图 {#map}

该图使用 10 × 16 英寸画布、144 dpi。脚本返回绘图对象，可传给 `ggplot2::ggsave()` 按相同尺寸保存。

{{example:case-sankey}}

<a id="reading"></a>

## 读图 {#reading}

输入含十条精选 HSP 和 16 条 CDS。带的填色表示一致性百分比，箭头标出 CDS 的链方向，序列配色帮助识别两列中同一 accession 的副本。另一列中的同一序列为这些特征提供了第二个显示位置。

<a id="inputs"></a>

## 来源记录与展示选择 {#inputs}

相邻行和跨行连接均来自共线性图的展示子集。两列的间距设为 `column_gap = 4.2`，序列与链图例共用。若缩小导出宽度，先检查两列外侧的标签空间。

[双列区间图读法](sankey.md) · [对照线性布局](case-collinear.md)
