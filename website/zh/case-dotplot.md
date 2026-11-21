四条全长参考序列组成 4 × 4 个有向序列对。图中完整保留随包比较数据中的 32 条 HSP。

<a id="map"></a>

## 绘制完整图 {#map}

该图使用 13 × 12 英寸画布、144 dpi。脚本返回绘图对象，可传给 `ggsave()` 按相同尺寸保存。

{{example:case-dotplot}}

<a id="reading"></a>

## 读图 {#reading}

从列标题读取 query accession，从行标题读取 subject accession。在对应格内，线色表示一致性，线宽表示以 bp 计的比对长度，斜率由成对端点的顺序决定。格内文字给出该有向序列对的输入 HSP 条数。

<a id="inputs"></a>

## 来源记录与展示选择 {#inputs}

六个无序序列组合都只在一个 query→subject 方向有记录。数据中没有自比对行或互换两端的副本。`flagship_dotplot_example$cell_labels` 在当前 4,000 bp 的轴间隙下提供全部 16 个单元格标签；若改变间隙或轴顺序，应重新确定标签位置。

[绘制一条 HSP](dotplot.md) · [点图坐标轴](../../man/coord_dotplot.Rd)
