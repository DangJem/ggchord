<a id="coordinates"></a>

## 核对第一个窗口 {#coordinates}

这份 pUC19 bedGraph 记录每个 100 bp 窗口的 GC 百分数，末窗较短。文件第一个区间是 `[0, 100)`。读取后，`chromStart, chromEnd` 保留文件边界，供绘图的坐标是 `start = 1, end = 100`。数值 60 表示 GC 为 60%。

{{example:installed-signal}}

<a id="display"></a>

## 选择信号的几何表示 {#display}

用 `geom_track_col()` 绘制覆盖实测区间的柱时，将区间中点映射到 `x`，一基闭区间的长度映射到 `width`，测量值映射到 `y`。如果用色带，可以在 `geom_track_tile()` 中设置每格的位置和宽度，再用 `fill = value` 映射颜色。在标度上注明测量单位。[叠加柱与折线](tracks.md)有使用这份数据的完整图。

`read_bigwig()` 可从二进制信号文件读取相同的核心列。根据输入信号说明确定单位。若要计算一组区间占窗口的比例，应把原始区间交给[覆盖统计](tracks-statistics.md)。
