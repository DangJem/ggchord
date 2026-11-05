共线性图默认只画相邻行的连接。用一条真实的第一行到第三行的比对，观察怎样显式显示跨行带。

## 从相邻连接开始 {#prepare}

{{code:frontlinearbase,frontlinks,frontlinearlink}}

{{figure:frontlinearlink}}

## 允许跨行连接 {#routes}

{{code:frontlinearcross}}

{{figure:frontlinearcross}}

相邻带默认使用细实线，非相邻行之间的带默认使用较粗虚线。普通线宽由 `linewidth` 控制；跨行线宽由 `nonadjacent_linewidth` 单独控制，默认 0.4。

第二个连接层只接收一条第一行到第三行的比对。`adjacent_only = FALSE` 让这条带可见，默认虚线轮廓提示它跨过中间一行。原来第一、二行的连接仍由第一个层绘制。

## 单独设置跨行边框 {#width}

{{code:frontcrosswidth}}
{{figure:frontcrosswidth}}

`nonadjacent_linewidth = 0.8` 加粗第二个图层中的跨行带，不改变第一个图层的相邻带，也不改变线型。

{{code:frontcrosssame}}
{{figure:frontcrosssame}}

设为 0 时，跨行带使用该图层的普通 `linewidth`，这里是 0.2。它仍是虚线。若使用 `aes(linewidth = ...)` 映射，实际线宽由映射及标度控制。

## 控制点与避障 {#clearance}

如果跨行带遮住了中间行，可以用 `link_ctrl_point` 将带引向侧边。控制点使用布局 x/y 坐标，连接端点仍使用 bp。多条跨行带的控制点，以及 `link_avoid = "smooth"` 的组合见[四序列综合案例](collinear.html#compose)。

## 比例尺与坐标零点 {#ruler}
默认标尺横向偏移为 0 mm，纵向留 6 mm。它对齐左侧序列包络。左对齐正向骨架的左端就是零点，居中或翻转后须从各行刻度读取零点。

## 继续深入 {#next}
局部比较优先用[focus/分段](views.html)，同源证据见[基因邻域](homology.html)，最终结果见[综合图](collinear.html)。GC等定量轨道留在[轨道专题](tracks.html)，首页为标签和连接留空间。
