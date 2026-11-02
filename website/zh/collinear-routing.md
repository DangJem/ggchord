相邻行连接与跨行连接有不同阅读负担。先分开数据，再分别控制；保留原始HSP端点，控制点只改变显示路径。

## 本页准备 {#prepare}
先解压[教学输入](../assets/tutorial-inputs.zip)。
{{prepare:sequencepalette,docstyle,showcasedata}}

## 明确长连接的控制点 {#routes}
{{code:showcasecollinear}}
{{figure:showcasecollinear}}
线性宣传图省略tail assembly protein与gp4 protein两条拥挤候选，它们仍保留在原始表和弦图中；其余16个特征逐一标注。相邻行保留筛选HSP，非相邻每对只选最长一条。`show_cross_controls` 与 `show_cross` 的行顺序一一对应；筛选或重新排序后应一同重建。控制点单位是导出布局x/y，不是bp。当前接口支持每条连接2个或4个数值；它塑造两条ribbon边的控制位置，并非全局路径规划器。

本例将短序列水平居中，使用骨架预布局导出的真实行位置、长度与偏移计算控制点，放在中间序列右端附近；不再把路径拉到远处。行距与seq_gap及最长骨架宽度有关，不能把bp直接填进控制点。修改行距、顺序或方向后，应重新检查实际路径；不能保证任意控制点都无自交。

## 局部缺口与虚线 {#clearance}
上下分链用 `position_strand()`；`link_avoid = "smooth"` 配合默认NULL gap，在箭头附近产生局部退让。跨行外框默认虚线；本例增大linewidth使其可见，不用显式linetype覆盖。此规则依据最终显示邻接，不表示反向比对。

## 比例尺为何曾经没对齐 {#ruler}
`geom_seq_scale()` 的默认offset是横向2 mm、纵向6 mm，因此原图左端略向右偏移是明确配置，而非bp长度计算出错。本例使用 `offset = unit(c(0, 10), "mm")`，水平起点对齐左侧序列包络，纵向留出标签空间。本例采用水平居中，各行源坐标0不再相互对齐；比例尺对齐左侧包络，反向行的源坐标0还会移到另一端。比例尺衡量距离，不代表某个源区间。

## 继续深入 {#next}
局部比较优先用[focus/分段](views.html)，同源证据见[基因邻域](homology.html)，最终结果见[2.6 综合图](collinear.html)。GC等定量轨道留在[轨道专题](tracks.html)，首页为标签和连接留空间。
