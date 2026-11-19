Fusion 可组合各自拥有独立坐标系的图。这里用圆图显示 MT108731.1 全长，在右侧用两条竖直骨架并列显示它和 MT118296.1。随后用一条实测 HSP 把圆上的区间连到第二条骨架。

<a id="draw"></a>

## 组合两幅独立构建的图 {#draw}

先分别构建圆形总览图和线性对照图。将总览图传给 `ggchord_fusion()`，再用 `fusion_panel()` 把第二幅图放在右侧。中央面板的标识为 `center`，给新面板取名 `detail`，后续图层就可以按这个标识指定位置。

{{example:fusion-panels}}

<a id="adapt"></a>

## 按面板标识连接一对区间 {#adapt}

继续使用这个查询方向的第一条 HSP。在 `fusion_ribbon()` 中，`qpanel` 将查询端放到圆图上，对应 MT108731.1 的 26,298–51,139 bp。`spanel` 将目标端放到右侧面板，对应 MT118296.1 的 7,121–31,959 bp。MT108731.1 也出现在右侧面板，但这条带的查询端连在圆图上。

{{example:first-fusion}}

<span id="shared-guides"></span>

接收端点的面板必须包含相应的参考序列及区间。要汇集各面板中相同的外置图例，可在 `ggchord_fusion()` 中设置 `guides = "collect"`，并为需要比较的测量值指定相同的标度范围和刻度。下一步可沿[六面板案例](case-fusion.md)从基因组总览读到局部比对，或[按最终宽度导出图形](output.md)。
