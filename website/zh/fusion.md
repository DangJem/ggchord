Fusion 可组合各自拥有独立坐标系的图。这里用圆图显示 MT108731.1 全长，在右侧用两条竖直骨架并列显示它和 MT118296.1。随后用一条实测 HSP 把圆上的区间连到第二条骨架。

<a id="draw"></a>

## 组合两幅独立构建的图 {#draw}

先分别构建圆形总览图和线性对照图，再把具名列表 `list(overview = overview, detail = detail)` 传给 `ggchord_fusion()`。默认网格按列表顺序排列两幅图。连接图层用这些名称指定面板路径。

{{example:fusion-panels}}

<a id="adapt"></a>

## 按面板标识连接一对区间 {#adapt}

继续使用这个查询方向的第一条 HSP。在传给 `fusion_ribbon()` 的关系表中加入 `qpanel_path` 和 `spanel_path` 两列。查询端连到 `overview` 中 MT108731.1 的 26,298–51,139 bp，目标端连到 `detail` 中 MT118296.1 的 7,121–31,959 bp。序列标识确定来源，面板路径确定显示位置。

{{example:first-fusion}}

<span id="shared-guides"></span>

接收端点的面板必须包含相应的参考序列及区间。要汇集各面板中相同的外置图例，可在 `ggchord_fusion()` 中设置 `guides = "collect"`，并为需要比较的测量值指定相同的标度范围和刻度。下一步可沿[六面板案例](case-fusion.md)从基因组总览读到局部比对，或[按最终宽度导出图形](output.md)。
