Fusion 把独立构建的图放在同一画布上。跨面板带将两个原始序列区间分别连接到指定面板，总览图和局部图可以采用不同的布局。

<a id="draw"></a>

## 把圆图连接到线性图 {#draw}

中心图只含 MT108731.1，右侧面板显示两条序列。两者间最短的已观测比对把 MT108731.1 的 20,376–20,677 bp 连接到 MT118296.1 的 403–704 bp。`qpanel` 将 query 端放到圆图，`spanel` 将 subject 端放到名为 `detail` 的面板。

{{example:first-fusion}}

<span id="shared-guides"></span>

<a id="adapt"></a>

## 给每一端指定明确的落点 {#adapt}

先核对两端都在指定面板的可见区间内，并使用相容的单位。`fusion_panel()` 的 side 决定子图摆放位置，子图自身的坐标系控制骨架方向。要把总览连接到来源窗口，可用局部放大函数创建并登记该窗口。

要集中放置外侧图例，可将 `guides = "collect"` 传给 `ggchord_fusion()`。当不同面板表示同一度量时，给它们设置相同的范围、刻度和样式。组合收集图例，各面板保留各自的标度。

[六面板比较图](case-fusion.md) · [跨面板连接](../../man/fusion_ribbon.Rd) · [来源窗口放大](../../man/magnify_ggchord.Rd)
