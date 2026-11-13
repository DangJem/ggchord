{{code:courseblockaf}}

{{prepare:frontbackbone}}

## 方向：让坐标从另一端开始 {#direction}

{{code:frontorientation}}

{{figure:frontorientation}}

`seq_orientation = -1` 反转该序列的显示方向。对照默认图，这条骨架的 0 刻度和箭头端互换。输入坐标和基因的生物学链仍来自数据表。三个向量都按序列 ID 命名，修改 `seq_order` 重排时仍能找到对应序列。

## 顺序：指定骨架的排列 {#order}

{{code:frontorder}}

{{figure:frontorder}}

`seq_order` 接收序列 ID 的排列。这里把 `OQ646790.1` 提到 `MT118296.1` 前面。重排后，连接表仍按原来的 ID 匹配。

接下来[改变连接带的路径与端点留白](multi-sequence-links.html#control)。

## 线性锚点 {#anchors}

[Collinear](collinear-layout.html#anchor) · [Radial](radial-synteny.html) · [Sankey](sankey.html)
