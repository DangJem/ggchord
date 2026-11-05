上一节调整了整张图的旋转和骨架间隙。本节仍用同样的四条序列，只改变 `MT118296.1`。每个例子都从默认设置重新开始，方便看出单个参数的作用。

## 半径：向内移动一条弧 {#shape}

{{code:frontradius}}

{{figure:frontradius}}

默认半径为 `1`。另外三条弧保留默认值，`MT118296.1` 用 `0.72`。对照[默认骨架](chord-start.html#draw)，看这条弧怎样向图心收进，序列末端仍在约 32,090 bp。半径是布局单位，读序列长度要看刻度。

## 曲率：把弧画成直线 {#curvature}

{{code:frontcurvature}}

{{figure:frontcurvature}}

`seq_curvature = 1` 是标准圆弧，`0` 把该骨架画成直线，负值会向相反方向弯。其余三条骨架保留标准圆弧。若同时调整半径和曲率，先分别观察两个参数，再把它们合用。

## 方向：让坐标从另一端开始 {#direction}

{{code:frontorientation}}

{{figure:frontorientation}}

`seq_orientation = -1` 反转该序列的显示方向。对照默认图，这条骨架的 0 刻度和箭头端互换。输入坐标和基因的生物学链仍来自数据表。三个向量都按序列 ID 命名，修改 `seq_order` 重排时仍能找到对应序列。

## 顺序：指定骨架的排列 {#order}

{{code:frontorder}}

{{figure:frontorder}}

`seq_order` 接收序列 ID 的排列。这里把 `OQ646790.1` 提到 `MT118296.1` 前面。重排后，连接表仍按原来的 ID 匹配。

接下来[改变连接带的路径与端点留白](chord-routing.html)。
