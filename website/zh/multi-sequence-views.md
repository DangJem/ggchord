<a id="comparison"></a>

## 固定输入，比较布局 {#comparison}

四个示例均使用两条完整参考，以及 MT108731.1 对 MT118296.1 的三条实测 HSP。每段代码都能独立运行。可比较每条连接带与骨架相接的位置。四种视图均保留两条参考的全长和 HSP 的端点顺序。

<a id="chord"></a>

### 弦图 {#chord}

沿圆周读取连接带两端。多对多关系较多时，可用这种布局查看整体结构。

{{example:views-chord}}

<a id="collinear"></a>

### 共线性图 {#collinear}

沿各行比较三条比对的位置和跨度。要翻转其中一条参考的显示方向，见[序列排列与显示方向](multi-sequence-layout.md)。

{{example:views-collinear}}

<a id="radial"></a>

### 径向共线性图 {#radial}

从中心向外沿各条辐条读取位置。加入更多参考时，辐条顺序决定哪些参考彼此相邻。

{{example:views-radial}}

<a id="sankey"></a>

### Sankey 双列图 {#sankey}

从左侧的 query 副本读向右侧的 subject 副本。带在两端的宽度分别反映该条序列上的观测跨度。需要区分 query 和 subject 角色时，可选择这种双列布局。

{{example:views-sankey}}

需要带有数据与绘图脚本的完整多序列案例时，可查阅[案例图库](gallery.md)。
