<a id="comparison"></a>

## 在四种坐标系中追踪同一条比对 {#comparison}

四幅图都使用 MT108731.1（64,323 bp）和 MT118296.1（32,090 bp）的全长，以及从前者到后者的第一条局部比对。每段示例均可独立运行，各图层使用相同的参考及比对区间，由各视图选择自己的间距和坐标系。

<a id="chord"></a>

### 在两段圆弧上找到端点区间 {#chord}

{{example:views-chord}}

较长的蓝色圆弧是 MT108731.1。带子从这条弧的 26,298–51,139 bp 连到较短橙色弧的 7,121–31,959 bp。骨架上的箭头标出显示方向。[逐层画出这幅弦图](chord.md)。

<a id="collinear"></a>

### 沿两条直线轴读取区间 {#collinear}

{{example:views-collinear}}

较短的序列现在位于较长序列下方。两行的左端都是第 1 bp。连接带的下端几乎延伸到该行骨架的右端。[加入另外两条比对](collinear.md#adapt)。

<a id="radial"></a>

### 从中心向外沿辐条读数 {#radial}

{{example:views-radial}}

现在两根辐条的第 1 bp 都朝向中央空隙。沿带子在上下两根辐条之间移动，可以找到相同的两个区间。[加入第三根辐条](radial-synteny.md#adapt)，可观察两条比对如何共用一条查询序列。

<a id="sankey"></a>

### 从左侧查询序列读到右侧目标序列 {#sankey}

{{example:views-sankey}}

每个登录号在左右两列各显示一次。唯一的输入记录从左上方的骨架连向右下方的骨架。下一步可[比较这对序列的三条连接带](sankey.md#adapt)，或[在点阵图的两条轴上读取比对坐标](dotplot.md)。
