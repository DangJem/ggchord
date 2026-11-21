灵活参数允许对不同的参考序列或显示实例分别设置。推荐用具名向量按序列标识匹配。在指定物理位置时，数值以毫米为单位，也可按参数说明提供绝对 grid 单位。API 参考中的「【支持灵活参数】」标记会链接到下面的对应小节。

<a id="named-label"></a>

## 移动一条噬菌体序列的标签 {#named-label}

两条骨架半径保持为 1。固定标签采用水平方向，为 MT108731.1 指定 18 mm 的具名偏移，再用 `side = "negative"` 指定向内侧。只有具名的标签向内移动。

{{example:flex-labels}}

<a id="scalar"></a>

## 改变一条骨架的半径 {#scalar}

<span id="recommendation"></span>

四条序列的默认半径均为 1。在 `layout_sequences()` 中设置 `radius`，即可让 MT118296.1 向内移动。半径取相对标准半径的比例。

{{example:flex-radius}}

| `radius` 的取值 | 受影响的参考序列 |
| --- | --- |
| `0.72` | 全部四条 |
| `c("MT118296.1" = 0.72)` | 仅 MT118296.1 |
| `c(1, 0.72, 1, 1)` | 第二条显示序列 |
| `list(source = c("MT118296.1" = 0.72))` | 仅 MT118296.1 |

<a id="reorder"></a>

## 按名称或显示顺序匹配 {#reorder}

<span id="vector"></span>
<span id="compatibility"></span>

两个面板使用同一个重排后的布局。左图的小半径按序列标识匹配到 MT118296.1。右图的小半径按第二个显示位置匹配到 MT108731.1。无名向量须为每个显示实例提供一个值。要让实例设置覆盖共用的来源设置，可以用 `list(source = ..., instance = ...)`。

{{example:flex-radius-order}}

<a id="strand"></a>

## 设置链偏移和箭头厚度 {#strand}

<span id="genes"></span>
<span id="lists"></span>
<span id="width"></span>

这幅弦图中，OR222515.1 的正链和负链 CDS 分别偏移 5 mm 和 −9 mm。对 OQ646790.1 指定单值 4 mm，两条链将等距反向偏移。示例读取原始链向，用两个图层分别给正链和负链设置 2 mm 和 5 mm 的箭头厚度。

{{example:flex-strand}}

<a id="invalid-radius"></a>

## 半径以链符号命名会怎样？ {#invalid-radius}

骨架属于参考序列。用 `+` 和 `-` 作半径名称时，布局会查找两个输入中没有的序列标识。示例捕获此报错。按链移动特征时，使用 [position_strand()](../../man/position_strand.Rd)。

{{example:flex-radius-error}}

<a id="scope"></a>

## 按任务查找参数 {#scope}

<span id="labels"></span>
<span id="check-result"></span>

| 任务 | 可查阅的参数 | 取值规则 |
| --- | --- | --- |
| 调整弧形参考 | [layout_sequences()](../../man/layout_sequences.Rd) · `radius`、`curvature`、`angular_gap` | 半径取正比例值，曲率可带符号，角间距以度计。 |
| 调整直线型参考 | [layout_sequences()](../../man/layout_sequences.Rd) · `orientation`、`row_gap`、`column_gap` | 方向取 1 或 −1。共线性图和 Sankey 图的行距、Sankey 图的列距以毫米计。 |
| 放置固定标签 | [label_layout_fixed()](../../man/label_layout_fixed.Rd) · `offset`、`along` | 垂直于参考的距离和沿参考的位移均以毫米计。在图层支持时可用具名值。 |
| 移动特征或分道 | [position_strand()](../../man/position_strand.Rd)、[position_offset()](../../man/position_offset.Rd) 和 [position_feature_stack()](../../man/position_feature_stack.Rd) | 按链偏移、共用偏移或指定道间距。 |
| 设置特征厚度 | [geom_feature()](../../man/geom_feature.Rd) 和 [geom_gene()](../../man/geom_gene.Rd) · `thickness` | 取以毫米计的正值。可按参考命名，也可显式指定实例。 |
| 调整定量轨带尺寸 | [track_band()](../../man/track_band.Rd) · `thickness`、`gap` | 厚度取正值，间隙可为零。两者均为物理尺寸。 |

基于数据的颜色用[视觉映射与标度](styling.md)设置。[API 参考](reference.md)列出各函数的参数和默认值。
