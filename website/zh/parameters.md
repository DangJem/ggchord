灵活参数可以在同一图层中为不同序列或链指定不同的几何值。先移动一个标签，再用同样的匹配方式调整骨架半径、重排序列，最后按链指定厚度或偏移。函数参考中参数旁的**灵活参数**链接可回到本页。

<a id="named-label"></a>

## 只移动一条噬菌体序列的标签 {#named-label}

保留两条骨架的半径 1。将 MT108731.1 的标签半径设为 0.72，将 MT118296.1 的设为 1。每个名称都匹配一条弧上的登录号。较长的蓝色圆弧位置不变，它的标签移向圆心。

{{example:flex-labels}}

<a id="scalar"></a>

## 改变一条序列的骨架半径 {#scalar}

<span id="recommendation"></span>

随包数据 `seq_data_example` 中的四条噬菌体基因组序列默认半径均为 1。给 `MT118296.1` 设置半径 0.72 即可把这条骨架向内移。名称按序列的 `accver` 匹配，其余三条仍为 1。

{{example:flex-radius}}

| `seq_radius` 的值 | 哪些骨架会改变？ |
| --- | --- |
| `0.72` | 全部四条 |
| `c("MT118296.1" = 0.72)` | 仅 MT118296.1 |
| `c(1, 0.72, 1, 1)` | 第二条显示序列 |
| `list("MT118296.1" = 0.72)` | 仅 MT118296.1 |

<a id="reorder"></a>

## 按名称还是按显示顺序匹配 {#reorder}

<span id="vector"></span>
<span id="compatibility"></span>

左右两图使用相同的四条序列和相同的重排次序。左图的小半径跟随名称 `MT118296.1`。右图的小半径跟随**第二个显示位置**，此处已换为 `MT108731.1`。两图的序列长度、颜色和序列标签均相同。

{{example:flex-radius-order}}

无名向量的长度须与显示序列数一致。用列表指定半径时，每个元素都是一个数。只改少数序列时，命名向量写法最短，重排后也能保留各条序列的设置。

<a id="strand"></a>

## 设置链偏移和箭头宽度 {#strand}

<span id="genes"></span>
<span id="lists"></span>
<span id="width"></span>

序列 `OR222515.1` 上的两条 CDS 来源链相反。给正链 CDS 设向外偏移 0.08，给负链 CDS 设向内偏移 −0.14。给 `OQ646790.1` 一个值 0.06，则它的正负链对称偏移。示例用 dplyr 筛选 CDS，并从记录的 `source_strand` 列取链向。

{{example:flex-strand}}

显式负链偏移自带符号。`position_strand()` 将未指定的序列或链的偏移补为 0。宽度改变箭头厚度。这里正链宽 0.03，负链宽 0.08，箭头端点仍取自 `start` 和 `end`。按序列**再按链**指定宽度时，每个命名元素都应提供一个链向量，例如 `list("OR222515.1" = c("+" = 0.03, "-" = 0.08))`。

<a id="invalid-radius"></a>

## 用链名指定半径会发生什么？ {#invalid-radius}

每条骨架对应一条序列。把 `+` 和 `-` 写成半径的名称，相当于请求两个输入中没有的序列标识。下面的示例捕获布局计算时的报错。要调整骨架半径，用 `accver` 作名称；要调整正负链偏移，使用 [position_strand()](../../man/position_strand.Rd)。

{{example:flex-radius-error}}

<a id="scope"></a>

## 按绘图任务查参数 {#scope}

<span id="labels"></span>
<span id="check-result"></span>

| 用途 | 参数入口 | 匹配规则 |
| --- | --- | --- |
| 调整序列骨架 | [geom_seq()](../../man/geom_seq.Rd) · `seq_radius`、`seq_curvature` | 两者默认均为 1。半径取有限正数，曲率取有限的有符号数。 |
| 设置方向与间隙 | [layout_sequences()](../../man/layout_sequences.Rd) · `orientation`、`gap` | 每条序列的方向取 1 或 −1。间隙的单位和默认值依视图而定。 |
| 放置序列名称 | [geom_seq_label()](../../man/geom_seq_label.Rd) · `labels`、`seq_label_radius`、`seq_label_rotation`、`seq_label_hjust`、`seq_label_vjust` | 名称按序列标识匹配。标签半径取相对该骨架的比例。 |
| 调整连接端点间隙 | [geom_link_ribbon()](../../man/geom_link_ribbon.Rd) 或 [geom_link_line()](../../man/geom_link_line.Rd) · `link_gap` | 每个端点按其所在序列取间隙值。 |
| 设置特征厚度 | [geom_feature()](../../man/geom_feature.Rd) · `width` | 可按序列、按链或同时按序列与链设置。 |
| 移动特征 | [position_strand()](../../man/position_strand.Rd) 或 [position_offset()](../../man/position_offset.Rd) · `offset` | 两条链分别移动用前者，每条序列的两条链共用偏移用后者。 |
| 固定标签位置 | [geom_feature_label()](../../man/geom_feature_label.Rd) 或 [geom_gene_label()](../../man/geom_gene_label.Rd) · `label_radial_offset`、`label_circum_offset`、`label_circum_limit`、`label_rotation` | 固定布局可按序列或链设置。自动布局自行决定标签位置。 |
| 调整定量轨道尺寸 | [track_band()](../../man/track_band.Rd) · `width`、`gap` | 取一个物理单位值，或用原生布局数值覆盖全部显示序列。宽度须为正，间隙可为零。 |

基于数据的颜色用[视觉映射与标度](styling.md)设置。函数签名、默认值及灵活参数标记见[API 参考](reference.md)。
