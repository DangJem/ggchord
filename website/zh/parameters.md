灵活的参数格式让支持它的参数从统一设置细化到局部控制。下面始终使用相同四条随包序列；前三幅图只比较 `seq_radius`，画布均为 7 × 7 英寸。需要 ggchord；真实 CDS 整理使用可选 dplyr。

{{code:courseblockac}}

## 一个值应用到所有序列 {#scalar}

{{code:parameteruniform}}
{{figure:parameteruniform}}

半径统一为 1。所有序列采用相同设置时，单值已经足够。

## 按名称只改一条序列 {#recommendation}

{{code:parameters}}
{{figure:parameterexample}}

只把 `MT118296.1` 的半径改为 0.72，其余保留本参数默认值 1。键是实际 `accver`；调整显示顺序后仍按 ID 匹配。完整配置中每条序列确有不同值时，保留完整命名向量。

## 按显示顺序分别设置 {#vector}

<span id="compatibility"></span>

{{code:parametervector}}
{{figure:parametervector}}

四个值对应处理 `seq_order` 后的显示顺序，长度须与显示序列数相同。当前第二条是 `MT118296.1`，所以本图与上一张相同。

### 重排时两种匹配的差异 {#reorder}

保持原向量，只改变 `seq_order`。命名值仍修改 `MT118296.1`；无名向量的第二个值改为作用在第二条显示序列 `MT108731.1`。

{{code:parameterreorder}}

```r
view_ggchord(p_parameter_named_order, width = 7, height = 7, units = "in")
view_ggchord(p_parameter_ordered_values, width = 7, height = 7, units = "in")
```

![重排后的命名值](../assets/parameternamedorder.png)
![重排后的顺序值](../assets/parameterorderedvalues.png)

## 准备真实正负链 CDS {#genes}

{{code:parametergenes}}
{{output:parametergenes}}

`OR222515.1` 同时有 2624–4195 的正链和 19246–20523 的负链记录；`OQ646790.1` 的 integrase 为正链。所有链取自 `source_strand`，不改变源坐标。

## 分别设置链的位置 {#strand}

{{code:parameterstrand}}
{{figure:parameterstrand}}

显式链值已经带符号：正链 0.08、负链 −0.08，不再次反转。单值 `position_strand(0.08)` 则按规则展开为正链 0.08、负链 −0.08。

### 同时按序列和链设置 {#lists}

实际目标是让 `OR222515.1` 正负链不对称偏移，另一条序列采用较小的对称偏移。此时列表能表达两层差异。

{{code:parameterlist}}

![按真实序列和链设置偏移](../assets/parameterlist.png)

`OR222515.1` 使用显式的 0.08/−0.14；`OQ646790.1` 的单值 0.06 按正负链对称展开。未指定序列或链在这个 Position 家族中补 0。该混合列表不能直接复制给 `width`。

## 分别设置链的宽度 {#width}

{{code:parameterwidth}}
{{figure:parameterwidth}}

正链宽度 0.03、负链 0.08，改变箭头厚度，区间长度仍来自 start/end。需要逐序列宽度列表时，各元素须满足本家族的链值规则，详见[几何参数](reference/geom_feature.html)。

## 按序列放置名称 {#labels}

{{code:parameterlabels}}
{{figure:parameterlabels}}

`MT108731.1` 与 `OQ646790.1` 的文字半径分别为骨架的 75% 和 115%，其余使用默认位置。显式 `seq_label_radius` 不适用于 Collinear。

## 参数家族与未指定项 {#scope}

| 接口 | 可用格式与匹配 | 边界与未指定项 |
| --- | --- | --- |
| `geom_seq()` | 单值、按 accver、按显示顺序及相应标量列表 | `seq_radius`、`seq_curvature`、`seq_gap`、`seq_orientation` 各有默认；效果与支持范围依 Coord |
| `geom_seq_label()` | `labels`、半径、旋转、size、hjust/vjust 按序列 | `seq_label_orientation` 是单个模式；显式半径与 Collinear 不兼容 |
| `geom_link_ribbon()` / `geom_link_line()` | `link_gap` 按端点所在序列 | `NULL` 保留自动行为；控制点不属于逐序列标量 |
| `geom_gene()` / `geom_feature()` | `width` 按序列或链，列表元素遵循链值规则 | 未指定项用本参数默认；不接受 Position 的任意混合列表 |
| `position_strand()` | `offset` 可单值、按序列、显式链值、逐序列混合列表 | 单值 x 展开为 +x/−x；显式链值不反转；未指定序列或链补 0 |
| `position_plasmid()` | `offset` 单值或逐序列标量配置 | 两链共用偏移，不支持逐链；未指定序列补 0 |
| 固定 `geom_gene_label()` / `geom_feature_label()` | `label_radial_offset`、`label_circum_offset`、`label_circum_limit`、`label_rotation` 按序列或链 | 不代表所有文本参数支持列表；按本家族链值规则解析 |
| `position_track()` | 单个绝对 `unit()`；或原生数值的逐序列 width/gap | 数值须覆盖全部显示序列；物理尺寸保留毫米含义 |

`seq_order` 是排列，`label_layout` 是布局模式。自动避让层不接收任意固定偏移。axis/theme、aesthetic、标度与 guide 不继承这套列表规则；手工颜色标度按**数据类别**匹配，见[配色说明](styling.html)。精确默认值和取值范围见[参数图谱](parameter-atlas.html)。

## 选择写法与核对结果 {#check-result}

统一设置用单值；少数局部差异优先用命名值；依赖显示顺序时用等长向量；两层差异确实必要且接口支持时才用列表。核对 ID、显示顺序、未指定项和真实结果。其他技巧页采用适合任务的最短表达，不按向量或列表的数量评价示例。
