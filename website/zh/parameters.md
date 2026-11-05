需要所有序列共用一个值、分别按序列指定，或按正负链设置时，可以沿用同一套写法。下面各图都使用四条随包序列。

## 一个值应用到所有序列 {#scalar}

{{code:parameteruniform}}
{{figure:parameteruniform}}

四个骨架间隙都设为 0.08。单个数适用于所有序列。

## 按名称只改一条序列 {#recommendation}

{{code:parameters}}
{{figure:parameterexample}}

命名向量把 MT118296.1 的半径设为 0.72，其余保留默认值 1。重排序列后，值仍按 ID 匹配。

## 按显示顺序分别设置 {#vector}

<span id="compatibility"></span>

{{code:parametervector}}
{{figure:parametervector}}

四个半径依次对应序列表中的四条记录。无名向量必须与序列数等长；若经常调整顺序，使用命名向量更容易核对。

## 为正负链准备两条真实 CDS {#genes}

{{code:parametergenes}}
{{output:parametergenes}}

integrase 是 OQ646790.1 上的正链 CDS，Int 是 OR222515.1 上的负链 CDS。保留原始链来演示下面两种设置。

## 分别设置链的位置 {#strand}

{{code:parameterstrand}}
{{figure:parameterstrand}}

正链偏移为 0.08，负链为 −0.08。位置由 `position_strand()` 控制，表中的坐标与链不变。

## 分别设置链的宽度 {#width}

{{code:parameterwidth}}
{{figure:parameterwidth}}

正链箭头宽度为 0.03，负链为 0.08。箭头宽度改变显示厚度，基因区间长度仍由 start/end 决定。

## 按序列放置名称 {#labels}

{{code:parameterlabels}}
{{figure:parameterlabels}}

MT108731.1 的文字半径为骨架半径的 75%，OQ646790.1 为 115%；另外两条使用默认位置。每个支持的参数独立接收命名值。

## 哪些参数支持 {#scope}
| 接口 | 支持参数 | 匹配范围 |
| --- | --- | --- |
| `geom_seq()` | `seq_radius`, `seq_curvature`, `seq_gap`, `seq_orientation` | 按序列；半径/曲率的效果取决于坐标模式 |
| `geom_seq_label()` | `labels`, `seq_label_radius`, `seq_label_rotation`, `size`, `seq_label_hjust`, `seq_label_vjust` | 按序列；`seq_label_orientation` 是单个模式选项 |
| `geom_link_ribbon()` / `geom_link_line()` | `link_gap` | 按端点所在序列；控制点不是逐序列标量 |
| `geom_gene()` / `geom_feature()` | `width` | 按序列，或分别按 `+`/`-` 链 |
| `position_strand()` / `position_plasmid()` | `offset` | 按序列或链；优先使用 Position 而非兼容 offset 入口 |
| 固定标签 `geom_gene_label()` / `geom_feature_label()` | `label_radial_offset`, `label_circum_offset`, `label_circum_limit`, `label_rotation` | 按序列或链；不代表所有文本参数都接受列表 |
| `position_track()` | 数值型 `width`, `gap` | 原生布局数值可按序列；物理 `unit()` 使用单个绝对尺寸 |

`seq_order` 本身就是序列排列，不是上述值的循环规则。自动避让标签的 `label_layout`
是布局模式；手动角度和偏移请使用固定标签层，不把它们塞进自动标签的 `...`。
axis/theme、aesthetic、标度、图例不继承这套列表扩展规则。手工标度本来就支持
按**数据类别**命名的颜色向量，那属于 ggplot2 的标度语义，见[配色教程](styling.html)。

## 完成后检查 {#check-result}
命名向量的键应对应实际 accver 或 +/− 链。重排后再次检查匹配；不要将列表规则推广到 theme、scale 或 aesthetic。

<div class="next-step">继续：<a href="styling.html">设置配色、图例和主题 →</a></div>
