<div class="task-intro"><p>按序列或链调整几何，保留旧版灵活参数的实用性，同时明确它不适用于所有参数。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。

{{prepare:setup,sequencepalette,docstyle}}

## 这个能力仍然保留 {#recommendation}
同一图里常要单独调整某条序列，或分别调整正负链。为这些几何参数保留简洁的
逐序列/逐链写法是合理的；没有必要为了多一条序列而复制整个图层。
推荐先用单值，确实需要差异时用按 `accver` 命名的向量；正负链使用 `+`/`-` 命名向量。

{{code:parameters}}

`radii` 的名称对应实际序列 ID，
`widths` 和位移对应链。重新排序时，命名值仍能找到原序列；无名值只按当前位置对应。
这个例子修改的是布局和几何，未把颜色或 theme 参数当成同一套扩展语法。

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

## 兼容写法与取舍 {#compatibility}
当前还兼容等长无名向量、按序列名的列表、按序号的列表及等长无名列表；gene/注释
的一项还可包含 `+`/`-` 两个值。它们是现有几何接口的兼容形式，不是任意参数的
通用容器。不同参数仍有自己的类型、长度、范围和单位检查。

建议不再扩展低频等价形式：位置列表在筛选、排序后容易错配，嵌套列表阅读成本也高。
本版不删除它们；原定 v0.16.6 的审计会核对是否需要收敛，并在决定破坏性变化前
说明迁移和删除版本。本文的建议不等于已决定弃用。

数值并不共用一个单位：`seq_gap` 是比例，旋转角度以度表示，物理轨道尺寸使用
明确的 `unit()`。`seq_backbone_gap`、`seq_backbone_width` 是单值样式参数，不在
上述列表规则中。共线性自动方向用 `coord_collinear(orientation = "auto")`；
`geom_seq(seq_orientation = ...)` 指定的是 +1/-1 数值方向。

## 完成后检查 {#check-result}
命名向量的键应对应实际 accver 或 +/− 链。重排后再次检查匹配；不要将列表规则推广到 theme、scale 或 aesthetic。

<div class="next-step">继续：<a href="styling.html">设置配色、图例和主题 →</a></div>
