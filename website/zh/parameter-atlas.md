本页按用途汇总当前完整函数参数，不替代学习课程。默认表达式来自实际 formals，说明来自同版本 Rd；NULL 可能在布局阶段解析为默认值。候选向量既可能供 match.arg 选择，也可能表示同时启用多个选项，须阅读说明。上游重导出遵循其自身帮助。

优先从[任务索引](task-index.html)进入教学，重要参数对照见对应课程，完整签名和额外参数转发见函数参考。以下按需展开；不会把上游 ... 的全部可能参数伪装成 ggchord 的独立 API。


## 坐标与全局视图 {#coordinates}

<span id="api-coord_chord"></span>
<details><summary><code>coord_chord()</code></summary>


[完整参考](reference/coord_chord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `rotation` | `0` | 数值/整数；范围见说明 | 全局布局相对 x 轴正方向逆时针旋转的角度，单位为度，默认 0。旋转整个布局，不改变生物学链或序列方向。 |
| `ratio` | `1` | 数值/整数；范围见说明 | 固定 y/x 纵横比，默认 1。 |
| `xlim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的用户坐标范围；显式范围优先于自动拟合范围。 |
| `ylim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的用户坐标范围；显式范围优先于自动拟合范围。 |
| `expand` | `FALSE` | 逻辑值 | 逻辑值，是否扩展坐标范围，默认 FALSE。自动拟合已包含少量安全边距；设为 TRUE 可请求 ggplot2 的额外坐标扩展。 |
| `clip` | `"off"` | 字符/枚举；取值见说明 | 是否将绘图裁剪到面板内，默认 "off"。 |
| `fit` | `c("labels", "geometry", "manual")` | 向量；候选/多选含义见说明 | 拟合策略："labels" 包含测量后的标签框，"geometry" 只拟合几何元素，"manual" 要求显式提供 xlim 和 ylim。 |

</details>

<span id="api-coord_circular"></span>
<details><summary><code>coord_circular()</code></summary>


[完整参考](reference/coord_circular.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `rotation` | `0` | 数值/整数；范围见说明 | 相对基因组原点位于 12 点钟的基准逆时针旋转的角度，单位为度。默认 0 将起点保持在 12 点钟；90 对应九点钟，-90 对应三点钟。 |
| `direction` | `c("clockwise", "counterclockwise")` | 向量；候选/多选含义见说明 | 基因组坐标增加的方向。 |
| `gap` | `0` | 数值/整数；范围见说明 | 开口角度，单位为度；零表示闭合圆。 |
| `ratio` | `1` | 数值/整数；范围见说明 | 固定 y/x 纵横比，默认 1。 |
| `xlim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的用户坐标范围；显式范围优先于自动拟合范围。 |
| `ylim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的用户坐标范围；显式范围优先于自动拟合范围。 |
| `expand` | `FALSE` | 逻辑值 | 逻辑值，是否扩展坐标范围，默认 FALSE。自动拟合已包含少量安全边距；设为 TRUE 可请求 ggplot2 的额外坐标扩展。 |
| `clip` | `"off"` | 字符/枚举；取值见说明 | 是否将绘图裁剪到面板内，默认 "off"。 |
| `fit` | `c("labels", "geometry", "manual")` | 向量；候选/多选含义见说明 | 拟合策略："labels" 包含测量后的标签框，"geometry" 只拟合几何元素，"manual" 要求显式提供 xlim 和 ylim。 |

</details>

<span id="api-coord_collinear"></span>
<details><summary><code>coord_collinear()</code></summary>


[完整参考](reference/coord_collinear.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `ratio` | `1` | 数值/整数；范围见说明 | 固定 y/x 纵横比，默认 1。 |
| `xlim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的用户坐标范围；显式范围优先于自动拟合范围。 |
| `ylim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的用户坐标范围；显式范围优先于自动拟合范围。 |
| `expand` | `FALSE` | 逻辑值 | 逻辑值，是否扩展坐标范围，默认 FALSE。自动拟合已包含少量安全边距；设为 TRUE 可请求 ggplot2 的额外坐标扩展。 |
| `clip` | `"off"` | 字符/枚举；取值见说明 | 是否将绘图裁剪到面板内，默认 "off"。 |
| `fit` | `c("labels", "geometry", "manual")` | 向量；候选/多选含义见说明 | 拟合策略："labels" 包含测量后的标签框，"geometry" 只拟合几何元素，"manual" 要求显式提供 xlim 和 ylim。 |
| `align` | `c("left", "center", "right")` | 向量；候选/多选含义见说明 | 显示对齐方式：默认 "left"，或 "center"、"right"；作用于整条序列或 anchor 指定区间。 |
| `anchor` | `NULL` | 显式窗口/锚点表或已记录的等价输入 | 可选锚点，可为按 accver 命名的数值向量，或含 accver 与 position 或 start/end 的数据框。每条显示序列恰好一个锚点；缺失或重复会报错。多拷贝目标须先明确选择。 |
| `orientation` | `c("manual", "auto")` | 向量；候选/多选含义见说明 | 默认 "manual" 使用 geom_seq() 提供的方向；"auto" 根据区间连接证据选择显示方向，记录支持及未解决情况。 |
| `min_support` | `0.1` | 数值/整数；范围见说明 | 自动改变方向所需的最小不平衡程度，范围 [0, 1]。正反支持的差值绝对值除以总和，默认 0.1；平局始终视为未解决。 |
| `segments` | `NULL` | 显式窗口/锚点表或已记录的等价输入 | 可选数据框，含唯一窗口 accver、row_id 及 orientation（1 或 -1），须声明全部显示 accession。每个窗口有独立局部坐标和可逆源元数据。 |
| `segment_gap` | `100` | 数值/整数；范围见说明 | 提供 segments 时的正显示间隙，单位为 原生序列单位，默认 100，不表示基因组连续。 |
| `seq_axis` | `c("horizontal", "vertical")` | 向量；候选/多选含义见说明 | "horizontal"（默认）或 "vertical"。纵向模式将序列局部几何投影为直立列，并独立保持文字可读方向，不改变源链。 |

</details>

<span id="api-coord_dotplot"></span>
<details><summary><code>coord_dotplot()</code></summary>


[完整参考](reference/coord_dotplot.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `x_seqs` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的唯一accession向量，按轴顺序排列。默认使用比对层中出现的query/subject accession；没有比对时使用全部序列。轴集合之外的原行会省略并报告；自比对和重复HSP原行保留。 |
| `y_seqs` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的唯一accession向量，按轴顺序排列。默认使用比对层中出现的query/subject accession；没有比对时使用全部序列。轴集合之外的原行会省略并报告；自比对和重复HSP原行保留。 |
| `x_orientation` | `1` | 数值/整数；范围见说明 | 显示方向为1或-1，可提供单值，或按全部已选accession命名的向量。输入不变。 |
| `y_orientation` | `1` | 数值/整数；范围见说明 | 显示方向为1或-1，可提供单值，或按全部已选accession命名的向量。输入不变。 |
| `gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 轴区段之间的非负原生序列单位间距。NULL使用最长序列的百分之五。 |
| `ratio` | `1` | 数值/整数；范围见说明 | 物理y单位与物理x单位的比例，默认1。 |
| `xlim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的全局拼接位置范围，仅缩放面板，不修改导出的HSP端点。 |
| `ylim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的全局拼接位置范围，仅缩放面板，不修改导出的HSP端点。 |
| `expand` | `FALSE` | 逻辑值 | 逻辑值，是否扩展面板范围，默认FALSE。 |
| `clip` | `"on"` | 字符/枚举；取值见说明 | 面板裁剪："on"为默认，或"off"。 |

</details>

<span id="api-coord_radial_synteny"></span>
<details><summary><code>coord_radial_synteny()</code></summary>


[完整参考](reference/coord_radial_synteny.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `start_angle` | `90` | 数值/整数；范围见说明 | 首条向外辐条的方向角，单位为度，从正x轴逆时针测量。默认90从12点钟开始。 |
| `clockwise` | `TRUE` | 逻辑值 | 逻辑值，后续辐条是否按顺时针排列，默认TRUE。 |
| `inner_radius` | `0.6` | 数值/整数；范围见说明 | 正的布局内孔半径，默认0.6。按真实比例展示时最长辐条长2。中心孔仅表示显示留白，不移除源位置。 |
| `ratio` | `1` | 数值/整数；范围见说明 | 固定y/x纵横比，默认1。 |
| `xlim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选用户范围，显式范围优先于自动拟合范围。 |
| `ylim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选用户范围，显式范围优先于自动拟合范围。 |
| `expand` | `FALSE` | 逻辑值 | 逻辑值，是否扩展坐标范围，默认FALSE。自动拟合已有小幅安全边距；设为TRUE可再请求ggplot2坐标扩展。 |
| `clip` | `"off"` | 字符/枚举；取值见说明 | 是否将绘制裁剪到面板，默认"off"。 |
| `fit` | `c("labels", "geometry", "manual")` | 向量；候选/多选含义见说明 | 拟合策略："labels"包含实测标签框，"geometry"只拟合几何元素，"manual"要求显式xlim及ylim。 |

</details>

<span id="api-coord_sankey"></span>
<details><summary><code>coord_sankey()</code></summary>


[完整参考](reference/coord_sankey.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `column_gap` | `3` | 数值/整数；范围见说明 | 两列骨架的正间距，使用布局单位。 |
| `self_links` | `FALSE` | 逻辑值 | 是否显示 query 和 subject 实体相同的关系。 |
| `ratio` | `1` | 数值/整数；范围见说明 | 固定 y/x 纵横比，默认 1。 |
| `xlim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选用户范围，显式范围优先于自动拟合范围。 |
| `ylim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选用户范围，显式范围优先于自动拟合范围。 |
| `expand` | `FALSE` | 逻辑值 | 是否扩展坐标范围，默认 FALSE。自动拟合已含小幅安全边距；设为 TRUE 可再请求 ggplot2 坐标扩展。 |
| `clip` | `"off"` | 字符/枚举；取值见说明 | 是否将绘图裁剪到 panel，默认 "off"。 |
| `fit` | `c("labels", "geometry", "manual")` | 向量；候选/多选含义见说明 | 拟合策略："labels" 包含实测标签框，"geometry" 只拟合几何，"manual" 要求显式 xlim 和 ylim。 |

</details>


## 几何与标签 {#layers}

<span id="api-geom_alignment_segment"></span>
<details><summary><code>geom_alignment_segment()</code></summary>


[完整参考](reference/geom_alignment_segment.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准ggplot2映射与图层数据，默认使用构造时的ribbon_data。几何角色为qaccver、saccver、qstart、qend、sstart和send。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准ggplot2映射与图层数据，默认使用构造时的ribbon_data。几何角色为qaccver、saccver、qstart、qend、sstart和send。 |
| `colour` | `"#59636D"` | 字符/枚举；取值见说明 | 固定样式。映射样式使用link_colour、link_alpha、link_linewidth和link_linetype，配合对应的scale_link函数。显式固定值优先。 |
| `alpha` | `0.6` | 数值/整数；范围见说明 | 固定样式。映射样式使用link_colour、link_alpha、link_linewidth和link_linetype，配合对应的scale_link函数。显式固定值优先。 |
| `linewidth` | `0.4` | 数值/整数；范围见说明 | 固定样式。映射样式使用link_colour、link_alpha、link_linewidth和link_linetype，配合对应的scale_link函数。显式固定值优先。 |
| `linetype` | `1` | 数值/整数；范围见说明 | 固定样式。映射样式使用link_colour、link_alpha、link_linewidth和link_linetype，配合对应的scale_link函数。显式固定值优先。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 仅支持"identity"；数据位置调整会改变比对端点。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准ggplot2图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准ggplot2图层参数。 |
| `allow_repeated` | `FALSE` | 逻辑值 | 是否允许显式标识的source_hsp_id记录重复展示，所有参与层均须选择允许。无声明ID的行不按值去重。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传递给ggplot2线段几何的其他参数。 |

</details>

<span id="api-geom_feature"></span>
<details><summary><code>geom_feature()</code></summary>


[完整参考](reference/geom_feature.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 可选映射。accver、start、end、strand 等角色 aesthetic 可重命名输入列；普通视觉映射在生成几何后计算。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 含 accver、start、end、strand 的数据框，可选 type、category、label。因子列 start/end 使用数值标签而非内部级别编码；坐标必须为有限数值。可提供 display_priority 或逻辑列 prioritized_display。也可用 directionality 或 direction 代替 strand，取值为 "forward"、"reverse"、"bidirectional"、"none"。 省略时使用构造函数的完整 feature_data（或兼容的 gene_data），其他图层不会从默认输入删除记录。 |
| `shape` | `"arrow"` | 字符/枚举；取值见说明 | 固定注释形状。显式 shape 优先于 aes(feature_shape = ...)：支持 "arrow"、"compact_arrow"、"promoter_arrow"、"primer_arrow"、"marker"、"block"、"capped_line"、"primer_arc"、"chevron"、"lollipop"，默认 "arrow"。使用 aes(feature_shape = type) 与 scale_feature_shape_manual() 将类别映射为形状。 |
| `width` | `NULL` | 几何数值/向量/已记录格式 | 可选数值、向量或列表，控制共享注释几何引擎的宽度。圆图中与标签组合时，基础宽度会在限定范围内响应最终字号。 |
| `arrow_head_length` | `0.04` | 数值/整数；范围见说明 | 序列局部坐标中的箭头尺寸，与 geom_gene() 共用。arrow_head_width 大于 1 时形成凸出的块状箭头。 |
| `arrow_head_width` | `1` | 数值/整数；范围见说明 | 序列局部坐标中的箭头尺寸，与 geom_gene() 共用。arrow_head_width 大于 1 时形成凸出的块状箭头。 |
| `arrow_head_style` | `c("shouldered", "flush", "triangle")` | 向量；候选/多选含义见说明 | 箭头构造方式："shouldered" 保留凸出箭头，"flush" 使箭头与主体齐平，"triangle" 绘制完整三角楔形。 |
| `short_feature` | `c("auto", "wedge", "block")` | 向量；候选/多选含义见说明 | 注释太短、容不下指定箭头时的处理：自动压缩为方向符号，或改用楔形、矩形。 |
| `segment_boundaries` | `TRUE` | 逻辑值 | 绘制分段注释携带的内部标记。连续源分段的连接处使用点线；分段级非实线 line_style 可覆盖该样式，cleavage_arrows 坐标可添加或合并生物学切割标记。这些标记不会产生额外箭头。 |
| `boundary_colour` | `"#666666"` | 字符/枚举；取值见说明 | 内部分段连接处的外观。 |
| `boundary_linewidth` | `0.25` | 数值/整数；范围见说明 | 内部分段连接处的外观。 |
| `boundary_linetype` | `"22"` | 字符/枚举；取值见说明 | 内部分段连接处的外观。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 注释位置，使用 "identity"、"strand"、"plasmid"、position_strand()、position_plasmid() 或 position_feature_stack()。 |
| `show.legend` | `TRUE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他固定多边形 aesthetic。 |

</details>

<span id="api-geom_feature_label"></span>
<details><summary><code>geom_feature_label()</code></summary>


[完整参考](reference/geom_feature_label.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入；通过 label 映射注释文字。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入；通过 label 映射注释文字。 |
| `label_orientation` | `c("feature", "horizontal", "radial", "tangent")` | 向量；候选/多选含义见说明 | 支持 "feature"、"horizontal"、"radial"、"tangent"。"feature" 使内部及相邻标签均沿注释切线；较紧凑的文字移到注释旁。 |
| `label_side` | `c("outside", "inside", "auto")` | 向量；候选/多选含义见说明 | 支持 "inside"、"outside" 或 "auto"。 |
| `label_overlap` | `c("hide", "nudge", "allow")` | 向量；候选/多选含义见说明 | 固定标签的碰撞策略。 |
| `label_wrap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的字符换行宽度。 |
| `label_rotation` | `NULL` | 几何数值/向量/已记录格式 | 可选数值、向量或列表，旋转角度，默认 0。 |
| `label_radial_offset` | `0.04` | 数值/整数；范围见说明 | 可选数值、向量或列表，径向偏移，默认 0.04。 |
| `label_circum_offset` | `NULL` | 几何数值/向量/已记录格式 | 可选数值、向量或列表，周向偏移，默认 0。 |
| `label_circum_limit` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选逻辑值、向量或列表，是否限制周向偏移，默认 TRUE。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标签锚点使用的注释 Position。在 coord_collinear() 中，省略时共享匹配注释层的轨道；显式 Position 只改变标签锚点，不移动多边形。若两组输入均有唯一 feature_id，标签可取子集或重排，但序列标识、区间和链须一致。没有 ID 时，继承要求区间行及顺序完全相同。多个匹配几何层需要显式标签 Position 或无歧义注释 ID。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他文字参数。 |

</details>

<span id="api-geom_feature_label_repel"></span>
<details><summary><code>geom_feature_label_repel()</code></summary>


[完整参考](reference/geom_feature_label_repel.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入；通过 label 映射注释文字。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入；通过 label 映射注释文字。 |
| `label_layout` | `c("feature", "radial", "auto", "callout")` | 向量；候选/多选含义见说明 | 布局支持 "feature"、"radial"、"auto"、"callout"。默认 "feature" 遵循质粒图习惯：先尝试内部，再放到固定半径的内部轨道。所有非外部文字均按选定轨道半径逐字排布，而不是整体旋转成直线；注释与标签之间隔着一个或多个轨道时绘制浅色引导线。需要不受限制的外部引出标签时选择其他模式。 |
| `label_wrap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 文字适配控制。显式 label_fit = "shrink" 优先尝试内部位置，不截断原文；字号下限由 label_min_size 决定。仍无法解决的标签保持请求字号，放在相邻/外侧或隐藏；其他模式维持已有行为。 |
| `label_fit` | `c("wrap", "none", "ellipsis", "auto", "shrink")` | 向量；候选/多选含义见说明 | 文字适配控制。显式 label_fit = "shrink" 优先尝试内部位置，不截断原文；字号下限由 label_min_size 决定。仍无法解决的标签保持请求字号，放在相邻/外侧或隐藏；其他模式维持已有行为。 |
| `label_max_lines` | `2L` | 数值/整数；范围见说明 | 文字适配控制。显式 label_fit = "shrink" 优先尝试内部位置，不截断原文；字号下限由 label_min_size 决定。仍无法解决的标签保持请求字号，放在相邻/外侧或隐藏；其他模式维持已有行为。 |
| `label_side` | `c("outside", "inside", "auto")` | 向量；候选/多选含义见说明 | 支持 "inside"、"outside" 或 "auto"。 |
| `max_overlaps` | `Inf` | R 对象；约束见说明 | 允许尚未解决的最大重叠数。 |
| `external` | `TRUE` | 逻辑值 | 内部和相邻位置都无法容纳后，是否允许放到圆形骨架外。FALSE 时隐藏无法解决的标签，不无限向内推。外部标签使用注释最终填色的浅色圆角框；与圆图酶切标签共存时，各自保留布局语义，但在坐标系共同管理的外部区域登记最终框，防止占用同一空间。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标签锚点使用的注释 Position。在 coord_collinear() 中，省略时共享匹配注释层的轨道；显式 Position 只改变标签锚点，不移动多边形。若两组输入均有唯一 feature_id，标签可取子集或重排，但序列标识、区间和链须一致。没有 ID 时，继承要求区间行及顺序完全相同。多个匹配几何层需要显式标签 Position 或无歧义注释 ID。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `label_min_size` | `1.5` | 数值/整数；范围见说明 | 显式缩小模式的最小字号，单位 mm，默认 1.5，不得超过请求字号。导出标签保留原文、请求/最终字号与适配原因。 |
| `label_segment_overlap` | `"fade"` | 字符/枚举；取值见说明 | 引导线被其他标签覆盖部分的绘制方式，默认 "fade"。"fade" 降低透明度，"clip" 删除覆盖部分，"show" 保留完整线。不会影响换行、省略或隐藏。 |
| `label_segment_overlap_alpha` | `0.18` | 数值/整数；范围见说明 | [0, 1] 内的数值，默认 0.18；仅在 label_segment_overlap = "fade" 时控制被覆盖引导线部分的相对透明度。 |
| `label_segment_linetype` | `"auto"` | 字符/枚举；取值见说明 | 字符或数值，默认 "auto"。"auto" 通常用实线，移到弧线另一侧的标签用虚线；其他合法 ggplot2 线型（如 "solid"、"dashed"、"dotted" 或数字虚线模式）应用于全部引导线。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他文字参数。 |

</details>

<span id="api-geom_feature_plasmid"></span>
<details><summary><code>geom_feature_plasmid()</code></summary>


[完整参考](reference/geom_feature_plasmid.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 可选映射。accver、start、end、strand 等角色 aesthetic 可重命名输入列；普通视觉映射在生成几何后计算。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 含 accver、start、end、strand 的数据框，可选 type、category、label。因子列 start/end 使用数值标签而非内部级别编码；坐标必须为有限数值。可提供 display_priority 或逻辑列 prioritized_display。也可用 directionality 或 direction 代替 strand，取值为 "forward"、"reverse"、"bidirectional"、"none"。 省略时使用构造函数的完整 feature_data（或兼容的 gene_data），其他图层不会从默认输入删除记录。 |
| `shape` | `"arrow"` | 字符/枚举；取值见说明 | 固定注释形状。显式 shape 优先于 aes(feature_shape = ...)：支持 "arrow"、"compact_arrow"、"promoter_arrow"、"primer_arrow"、"marker"、"block"、"capped_line"、"primer_arc"、"chevron"、"lollipop"，默认 "arrow"。使用 aes(feature_shape = type) 与 scale_feature_shape_manual() 将类别映射为形状。 |
| `width` | `0.058` | 数值/整数；范围见说明 | 轨道宽度，默认对应桌面质粒编辑器常用的紧凑注释带。 |
| `arrow_head_length` | `0.068` | 数值/整数；范围见说明 | 序列局部坐标中的箭头尺寸，与 geom_gene() 共用。arrow_head_width 大于 1 时形成凸出的块状箭头。 |
| `arrow_head_width` | `1.42` | 数值/整数；范围见说明 | 序列局部坐标中的箭头尺寸，与 geom_gene() 共用。arrow_head_width 大于 1 时形成凸出的块状箭头。 |
| `arrow_head_style` | `"shouldered"` | 字符/枚举；取值见说明 | 箭头构造方式："shouldered" 保留凸出箭头，"flush" 使箭头与主体齐平，"triangle" 绘制完整三角楔形。 |
| `short_feature` | `"auto"` | 字符/枚举；取值见说明 | 注释太短、容不下指定箭头时的处理：自动压缩为方向符号，或改用楔形、矩形。 |
| `segment_boundaries` | `TRUE` | 逻辑值 | 绘制分段注释携带的内部标记。连续源分段的连接处使用点线；分段级非实线 line_style 可覆盖该样式，cleavage_arrows 坐标可添加或合并生物学切割标记。这些标记不会产生额外箭头。 |
| `boundary_colour` | `"#666666"` | 字符/枚举；取值见说明 | 内部分段连接处的外观。 |
| `boundary_linewidth` | `0.25` | 数值/整数；范围见说明 | 内部分段连接处的外观。 |
| `boundary_linetype` | `"22"` | 字符/枚举；取值见说明 | 内部分段连接处的外观。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | 注释位置。默认 NULL 在质粒骨架内侧将重叠注释堆叠到紧凑轨道。独立标签层应复用同一个显式 position_feature_stack() 对象，使多边形与标签同轨。可选数值列 display_priority（较大优先）或逻辑列 prioritized_display 可覆盖重叠注释默认的最长优先顺序。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他固定多边形 aesthetic。 |

</details>

<span id="api-geom_gene"></span>
<details><summary><code>geom_gene()</code></summary>


[完整参考](reference/geom_gene.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 默认 NULL，使用预计算数据。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选图层表或数据函数，优先于构造数据。使用统一 feature_data 时，默认只选择明确分类的 gene/CDS 记录，见 ggchord()。旧 gene_data 保留全部行。 |
| `width` | `NULL` | 几何数值/向量/已记录格式 | 可选数值、向量或列表，基因箭头宽度，默认 0.05。 |
| `arrow_head_length` | `0.04` | 数值/整数；范围见说明 | 序列局部坐标中的共享箭头尺寸。 |
| `arrow_head_width` | `1` | 数值/整数；范围见说明 | 序列局部坐标中的共享箭头尺寸。 |
| `arrow_head_style` | `c("shouldered", "flush", "triangle")` | 向量；候选/多选含义见说明 | 箭头形状：凸肩、齐平或完整三角楔形。 |
| `short_feature` | `c("auto", "wedge", "block")` | 向量；候选/多选含义见说明 | 基因太短、容不下指定箭头时，自动选择楔形/矩形，或明确使用其中一种。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 注释位置，使用 "identity"、"strand"、"plasmid"、position_strand()、position_plasmid() 或 position_feature_stack()。 |
| `show.legend` | `TRUE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_polygon() 的其他参数。 |

</details>

<span id="api-geom_gene_label"></span>
<details><summary><code>geom_gene_label()</code></summary>


[完整参考](reference/geom_gene_label.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 默认 NULL，使用预计算数据。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选图层表或数据函数，默认数据按 ggchord() 所述选择 gene/CDS。 |
| `label_orientation` | `"horizontal"` | 字符/枚举；取值见说明 | 文字相对序列路径的方向，默认 "horizontal"，支持 "radial"、"tangent" 或 "horizontal"。 |
| `label_side` | `"outside"` | 字符/枚举；取值见说明 | 标签侧向，默认 "outside"。可放在弦图外侧、按链选择 "auto" 侧，或强制 "inside"。 |
| `label_overlap` | `"hide"` | 字符/枚举；取值见说明 | 固定标签碰撞策略，默认 "hide"。"hide" 省略后出现的冲突标签，"nudge" 小幅移动位置，"allow" 保留全部请求位置。因此 "hide" 模式下可通过输入行顺序设置优先级。 |
| `label_rotation` | `NULL` | 几何数值/向量/已记录格式 | 可选数值、向量或列表，旋转角度，默认 0。 |
| `label_radial_offset` | `0.04` | 数值/整数；范围见说明 | 可选数值、向量或列表，径向偏移，默认 0.04。 |
| `label_circum_offset` | `NULL` | 几何数值/向量/已记录格式 | 可选数值、向量或列表，周向偏移，默认 0。 |
| `label_circum_limit` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选逻辑值、向量或列表，是否限制周向偏移，默认 TRUE。 |
| `label_wrap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 数值或 NULL，默认 NULL。按指定字符数（如 15）换行，使长注释更窄、减少重叠。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 位置调整。标签使用独立数据或没有几何层时，应提供与基因/注释层相同的 position_feature_stack()。在 coord_collinear() 中，未显式提供 Position 的标签共享匹配几何层的轨道；显式标签 Position 不移动几何。唯一 feature_id 列允许取子集或重排，但标识、区间和链必须一致。多个匹配几何层需要显式标签 Position 或无歧义注释 ID。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_text() 的其他参数。 |

</details>

<span id="api-geom_gene_label_repel"></span>
<details><summary><code>geom_gene_label_repel()</code></summary>


[完整参考](reference/geom_gene_label_repel.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 默认 NULL，使用预计算数据。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选图层表或数据函数，默认数据按 ggchord() 所述选择 gene/CDS。 |
| `label_layout` | `"radial"` | 字符/枚举；取值见说明 | 标签布局，默认 "radial"。"radial" 沿序列相对轮廓排列水平外部文字，引导线沿法线直出；移位时先沿法线短出，再用长线连接。"auto" 在左右共享竖列，其余位置使用径向布局。"arc" 沿序列切线旋转文字。旧 "aligned" 模式已删除，侧列请使用 "auto"。 |
| `label_wrap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 数值或 NULL，默认 NULL。按指定字符数（如 15）对长基因注释换行，覆盖 label_fit。 |
| `label_fit` | `"wrap"` | 字符/枚举；取值见说明 | 文字适配方式，默认 "wrap"，在确定性布局前处理拥挤或过宽的初始框。"wrap" 换行保留完整注释，"ellipsis" 缩成一行，"auto" 换行且仅超过 label_max_lines 时加省略号，"none" 关闭自适应。显式 "shrink" 优先内部放置，完整文字最多缩小到 label_min_size。仍无法容纳时，按现有策略以请求字号移到旁边/外侧或隐藏。该选项使用注释布局，而不是 label_layout 选择的外部布局；其他适配方式不变。 |
| `label_max_lines` | `2L` | 数值/整数；范围见说明 | 正整数，默认 2；自适应 "wrap" 或 "auto" 允许的最大行数。 |
| `label_side` | `"outside"` | 字符/枚举；取值见说明 | 标签位于弧线哪一侧，默认 "outside"。"auto" 保持按链放置；"outside" 将本会位于弦图内部、可能覆盖连接的标签移到各弧外侧；"inside" 相反。换侧标签使用虚线引导，见 label_segment_linetype。 |
| `max_overlaps` | `Inf` | R 对象；约束见说明 | 数值，默认 Inf。避让后仍与超过此数量的其他标签重叠时隐藏，类似 ggrepel 的去拥挤策略；有限值可简化拥挤图。 |
| `label_segment_overlap` | `"fade"` | 字符/枚举；取值见说明 | 引导线被其他标签覆盖部分的绘制方式，默认 "fade"。"fade" 降低透明度，"clip" 删除覆盖部分，"show" 保留完整线。不会影响换行、省略或隐藏。 |
| `label_segment_overlap_alpha` | `0.18` | 数值/整数；范围见说明 | [0, 1] 内的数值，默认 0.18；仅在 label_segment_overlap = "fade" 时控制被覆盖引导线部分的相对透明度。 |
| `label_segment_linetype` | `"auto"` | 字符/枚举；取值见说明 | 字符或数值，默认 "auto"。"auto" 通常用实线，移到弧线另一侧的标签用虚线；其他合法 ggplot2 线型（如 "solid"、"dashed"、"dotted" 或数字虚线模式）应用于全部引导线。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 位置调整。标签使用独立数据或没有几何层时，提供与基因/注释层相同的 position_feature_stack()。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `label_min_size` | `1.5` | 数值/整数；范围见说明 | 显式 label_fit = "shrink" 的最小字号，单位 mm，默认 1.5，不得超过请求字号。导出标签几何保留原文、请求/最终字号及适配原因。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他文字样式及绘制选项。自动布局控制旋转、偏移、方向和碰撞策略；手动 label_rotation、label_radial_offset、label_circum_offset、label_circum_limit、label_orientation、label_overlap 应使用 geom_gene_label()，本层会拒绝这些参数。 |

</details>

<span id="api-geom_link_line"></span>
<details><summary><code>geom_link_line()</code></summary>


[完整参考](reference/geom_link_line.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准 ggplot2 映射及图层数据，支持数据函数。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准 ggplot2 映射及图层数据，支持数据函数。 |
| `link_type` | `c("curve", "straight")` | 向量；候选/多选含义见说明 | 连接形状：默认 "curve" 曲线，或 "straight" 直线。 |
| `link_gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 显式端点间距。默认 link_avoid = "none" 时，NULL 使用固定间距 0.035；将 link_avoid 设为 "smooth" 或 "uniform" 可启用考虑图形实体的间距。 |
| `link_ctrl_point` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 贝塞尔控制点：两个数、四个数或逐输入行列表。 |
| `link_avoid` | `c("none", "smooth", "uniform")` | 向量；候选/多选含义见说明 | 避障方式："none"（默认）、"smooth" 或 "uniform"；显式间距优先。 |
| `link_branch` | `c("none", "query", "subject")` | 向量；候选/多选含义见说明 | 共享相同 query 或 subject 端点："none"、"query"、"subject"。 |
| `link_branch_fraction` | `0.15` | 数值/整数；范围见说明 | 共享段相对最短成员的长度比例，范围 (0, 0.5]。 |
| `arrow` | `NULL` | 可空 R 对象；具体类型/继承见说明 | grid 箭头定义或 NULL。 |
| `colour` | `"#59636D"` | 字符/枚举；取值见说明 | 固定样式；映射使用 link_* aesthetic。显式固定值优先于对应映射。要保留映射请省略固定参数；inherit.aes = TRUE 时也包括继承的映射。 |
| `alpha` | `0.6` | 数值/整数；范围见说明 | 固定样式；映射使用 link_* aesthetic。显式固定值优先于对应映射。要保留映射请省略固定参数；inherit.aes = TRUE 时也包括继承的映射。 |
| `linewidth` | `0.4` | 数值/整数；范围见说明 | 固定样式；映射使用 link_* aesthetic。显式固定值优先于对应映射。要保留映射请省略固定参数；inherit.aes = TRUE 时也包括继承的映射。 |
| `linetype` | `1` | 数值/整数；范围见说明 | 固定样式；映射使用 link_* aesthetic。显式固定值优先于对应映射。要保留映射请省略固定参数；inherit.aes = TRUE 时也包括继承的映射。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `between` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 两个显示实例名称或声明顺序序号。 |
| `allow_repeated` | `FALSE` | 逻辑值 | 允许显式声明的 source_hsp_id 被重复显示。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他路径参数。 |

</details>

<span id="api-geom_link_ribbon"></span>
<details><summary><code>geom_link_ribbon()</code></summary>


[完整参考](reference/geom_link_ribbon.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 默认 NULL，使用预计算数据。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 默认 NULL，从布局自动获取。 |
| `link_ctrl_point` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的贝塞尔边界控制向量/列表。两个数指定共享控制点，四个数分别指定两条边界的控制点；逐行列表可含成对坐标向量。NULL 在弦图中使用原点 c(0, 0)，在线性图中由序列端点推导控制点。共享分支时，控制点影响远端形状，汇合处跟随共享主干。 |
| `link_gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值/向量，控制序列与连接端点的间距。默认 link_avoid = "uniform" 时，NULL 为每个端点统一选择考虑障碍的间距；"smooth" 仅在 geom_gene()、geom_feature() 多边形附近局部调整，"none" 则使用固定 0.035。忽略文字与引导线。提供数值会关闭自动规则并精确使用该间距。线性布局会缩减自动基础间距以适应窄通道；无法容纳的显式间距或注释避让距离会报错。 |
| `link_avoid` | `c("uniform", "none", "smooth")` | 向量；候选/多选含义见说明 | 避障方式："uniform"（默认）、"none" 或 "smooth"。平滑方式局部调整，统一方式在整个端点采用同一间距。 |
| `link_branch` | `c("none", "query", "subject")` | 向量；候选/多选含义见说明 | 共享相同 query 或 subject 端点："none"、"query"、"subject"。 |
| `link_branch_fraction` | `0.15` | 数值/整数；范围见说明 | 共享段相对最短成员的长度比例，范围 (0, 0.5]。 |
| `fill` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的固定 Ribbon 填充色。NULL 使用默认 Identity 填充。 |
| `alpha` | `0.42` | 数值/整数；范围见说明 | 标准固定样式。显式值优先于对应映射；省略的默认值不会阻止用户映射。例如映射 ribbon_alpha 时应省略 alpha。 |
| `colour` | `"#59636D"` | 字符/枚举；取值见说明 | 标准固定样式。显式值优先于对应映射；省略的默认值不会阻止用户映射。例如映射 ribbon_alpha 时应省略 alpha。 |
| `linewidth` | `0.08` | 数值/整数；范围见说明 | 标准固定样式。显式值优先于对应映射；省略的默认值不会阻止用户映射。例如映射 ribbon_alpha 时应省略 alpha。 |
| `linetype` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 固定边框线型。默认 NULL 对相邻显示序列使用实线，对非相邻序列使用虚线。共线性边框宽度由 nonadjacent_linewidth 控制。弦图首尾视为相邻。相邻关系按最终显示顺序确定，与比对方向无关。显式固定线型或 aes(ribbon_linetype = ...)（也支持 aes(linetype = ...)）覆盖自动线型。数据映射应搭配对应标度，自动线型不会训练标度；inherit.aes = TRUE 时使用全局线型映射。 |
| `adjacent_only` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 是否只连接相邻序列。NULL 在 coord_collinear() 中相当于 TRUE，在 coord_chord() 中相当于 FALSE。设为 FALSE 可保留跨轨道连接；弦图首尾也是相邻。 |
| `pairs` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数据框，以 qaccver 和 saccver 选择无序比较对。adjacent_only = NULL 时，显式比较对优先于坐标系默认值。export_ggchord_layout() 报告省略行。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `TRUE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `nonadjacent_linewidth` | `0.4` | 数值/整数；范围见说明 | 共线性布局 coord_collinear() 中非相邻显示行的固定边框线宽，默认 0.4 mm。这些连接带默认使用虚线。设为 0 时，所有连接带沿用普通 linewidth。此参数只控制线宽，不改变线型；显式 aes(linewidth = ...) 映射优先。 |
| `between` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 两个显示实例名称或声明序号，用于确定唯一路由。 |
| `allow_repeated` | `FALSE` | 逻辑值 | 显式允许声明的 source_hsp_id 出现在多个显示关系中。 |
| `tip` | `c("none", "subject")` | 向量；候选/多选含义见说明 | 默认 none 保留完整区间端面；也可选择 subject 端收尖。肩部在通道 88% 处保留平移后的目标端点跨度。尖点是端面中点的装饰标记，不代表单个匹配碱基或生物转移方向。退化端面会警告并恢复完整端面。源区间保持不变。 |
| `gradient` | `FALSE` | 逻辑值 | 逻辑值。显式启用端点序列色插值，默认 FALSE。显式 fill、映射和标度仍然优先。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_polygon() 的其他参数。 |

</details>

<span id="api-geom_primer"></span>
<details><summary><code>geom_primer()</code></summary>


[完整参考](reference/geom_primer.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入。 |
| `side` | `c("between", "outside", "inside")` | 向量；候选/多选含义见说明 | 引物弧所在骨架轨道。"between" 使用中心线；"outside"、"inside" 使用 arc_offset。 |
| `arc_offset` | `0.045` | 数值/整数；范围见说明 | side 为 "outside" 或 "inside" 时相对骨架中心的距离。 |
| `arc_linewidth` | `0.026` | 数值/整数；范围见说明 | 引物弧的径向宽度。 |
| `colour` | `"#A020F0"` | 字符/枚举；取值见说明 | 引物弧颜色。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | 标准图层参数。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_feature() 的其他参数。 |

</details>

<span id="api-geom_primer_label_repel"></span>
<details><summary><code>geom_primer_label_repel()</code></summary>


[完整参考](reference/geom_primer_label_repel.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入。 |
| `side` | `c("outside", "inside")` | 向量；候选/多选含义见说明 | 引物引出标签相对于骨架所在侧。"inside" 保留内侧标签布局，包括图中共同存在外侧酶切或注释标签时。默认为 "outside"。 |
| `show_location` | `TRUE` | 逻辑值 | 是否在标签中包含结合范围。 |
| `colour` | `"#A020F0"` | 字符/枚举；取值见说明 | 引物弧颜色。 |
| `max_overlaps` | `Inf` | R 对象；约束见说明 | 允许尚未解决的最大重叠数。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | 标准图层参数。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_feature() 的其他参数。 |

</details>

<span id="api-geom_residue"></span>
<details><summary><code>geom_residue()</code></summary>


[完整参考](reference/geom_residue.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准 ggplot2 图层映射与数据。必需美学为 x、y 和 label；支持普通 colour、fill、alpha、size、family 和 fontface 映射与标度。size 的单位为 mm。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准 ggplot2 图层映射与数据。必需美学为 x、y 和 label；支持普通 colour、fill、alpha、size、family 和 fontface 映射与标度。size 的单位为 mm。 |
| `overflow` | `c("tile", "error")` | 向量；候选/多选含义见说明 | 字母无法容纳时使用 "tile"（默认）或 "error"。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `TRUE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给图层的固定文字/色格美学及 na.rm。 |

</details>

<span id="api-geom_restriction_site"></span>
<details><summary><code>geom_restriction_site()</code></summary>


[完整参考](reference/geom_restriction_site.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入。 |
| `label` | `TRUE` | 逻辑值 | 是否绘制标签。 |
| `label_style` | `c("enzyme_position", "enzyme", "position")` | 向量；候选/多选含义见说明 | 标签显示酶及位置、仅酶名或仅位置。 |
| `label_order` | `c("auto", "enzyme_position", "position_enzyme")` | 向量；候选/多选含义见说明 | 自动镜像顺序、酶名优先或位置优先。 |
| `label_side` | `c("outside", "inside", "auto")` | 向量；候选/多选含义见说明 | 标签放在外侧、内侧或自动选择。 |
| `leader` | `c("radial", "elbow", "straight")` | 向量；候选/多选含义见说明 | 引导线路径。默认 "radial" 测量实际文字框，分开位点锚点与最终径向位置，使用直连，或独立径向短线加扇形连接。所有面向酶位点的文字边缘按基因组顺序排列在同一紧凑圆形轮廓。密集侧向簇使用等距轮廓扇形，稀疏簇保留自然角度。旧 "trunk" 别名在 v0.16.1 删除，改用 "radial"。"straight" 从刻度尖端直连最终文字边缘，不绕行序列或其他线；"elbow" 使用一个折点，尽可能先沿局部法线出发。所有路径共享标签位置及分组。 |
| `min_label_gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 标签槽之间的最小基因组比例，NULL 根据实际标签和设备推导紧凑默认值。 |
| `tick_length` | `0.04` | 数值/整数；范围见说明 | 局部法线距离。 |
| `label_offset` | `0.18` | 数值/整数；范围见说明 | 局部法线距离。 |
| `colour` | `"#202020"` | 字符/枚举；取值见说明 | 固定外观。默认 fontface = NULL 时，仅在该酶对完整输入切割一次时加粗酶名，坐标仍为常规字重；显式 fontface 关闭此自动区分。 |
| `linewidth` | `0.28` | 数值/整数；范围见说明 | 固定外观。默认 fontface = NULL 时，仅在该酶对完整输入切割一次时加粗酶名，坐标仍为常规字重；显式 fontface 关闭此自动区分。 |
| `label_size` | `3.2` | 数值/整数；范围见说明 | 固定外观。默认 fontface = NULL 时，仅在该酶对完整输入切割一次时加粗酶名，坐标仍为常规字重；显式 fontface 关闭此自动区分。 |
| `fontface` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 固定外观。默认 fontface = NULL 时，仅在该酶对完整输入切割一次时加粗酶名，坐标仍为常规字重；显式 fontface 关闭此自动区分。 |
| `family` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 固定外观。默认 fontface = NULL 时，仅在该酶对完整输入切割一次时加粗酶名，坐标仍为常规字重；显式 fontface 关闭此自动区分。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准图层参数。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他几何参数。 |

</details>

<span id="api-geom_seq"></span>
<details><summary><code>geom_seq()</code></summary>


[完整参考](reference/geom_seq.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 默认 NULL，使用预计算数据。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 默认 NULL，从布局自动获取。 |
| `seq_order` | `NULL` | 几何数值/向量/已记录格式 | 显示顺序中的来源 ID，重复值创建显示实例。可使用 c(A_high='A', B='B', A_low='A') 这样的名称稳定标识实例。重复实体的 HSP 要求每个连接层提供唯一 between 路由。 |
| `seq_orientation` | `NULL` | 几何数值/向量/已记录格式 | 可选数值 1 或 -1，序列方向，默认 1。 |
| `seq_gap` | `NULL` | 几何数值/向量/已记录格式 | 可选数值，序列间隙比例；弦图默认 0.03。在 coord_collinear() 中表示相邻骨架间距相对最长骨架宽度的正比例，默认 0.3；线性布局不使用最后一条序列的 gap。分段线性图以最长显示行衡量，同一行窗口须使用相同 gap，只有行间间隙影响纵向位置。可按序列 ID 命名以指定其后的间隙，未列出的序列使用相应坐标默认值。 |
| `seq_radius` | `NULL` | 几何数值/向量/已记录格式 | 可选正数，序列弧半径，默认 1.0。 |
| `seq_curvature` | `NULL` | 几何数值/向量/已记录格式 | 可选有符号弯曲程度：0 为直线，1 为标准弧，负值向相反方向弯曲，默认 1.0。接受有限正负值且不截断，绝对值大于 1 会增强弯曲。 |
| `seq_style` | `c("auto", "single", "double", "band")` | 向量；候选/多选含义见说明 | 骨架形状。"auto" 在 coord_circular() 中使用双线，在 coord_chord() 中使用单线骨架。 |
| `seq_backbone_gap` | `0.025` | 数值/整数；范围见说明 | 双线骨架之间的间距。 |
| `seq_backbone_width` | `0.035` | 数值/整数；范围见说明 | 带状骨架的宽度。 |
| `linewidth` | `0.9` | 数值/整数；范围见说明 | 骨架线宽，默认 0.9。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `TRUE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `length_mode` | `c("proportional", "equal")` | 向量；候选/多选含义见说明 | "proportional"（默认）保留源长度比例；"equal" 为各显示序列分配相同角度或直线跨度，不修改源长度。多序列 equal 模式不能使用共享绝对比例尺，局部轴保留源单位。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他固定路径 aesthetic。colour（或 color）设置 seq_colour。显式固定值优先于对应映射；要保留映射请省略固定参数。 |

</details>

<span id="api-geom_seq_break"></span>
<details><summary><code>geom_seq_break()</code></summary>


[完整参考](reference/geom_seq_break.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准映射及序列数据，默认使用图的序列表。截断来自 .truncated_left/.truncated_right，可用时以源区间和长度元数据补充。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准映射及序列数据，默认使用图的序列表。截断来自 .truncated_left/.truncated_right，可用时以源区间和长度元数据补充。 |
| `break_width` | `0.02` | 数值/整数；范围见说明 | 整个断裂标记的宽高，单位为线性布局距离，不是 mm。 |
| `break_height` | `0.04` | 数值/整数；范围见说明 | 整个断裂标记的宽高，单位为线性布局距离，不是 mm。 |
| `linewidth` | `0.4` | 数值/整数；范围见说明 | 固定线外观；也可映射普通视觉 aesthetic 并使用 ggplot2 标度。 |
| `colour` | `"#3D4348"` | 字符/枚举；取值见说明 | 固定线外观；也可映射普通视觉 aesthetic 并使用 ggplot2 标度。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::geom_path() 的其他参数。 |

</details>

<span id="api-geom_seq_center_label"></span>
<details><summary><code>geom_seq_center_label()</code></summary>


[完整参考](reference/geom_seq_center_label.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入。 |
| `name` | `TRUE` | 逻辑值 | 是否显示序列名；字符值可覆盖名称。 |
| `length` | `TRUE` | 逻辑值 | 是否显示格式化后的序列长度。 |
| `separator` | `"\n"` | 字符/枚举；取值见说明 | 名称与长度的分隔符。 |
| `name_style` | `list()` | 命名列表 | 名称和长度 grob 的固定文字样式命名列表，支持 colour、size、family、fontface、lineheight。 |
| `length_style` | `list()` | 命名列表 | 名称和长度 grob 的固定文字样式命名列表，支持 colour、size、family、fontface、lineheight。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_text() 的其他参数。 |

</details>

<span id="api-geom_seq_label"></span>
<details><summary><code>geom_seq_label()</code></summary>


[完整参考](reference/geom_seq_label.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 默认 NULL，使用预计算数据。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 默认 NULL，从布局自动获取。 |
| `seq_label_radius` | `1` | 数值/整数；范围见说明 | 可选数值/向量，标签径向位置相对序列弧半径的倍数：1 在弧上，> 1 在外侧（远离弦图中心），< 1 在内侧，默认 1。不能与 coord_collinear() 同时显式提供；该布局将序列名放在水平轨道左端旁。 |
| `seq_label_rotation` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值/向量，在沿弧方向上追加的旋转角度，默认 NULL（0 度）；seq_label_orientation = "horizontal" 时忽略。 |
| `labels` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选字符向量，覆盖名称文字，默认使用 geom_seq() 中的标签或序列 ID。 |
| `seq_label_orientation` | `c("arc", "horizontal")` | 向量；候选/多选含义见说明 | 文字方向，默认 "arc"。"arc" 沿序列弧旋转并保持可读，"horizontal" 水平绘制并向远离弦图中心的方向延伸。 |
| `seq_label_hjust` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值/向量，水平对齐，默认 NULL。沿弧方向默认 -0.2，使文字位于序列内侧附近；seq_label_orientation = "horizontal" 时自动选择，使文字向中心外延伸。 |
| `seq_label_vjust` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值/向量，垂直对齐，默认 NULL（0.5）。 |
| `check_overlap` | `FALSE` | 逻辑值 | 逻辑值，默认 FALSE。TRUE 时跳过与已绘标签重叠的文字，对应 ggplot2 的 geom_text() 选项。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_text() 的其他参数。 |

</details>

<span id="api-geom_seq_region"></span>
<details><summary><code>geom_seq_region()</code></summary>


[完整参考](reference/geom_seq_region.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 默认 NULL，使用预计算数据。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 数据框，至少含 accver、start、end，可选 label、category、color。因子 start/end 使用数值标签而非级别编码。 |
| `fill` | `"#F59E0B"` | 字符/枚举；取值见说明 | 区域默认填色，默认 "#F59E0B"。 |
| `colour` | `"#B45309"` | 字符/枚举；取值见说明 | 轮廓颜色，默认 "#B45309"。 |
| `alpha` | `0.25` | 数值/整数；范围见说明 | 0–1 数值，区域透明度，默认 0.25。 |
| `region_width` | `0.08` | 数值/整数；范围见说明 | 区域带宽，单位为弦图半径，默认 0.08。 |
| `region_offset` | `0` | 数值/整数；范围见说明 | 相对序列弧的径向偏移，默认 0。 |
| `region_side` | `c("inside", "outside", "auto")` | 向量；候选/多选含义见说明 | 侧向取 "inside"、"outside" 或 "auto"；"auto" 在可能时放在弦图内侧。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_polygon() 的其他参数。 |

</details>

<span id="api-geom_seq_scale"></span>
<details><summary><code>geom_seq_scale()</code></summary>


[完整参考](reference/geom_seq_scale.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入。可用 scale_length 映射每个尺的正整数长度，单位由序列表声明，label 提供文字。多个尺须在 data 中显式提供 x、y 布局位置或相应角色映射，避免意外重叠。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入。可用 scale_length 映射每个尺的正整数长度，单位由序列表声明，label 提供文字。多个尺须在 data 中显式提供 x、y 布局位置或相应角色映射，避免意外重叠。 |
| `length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选固定正整数长度，使用序列表的单位，不可与 scale_length 映射并用。NULL 使用输入长度或自动参照长度。 |
| `label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选固定字符标签，或接收所选长度的函数。NULL 使用输入标签或带 bp/aa 单位的格式化文字。 |
| `location` | `c("bottomleft", "bottomright", "topleft", "topright")` | 向量；候选/多选含义见说明 | 自动放置所依据的序列外包区域角落。 |
| `offset` | `grid::unit(c(0, 6), "mm")` | grid 物理单位对象 | 两个非负距离，先水平后垂直。支持绝对 grid 单位；普通数值为布局距离。默认水平 0 mm，垂直 6 mm。 |
| `tick_height` | `grid::unit(1.5, "mm")` | grid 物理单位对象 | 端刻度及尺与文字间距，支持非负绝对 grid 单位或布局数值。有渲染环境时做物理换算；无设备构建保留明确的估计测量状态。 |
| `label_gap` | `grid::unit(1.5, "mm")` | grid 物理单位对象 | 端刻度及尺与文字间距，支持非负绝对 grid 单位或布局数值。有渲染环境时做物理换算；无设备构建保留明确的估计测量状态。 |
| `linewidth` | `0.4` | 数值/整数；范围见说明 | 固定线色、线宽及字号，文字使用 ggplot2 的 mm。未作为固定值提供时，也可使用普通映射和标度。 |
| `colour` | `"#3D4348"` | 字符/枚举；取值见说明 | 固定线色、线宽及字号，文字使用 ggplot2 的 mm。未作为固定值提供时，也可使用普通映射和标度。 |
| `size` | `3` | 数值/整数；范围见说明 | 固定线色、线宽及字号，文字使用 ggplot2 的 mm。未作为固定值提供时，也可使用普通映射和标度。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::layer() 的其他固定 aesthetic。 |

</details>

<span id="api-geom_track"></span>
<details><summary><code>geom_track()</code></summary>


[完整参考](reference/geom_track.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入。提供 accver，以及 x/pos 源位置或 start/end 区间；给定数值映射到 y/track_value。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入。提供 accver，以及 x/pos 源位置或 start/end 区间；给定数值映射到 y/track_value。 |
| `geom` | `"point"` | 字符/枚举；取值见说明 | 原生几何与统计名称，bar 默认 count，col 及其他几何默认 identity。 |
| `stat` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 原生几何与统计名称，bar 默认 count，col 及其他几何默认 identity。 |
| `track` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 物理 position_track() 声明。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | identity、stack、dodge 等数据 Position；旧 position=position_track(...) 写法仍声明物理 band。 |
| `bounds` | `c("error", "clip")` | 向量；候选/多选含义见说明 | 坐标报错或裁剪策略。 |
| `baseline` | `0` | 数值/整数；范围见说明 | 给定数值柱的基线，非零值要求 position='identity'，避免堆叠改变其含义。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准图层控制参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准图层控制参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 原生 geom/stat 参数和固定样式。 |

</details>

<span id="api-geom_track_bar"></span>
<details><summary><code>geom_track_bar()</code></summary>


[完整参考](reference/geom_track.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入。提供 accver，以及 x/pos 源位置或 start/end 区间；给定数值映射到 y/track_value。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入。提供 accver，以及 x/pos 源位置或 start/end 区间；给定数值映射到 y/track_value。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | identity、stack、dodge 等数据 Position；旧 position=position_track(...) 写法仍声明物理 band。 |
| `baseline` | `0` | 数值/整数；范围见说明 | 给定数值柱的基线，非零值要求 position='identity'，避免堆叠改变其含义。 |
| `fill` | `"#86AFB9"` | 字符/枚举；取值见说明 | 计数快捷层的固定样式。 |
| `colour` | `NA` | R 对象；约束见说明 | 计数快捷层的固定样式。 |
| `alpha` | `1` | 数值/整数；范围见说明 | 计数快捷层的固定样式。 |
| `linewidth` | `0.25` | 数值/整数；范围见说明 | 计数快捷层的固定样式。 |
| `linetype` | `1` | 数值/整数；范围见说明 | 计数快捷层的固定样式。 |
| `bounds` | `c("error", "clip")` | 向量；候选/多选含义见说明 | 坐标报错或裁剪策略。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准图层控制参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准图层控制参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 原生 geom/stat 参数和固定样式。 |

</details>

<span id="api-geom_track_col"></span>
<details><summary><code>geom_track_col()</code></summary>


[完整参考](reference/geom_track.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层输入。提供 accver，以及 x/pos 源位置或 start/end 区间；给定数值映射到 y/track_value。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层输入。提供 accver，以及 x/pos 源位置或 start/end 区间；给定数值映射到 y/track_value。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | identity、stack、dodge 等数据 Position；旧 position=position_track(...) 写法仍声明物理 band。 |
| `track` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 物理 position_track() 声明。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 原生 geom/stat 参数和固定样式。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准图层控制参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准图层控制参数。 |

</details>

<span id="api-geom_track_line"></span>
<details><summary><code>geom_track_line()</code></summary>


[完整参考](reference/geom_track_tile.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准 ggplot2 映射与图层输入。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准 ggplot2 映射与图层输入。 |
| `position` | `必需/互斥入口，见说明` | Position 对象或已记录的字符串 | 每层必须提供 position_track() 声明。 |
| `colour` | `"#387D88"` | 字符/枚举；取值见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `alpha` | `1` | 数值/整数；范围见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `linewidth` | `0.4` | 数值/整数；范围见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `linetype` | `1` | 数值/整数；范围见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `closed` | `FALSE` | 逻辑值 | 是否显式闭合完整圆形折线。 |
| `bounds` | `c("error", "clip")` | 向量；候选/多选含义见说明 | 基因组坐标策略为 "error"（默认）或 "clip"。裁剪会报告，不改变源值。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他几何参数。 |

</details>

<span id="api-geom_track_tile"></span>
<details><summary><code>geom_track_tile()</code></summary>


[完整参考](reference/geom_track_tile.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准 ggplot2 映射与图层输入。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准 ggplot2 映射与图层输入。 |
| `position` | `必需/互斥入口，见说明` | Position 对象或已记录的字符串 | 每层必须提供 position_track() 声明。 |
| `fill` | `"#86AFB9"` | 字符/枚举；取值见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `colour` | `NA` | R 对象；约束见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `alpha` | `1` | 数值/整数；范围见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `linewidth` | `0.25` | 数值/整数；范围见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `linetype` | `1` | 数值/整数；范围见说明 | 固定样式；填色和颜色映射使用专用轨道角色及作用域标度。 |
| `bounds` | `c("error", "clip")` | 向量；候选/多选含义见说明 | 基因组坐标策略为 "error"（默认）或 "clip"。裁剪会报告，不改变源值。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他几何参数。 |

</details>


## 局部放置 {#positions}

<span id="api-position_feature_stack"></span>
<details><summary><code>position_feature_stack()</code></summary>


[完整参考](reference/position_feature_stack.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `spacing` | `0.1` | 数值/整数；范围见说明 | 相邻注释轨道的最小正径向距离。圆形质粒布局会加宽首个过渡作为主要标签通道，在更深重叠轨道间至少保留一条标签轨道。 |
| `side` | `c("strand", "outside", "inside")` | 向量；候选/多选含义见说明 | 为兼容保留的历史侧向选项，用于 base_position 为 NULL 时；历史名称保持原几何含义。 |
| `base_position` | `NULL` | Position 对象或已记录的字符串 | 可选 identity、strand 或 plasmid Position，在重叠检测和轨道分配前应用。 |

</details>

<span id="api-position_plasmid"></span>
<details><summary><code>position_plasmid()</code></summary>


[完整参考](reference/position_plasmid.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `offset` | `-0.1` | 数值/整数；范围见说明 | 有符号局部法线偏移，负值在圆形序列内侧，正值在外侧。 |

</details>

<span id="api-position_strand"></span>
<details><summary><code>position_strand()</code></summary>


[完整参考](reference/position_strand.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `offset` | `0.1` | 数值/整数；范围见说明 | 有限单值或灵活的逐序列/逐链参数。 |

</details>

<span id="api-position_track"></span>
<details><summary><code>position_track()</code></summary>


[完整参考](reference/position_track.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `track` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 一个稳定非空字符 ID。 |
| `side` | `c("outside", "inside")` | 向量；候选/多选含义见说明 | 支持 "outside" 或 "inside"。 |
| `order` | `1` | 数值/整数；范围见说明 | 有限排序优先级，靠近骨架的带优先；相同值按 ID 排序。 |
| `width` | `grid::unit(4, "mm")` | grid 物理单位对象 | 正带宽和非负前置间距。绝对 grid 单位（mm、cm、英寸、点）根据渲染环境换算；数值保留原生布局单位，逐序列数值向量/列表遵循常用序列参数格式。 |
| `gap` | `grid::unit(1, "mm")` | grid 物理单位对象 | 正带宽和非负前置间距。绝对 grid 单位（mm、cm、英寸、点）根据渲染环境换算；数值保留原生布局单位，逐序列数值向量/列表遵循常用序列参数格式。 |

</details>


## 连接统计 {#statistics}

<span id="api-stat_ribbon_bundle"></span>
<details><summary><code>stat_ribbon_bundle()</code></summary>


[完整参考](reference/stat_ribbon_bundle.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层映射与数据参数。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层映射与数据参数。 |
| `bins` | `80L` | 数值/整数；范围见说明 | 传给 bundle_ggchord_ribbons()。 |
| `min_bundle` | `2L` | 数值/整数；范围见说明 | 传给 bundle_ggchord_ribbons()。 |
| `weight` | `c("length", "count", "pident")` | 向量；候选/多选含义见说明 | 传给 bundle_ggchord_ribbons()。 |
| `group_by` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 传给 bundle_ggchord_ribbons()。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `TRUE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_link_ribbon() 的连接几何与外观参数。 |

</details>

<span id="api-stat_ribbon_density"></span>
<details><summary><code>stat_ribbon_density()</code></summary>


[完整参考](reference/stat_ribbon_density.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `mapping` | `NULL` | aes 映射 | 标准图层映射与数据参数。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数 | 标准图层映射与数据参数。 |
| `bins` | `80L` | 数值/整数；范围见说明 | 每条序列归一化中点分箱的正整数数量。 |
| `weight` | `c("length", "count", "pident")` | 向量；候选/多选含义见说明 | 密度权重："length"、"count" 或 "pident"。 |
| `group_by` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选列，定义独立密度组。 |
| `position` | `"identity"` | 字符/枚举；取值见说明 | 标准 ggplot2 图层参数。 |
| `show.legend` | `TRUE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `inherit.aes` | `FALSE` | 逻辑值 | 标准 ggplot2 图层参数。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 geom_link_ribbon() 的连接几何与外观参数。 |

</details>


## 角色标度 {#scales}

<span id="api-scale_feature_fill_manual"></span>
<details><summary><code>scale_feature_fill_manual()</code></summary>


[完整参考](reference/scale_gene_fill_manual.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::scale_fill_manual() 的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 将数据值对应到填色的一组值。 |
| `limits` | `NULL` | 向量/标度规则；按标度类型解释 | 可选显示顺序，命名注释值默认按名称顺序。 |

</details>

<span id="api-scale_feature_fill_plasmid"></span>
<details><summary><code>scale_feature_fill_plasmid()</code></summary>


[完整参考](reference/scale_feature_fill_plasmid.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他标度参数。 |
| `limits` | `NULL` | 向量/标度规则；按标度类型解释 | 可选注释类型顺序。 |

</details>

<span id="api-scale_feature_shape_manual"></span>
<details><summary><code>scale_feature_shape_manual()</code></summary>


[完整参考](reference/scale_feature_shape_manual.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::discrete_scale() 的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 注释几何名称的命名或无名字符向量。 |
| `name` | `"Feature"` | 字符/枚举；取值见说明 | 标度及图例标题，默认 "Feature"。 |
| `limits` | `NULL` | 向量/标度规则；按标度类型解释 | 可选类别顺序，命名 values 默认按名称顺序，以便形状与填色图例合并。 |
| `guide` | `"none"` | 字符/枚举；取值见说明 | 默认 "none"，因为 geom_feature() 在两种 aesthetic 表达相同类别时将形状轮廓合并到填色图例；提供 guide_ggchord_legend() 可单独显示形状图例。 |

</details>

<span id="api-scale_feature_shape_plasmid"></span>
<details><summary><code>scale_feature_shape_plasmid()</code></summary>


[完整参考](reference/scale_feature_fill_plasmid.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他标度参数。 |
| `limits` | `NULL` | 向量/标度规则；按标度类型解释 | 可选注释类型顺序。 |
| `guide` | `"none"` | 字符/枚举；取值见说明 | 形状预设使用的图例。 |

</details>

<span id="api-scale_gene_fill_manual"></span>
<details><summary><code>scale_gene_fill_manual()</code></summary>


[完整参考](reference/scale_gene_fill_manual.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::scale_fill_manual() 的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 将数据值对应到填色的一组值。 |

</details>

<span id="api-scale_link_alpha_continuous"></span>
<details><summary><code>scale_link_alpha_continuous()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |

</details>

<span id="api-scale_link_alpha_manual"></span>
<details><summary><code>scale_link_alpha_manual()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_link_color_gradientn"></span>
<details><summary><code>scale_link_color_gradientn()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |

</details>

<span id="api-scale_link_color_manual"></span>
<details><summary><code>scale_link_color_manual()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_link_colour_gradientn"></span>
<details><summary><code>scale_link_colour_gradientn()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |

</details>

<span id="api-scale_link_colour_manual"></span>
<details><summary><code>scale_link_colour_manual()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_link_linetype_manual"></span>
<details><summary><code>scale_link_linetype_manual()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_link_linewidth_continuous"></span>
<details><summary><code>scale_link_linewidth_continuous()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |

</details>

<span id="api-scale_link_linewidth_manual"></span>
<details><summary><code>scale_link_linewidth_manual()</code></summary>


[完整参考](reference/scale_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_region_fill_manual"></span>
<details><summary><code>scale_region_fill_manual()</code></summary>


[完整参考](reference/scale_gene_fill_manual.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::scale_fill_manual() 的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 将数据值对应到填色的一组值。 |

</details>

<span id="api-scale_ribbon_alpha_continuous"></span>
<details><summary><code>scale_ribbon_alpha_continuous()</code></summary>


[完整参考](reference/scale_ribbon_alpha_continuous.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `range` | `NULL` | 向量/标度规则；按标度类型解释 | 连续标度的输出透明度范围；NULL 使用对应 ggplot2 默认范围。 |

</details>

<span id="api-scale_ribbon_alpha_manual"></span>
<details><summary><code>scale_ribbon_alpha_manual()</code></summary>


[完整参考](reference/scale_ribbon_alpha_continuous.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_ribbon_color_manual"></span>
<details><summary><code>scale_ribbon_color_manual()</code></summary>


[完整参考](reference/scale_ribbon_alpha_continuous.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_ribbon_colour_manual"></span>
<details><summary><code>scale_ribbon_colour_manual()</code></summary>


[完整参考](reference/scale_ribbon_alpha_continuous.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_ribbon_fill_gradientn"></span>
<details><summary><code>scale_ribbon_fill_gradientn()</code></summary>


[完整参考](reference/scale_ribbon_fill_stepsn.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `colours` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 渐变颜色。 |
| `values` | `NULL` | 向量/标度规则；按标度类型解释 | 手动值；对应 ggplot2 标度支持时，也可为渐变颜色位置。 |
| `colors` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 渐变颜色。 |
| `guide` | `ggplot2::guide_colourbar(available_aes = "ribbon_fill")` | guide 对象或名称 | 图例对象，默认接受 ribbon_fill aesthetic。 |

</details>

<span id="api-scale_ribbon_fill_identity"></span>
<details><summary><code>scale_ribbon_fill_identity()</code></summary>


[完整参考](reference/scale_ribbon_fill_stepsn.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |

</details>

<span id="api-scale_ribbon_fill_manual"></span>
<details><summary><code>scale_ribbon_fill_manual()</code></summary>


[完整参考](reference/scale_ribbon_fill_stepsn.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动值；对应 ggplot2 标度支持时，也可为渐变颜色位置。 |

</details>

<span id="api-scale_ribbon_fill_stepsn"></span>
<details><summary><code>scale_ribbon_fill_stepsn()</code></summary>


[完整参考](reference/scale_ribbon_fill_stepsn.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `colours` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 渐变颜色。 |
| `values` | `NULL` | 向量/标度规则；按标度类型解释 | 手动值；对应 ggplot2 标度支持时，也可为渐变颜色位置。 |
| `colors` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 渐变颜色。 |
| `guide` | `ggplot2::guide_coloursteps(available_aes = "ribbon_fill")` | guide 对象或名称 | 图例对象，默认接受 ribbon_fill aesthetic。 |

</details>

<span id="api-scale_ribbon_linetype_manual"></span>
<details><summary><code>scale_ribbon_linetype_manual()</code></summary>


[完整参考](reference/scale_ribbon_alpha_continuous.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 标度的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动 aesthetic 值。 |

</details>

<span id="api-scale_seq_color_manual"></span>
<details><summary><code>scale_seq_color_manual()</code></summary>


[完整参考](reference/scale_seq_colour_manual.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::scale_colour_manual() 的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 将数据值对应到 aesthetic 的一组值。 |

</details>

<span id="api-scale_seq_colour_manual"></span>
<details><summary><code>scale_seq_colour_manual()</code></summary>


[完整参考](reference/scale_seq_colour_manual.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::scale_colour_manual() 的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 将数据值对应到 aesthetic 的一组值。 |

</details>

<span id="api-scale_seq_position_continuous"></span>
<details><summary><code>scale_seq_position_continuous()</code></summary>


[完整参考](reference/scale_seq_position_continuous.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `name` | `ggplot2::waiver()` | R 对象；约束见说明 | 标度名称，默认禁用位置图例。 |
| `breaks` | `ggplot2::waiver()` | 向量/标度规则；按标度类型解释 | 标准连续标度控制项。 |
| `minor_breaks` | `ggplot2::waiver()` | 向量/标度规则；按标度类型解释 | 标准连续标度控制项。 |
| `labels` | `ggplot2::waiver()` | R 对象；约束见说明 | 标准连续标度控制项。 |
| `limits` | `NULL` | 向量/标度规则；按标度类型解释 | 标准连续标度控制项。 |
| `expand` | `ggplot2::waiver()` | R 对象；约束见说明 | 标准连续标度控制项。 |
| `oob` | `scales::censor` | R 对象；约束见说明 | 标准连续标度控制项。 |
| `transform` | `"identity"` | 字符/枚举；取值见说明 | 标准连续标度控制项。 |
| `n.breaks` | `6` | 数值/整数；范围见说明 | 标准连续标度控制项。 |
| `coordinates` | `c("relative", "source")` | 向量；候选/多选含义见说明 | 在线性图中使用局部相对坐标，或 focus_ggchord_data() 保留的源坐标。只影响轴刻度和标签，不改变几何。所有模式均将变换空间的刻度反变换回基因组位置；显式 breaks = NULL、minor_breaks = NULL 分别关闭对应刻度。 |
| `sequences` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选字符向量，以显示 accver 选择序列轴。NULL 保留全部，character() 隐藏全部。局部/分段视图使用窗口 ID，不是父源 accession；未知或未显示 ID 会报错。同时控制轴线、刻度和文字，不删除序列几何。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::continuous_scale() 的其他参数。 |

</details>

<span id="api-scale_seq_ring_manual"></span>
<details><summary><code>scale_seq_ring_manual()</code></summary>


[完整参考](reference/scale_seq_ring_manual.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 ggplot2::discrete_scale() 的参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 正半径的命名或无名数值向量。 |
| `name` | `"Ring"` | 字符/枚举；取值见说明 | 标度标题。默认隐藏图例，因为该标度控制几何而非可见图例属性。 |
| `limits` | `NULL` | 向量/标度规则；按标度类型解释 | 可选环顺序。 |
| `guide` | `"none"` | 字符/枚举；取值见说明 | 图例定义，默认 "none"。 |

</details>

<span id="api-scale_track_colour_gradientn"></span>
<details><summary><code>scale_track_colour_gradientn()</code></summary>


[完整参考](reference/scale_track_fill_gradientn.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `track` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 轨道 ID。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 标度参数。 |
| `colours` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 渐变颜色。 |
| `guide` | `NULL` | guide 对象或名称 | 可选标准填色/颜色图例。 |

</details>

<span id="api-scale_track_colour_manual"></span>
<details><summary><code>scale_track_colour_manual()</code></summary>


[完整参考](reference/scale_track_fill_gradientn.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `track` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 轨道 ID。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 标度参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动类别颜色。 |

</details>

<span id="api-scale_track_fill_gradientn"></span>
<details><summary><code>scale_track_fill_gradientn()</code></summary>


[完整参考](reference/scale_track_fill_gradientn.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `track` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 轨道 ID。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 标度参数。 |
| `colours` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 渐变颜色。 |
| `guide` | `NULL` | guide 对象或名称 | 可选标准填色/颜色图例。 |

</details>

<span id="api-scale_track_fill_manual"></span>
<details><summary><code>scale_track_fill_manual()</code></summary>


[完整参考](reference/scale_track_fill_gradientn.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `track` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 轨道 ID。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 标度参数。 |
| `values` | `必需/互斥入口，见说明` | 向量/标度规则；按标度类型解释 | 手动类别颜色。 |

</details>

<span id="api-scale_track_value_continuous"></span>
<details><summary><code>scale_track_value_continuous()</code></summary>


[完整参考](reference/scale_track_value_continuous.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `track` | `必需/互斥入口，见说明` | R 对象；约束见说明 | position_track() 中的轨道 ID。 |
| `name` | `track` | R 对象；约束见说明 | 轴标题，默认轨道 ID。 |
| `limits` | `NULL` | 向量/标度规则；按标度类型解释 | 标准连续标度控制项。 |
| `breaks` | `ggplot2::waiver()` | 向量/标度规则；按标度类型解释 | 标准连续标度控制项。 |
| `labels` | `ggplot2::waiver()` | R 对象；约束见说明 | 标准连续标度控制项。 |
| `transform` | `"identity"` | 字符/枚举；取值见说明 | 标准连续标度控制项。 |
| `oob` | `scales::censor` | R 对象；约束见说明 | 标准连续标度控制项。 |
| `unit` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数据单位标签，保留在轨道导出中。 |
| `guide` | `c("axis", "none")` | 向量；候选/多选含义见说明 | 数值轨道轴使用 "axis" 或 "none"。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 其他 ggplot2 连续标度参数。 |

</details>


## 图例与 guide {#guides}

<span id="api-guide_ggchord_colorbar"></span>
<details><summary><code>guide_ggchord_colorbar()</code></summary>


[完整参考](reference/guide_ggchord_legend.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `title` | `ggplot2::waiver()` | R 对象；约束见说明 | 图例标题。 |
| `theme` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选图例专属主题。 |
| `nbin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 颜色条分箱数。 |
| `display` | `"raster"` | 字符/枚举；取值见说明 | 颜色条显示方式。 |
| `alpha` | `NA` | R 对象；约束见说明 | 颜色条透明度覆盖值。 |
| `draw.ulim` | `TRUE` | 逻辑值 | 是否绘制上/下限刻度。 |
| `draw.llim` | `TRUE` | 逻辑值 | 是否绘制上/下限刻度。 |
| `angle` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 刻度标签角度。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | 图例位置和方向。 |
| `direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 图例位置和方向。 |
| `reverse` | `FALSE` | 逻辑值 | 是否反转图例键顺序。 |
| `order` | `0` | 数值/整数；范围见说明 | 图例顺序。 |
| `available_aes` | `c("colour", "color", "fill", "ribbon_fill", "gene_fill", "feature_fill",      "region_fill")` | 向量；候选/多选含义见说明 | 图例接受的 aesthetic。 |
| `size_scale` | `1` | 数值/整数；范围见说明 | 图例键、间距及文字的正缩放系数。自动图例由当前导出设备推导，直接调用默认 1。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 图例的其他参数。 |

</details>

<span id="api-guide_ggchord_colourbar"></span>
<details><summary><code>guide_ggchord_colourbar()</code></summary>


[完整参考](reference/guide_ggchord_legend.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `title` | `ggplot2::waiver()` | R 对象；约束见说明 | 图例标题。 |
| `theme` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选图例专属主题。 |
| `nbin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 颜色条分箱数。 |
| `display` | `"raster"` | 字符/枚举；取值见说明 | 颜色条显示方式。 |
| `alpha` | `NA` | R 对象；约束见说明 | 颜色条透明度覆盖值。 |
| `draw.ulim` | `TRUE` | 逻辑值 | 是否绘制上/下限刻度。 |
| `draw.llim` | `TRUE` | 逻辑值 | 是否绘制上/下限刻度。 |
| `angle` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 刻度标签角度。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | 图例位置和方向。 |
| `direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 图例位置和方向。 |
| `reverse` | `FALSE` | 逻辑值 | 是否反转图例键顺序。 |
| `order` | `0` | 数值/整数；范围见说明 | 图例顺序。 |
| `available_aes` | `c("colour", "color", "fill", "ribbon_fill", "gene_fill", "feature_fill",      "region_fill")` | 向量；候选/多选含义见说明 | 图例接受的 aesthetic。 |
| `size_scale` | `1` | 数值/整数；范围见说明 | 图例键、间距及文字的正缩放系数。自动图例由当前导出设备推导，直接调用默认 1。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 图例的其他参数。 |

</details>

<span id="api-guide_ggchord_legend"></span>
<details><summary><code>guide_ggchord_legend()</code></summary>


[完整参考](reference/guide_ggchord_legend.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `title` | `ggplot2::waiver()` | R 对象；约束见说明 | 图例标题。 |
| `theme` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选图例专属主题。 |
| `position` | `NULL` | Position 对象或已记录的字符串 | 图例位置和方向。 |
| `direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 图例位置和方向。 |
| `override.aes` | `list()` | 命名列表 | 覆盖图例键 aesthetic 的列表。 |
| `nrow` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 图例键排列的行列数。 |
| `ncol` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 图例键排列的行列数。 |
| `reverse` | `FALSE` | 逻辑值 | 是否反转图例键顺序。 |
| `order` | `0` | 数值/整数；范围见说明 | 图例顺序。 |
| `size_scale` | `1` | 数值/整数；范围见说明 | 图例键、间距及文字的正缩放系数。自动图例由当前导出设备推导，直接调用默认 1。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给对应 ggplot2 图例的其他参数。 |

</details>


## 主题与物理外观 {#themes}

<span id="api-theme_ggchord"></span>
<details><summary><code>theme_ggchord()</code></summary>


[完整参考](reference/theme_ggchord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `axis` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示自动序列轴，ggplot2::element_blank() 隐藏整组轴装饰。 |
| `axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.minor.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴文字元素。 |
| `axis.gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.minor.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.offset` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.orientation` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 取"horizontal"、"parallel"、"perpendicular"或有限数值角度。圆形布局的 parallel 文字逐字符沿圆弧排列；水平和数值角度保留直线文字。 |
| `axis.text.check.overlap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 是否隐藏重叠坐标轴标签，默认启用。 |
| `seq.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `feature.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `feature.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `restriction.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `restriction.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `seq.center.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.name` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `legend.seq` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.seq.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.seq.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.seq.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.seq.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.seq.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.seq.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.seq.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.seq.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.seq.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.seq.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.link` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.ribbon` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.ribbon.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.ribbon.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.ribbon.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.ribbon.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.ribbon.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.ribbon.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.ribbon.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.ribbon.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.ribbon.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.ribbon.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.ribbon.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.ribbon.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条刻度长度。 |
| `legend.ribbon.axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.gene` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.gene.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.gene.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.gene.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.gene.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.gene.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.gene.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.gene.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.gene.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.gene.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.gene.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.feature` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.feature.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.feature.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.feature.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.feature.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.feature.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.feature.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.feature.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.feature.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.feature.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.feature.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.region` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.region.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.region.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.region.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.region.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.region.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.region.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.region.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.region.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.region.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.region.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 主题参数，包括公共文字、图及图例设置；公共图例设置也提供角色继承。已删除的下划线风格参数会给出迁移错误。 |

</details>

<span id="api-theme_ggchord_dark"></span>
<details><summary><code>theme_ggchord_dark()</code></summary>


[完整参考](reference/theme_ggchord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `axis` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示自动序列轴，ggplot2::element_blank() 隐藏整组轴装饰。 |
| `axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.minor.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴文字元素。 |
| `axis.gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.minor.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.offset` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.orientation` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 取"horizontal"、"parallel"、"perpendicular"或有限数值角度。圆形布局的 parallel 文字逐字符沿圆弧排列；水平和数值角度保留直线文字。 |
| `axis.text.check.overlap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 是否隐藏重叠坐标轴标签，默认启用。 |
| `seq.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `feature.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `feature.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `restriction.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `restriction.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `seq.center.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.name` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `legend.seq` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.seq.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.seq.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.seq.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.seq.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.seq.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.seq.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.seq.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.seq.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.seq.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.seq.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.link` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.ribbon` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.ribbon.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.ribbon.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.ribbon.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.ribbon.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.ribbon.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.ribbon.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.ribbon.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.ribbon.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.ribbon.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.ribbon.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.ribbon.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.ribbon.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条刻度长度。 |
| `legend.ribbon.axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.gene` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.gene.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.gene.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.gene.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.gene.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.gene.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.gene.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.gene.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.gene.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.gene.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.gene.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.feature` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.feature.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.feature.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.feature.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.feature.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.feature.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.feature.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.feature.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.feature.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.feature.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.feature.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.region` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.region.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.region.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.region.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.region.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.region.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.region.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.region.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.region.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.region.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.region.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 主题参数，包括公共文字、图及图例设置；公共图例设置也提供角色继承。已删除的下划线风格参数会给出迁移错误。 |

</details>

<span id="api-theme_ggchord_minimal"></span>
<details><summary><code>theme_ggchord_minimal()</code></summary>


[完整参考](reference/theme_ggchord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `axis` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示自动序列轴，ggplot2::element_blank() 隐藏整组轴装饰。 |
| `axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.minor.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴文字元素。 |
| `axis.gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.minor.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.offset` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.orientation` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 取"horizontal"、"parallel"、"perpendicular"或有限数值角度。圆形布局的 parallel 文字逐字符沿圆弧排列；水平和数值角度保留直线文字。 |
| `axis.text.check.overlap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 是否隐藏重叠坐标轴标签，默认启用。 |
| `seq.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `feature.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `feature.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `restriction.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `restriction.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `seq.center.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.name` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `legend.seq` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.seq.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.seq.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.seq.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.seq.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.seq.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.seq.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.seq.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.seq.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.seq.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.seq.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.link` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.ribbon` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.ribbon.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.ribbon.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.ribbon.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.ribbon.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.ribbon.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.ribbon.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.ribbon.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.ribbon.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.ribbon.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.ribbon.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.ribbon.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.ribbon.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条刻度长度。 |
| `legend.ribbon.axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.gene` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.gene.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.gene.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.gene.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.gene.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.gene.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.gene.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.gene.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.gene.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.gene.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.gene.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.feature` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.feature.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.feature.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.feature.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.feature.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.feature.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.feature.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.feature.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.feature.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.feature.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.feature.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.region` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.region.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.region.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.region.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.region.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.region.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.region.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.region.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.region.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.region.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.region.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 主题参数，包括公共文字、图及图例设置；公共图例设置也提供角色继承。已删除的下划线风格参数会给出迁移错误。 |

</details>

<span id="api-theme_ggchord_plasmid"></span>
<details><summary><code>theme_ggchord_plasmid()</code></summary>


[完整参考](reference/theme_ggchord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `axis` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示自动序列轴，ggplot2::element_blank() 隐藏整组轴装饰。 |
| `axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.minor.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴文字元素。 |
| `axis.gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.minor.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.offset` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.orientation` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 取"horizontal"、"parallel"、"perpendicular"或有限数值角度。圆形布局的 parallel 文字逐字符沿圆弧排列；水平和数值角度保留直线文字。 |
| `axis.text.check.overlap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 是否隐藏重叠坐标轴标签，默认启用。 |
| `seq.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `feature.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `feature.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `restriction.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `restriction.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `seq.center.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.name` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `legend.seq` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.seq.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.seq.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.seq.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.seq.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.seq.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.seq.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.seq.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.seq.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.seq.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.seq.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.link` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.ribbon` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.ribbon.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.ribbon.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.ribbon.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.ribbon.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.ribbon.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.ribbon.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.ribbon.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.ribbon.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.ribbon.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.ribbon.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.ribbon.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.ribbon.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条刻度长度。 |
| `legend.ribbon.axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.gene` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.gene.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.gene.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.gene.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.gene.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.gene.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.gene.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.gene.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.gene.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.gene.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.gene.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.feature` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.feature.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.feature.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.feature.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.feature.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.feature.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.feature.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.feature.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.feature.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.feature.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.feature.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.region` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.region.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.region.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.region.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.region.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.region.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.region.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.region.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.region.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.region.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.region.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 主题参数，包括公共文字、图及图例设置；公共图例设置也提供角色继承。已删除的下划线风格参数会给出迁移错误。 |

</details>

<span id="api-theme_ggchord_publication"></span>
<details><summary><code>theme_ggchord_publication()</code></summary>


[完整参考](reference/theme_ggchord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `axis` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示自动序列轴，ggplot2::element_blank() 隐藏整组轴装饰。 |
| `axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.minor.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴线元素。 |
| `axis.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 坐标轴文字元素。 |
| `axis.gap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.minor.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.offset` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 控制轴间距的单个绝对 grid::unit() 值。负 axis.gap 将圆图轴放在骨架内侧，正值放在外侧。零将刻度附着于骨架内缘，是圆图和质粒主题的默认值。圆图主/次刻度默认长 1.5/0.7 mm，质粒主题为 2.2/0.8 mm；对应默认线宽为 0.45/0.3 mm 和 0.6/0.35 mm。显式长度及线条样式优先。 |
| `axis.text.orientation` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 取"horizontal"、"parallel"、"perpendicular"或有限数值角度。圆形布局的 parallel 文字逐字符沿圆弧排列；水平和数值角度保留直线文字。 |
| `axis.text.check.overlap` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 是否隐藏重叠坐标轴标签，默认启用。 |
| `seq.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `gene.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `feature.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `feature.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `restriction.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `restriction.label.segment` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 基因、注释及酶切引导线的线元素。 |
| `seq.center.label` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.name` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `seq.center.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 序列、基因、通用注释、酶切及圆图中心标签的文字元素。 |
| `legend.seq` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.seq.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.seq.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.seq.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.seq.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.seq.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.seq.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.seq.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.seq.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.seq.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.seq.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.seq.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.link` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.link.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 点连接图例外观，遵循相同角色专属规则。 |
| `legend.ribbon` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.ribbon.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.ribbon.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.ribbon.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.ribbon.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.ribbon.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.ribbon.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.ribbon.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.ribbon.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.ribbon.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.ribbon.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.ribbon.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.ribbon.ticks` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.ribbon.ticks.length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条刻度长度。 |
| `legend.ribbon.axis.line` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 区间连接颜色条线元素。 |
| `legend.gene` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.gene.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.gene.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.gene.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.gene.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.gene.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.gene.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.gene.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.gene.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.gene.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.gene.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.gene.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.feature` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.feature.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.feature.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.feature.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.feature.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.feature.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.feature.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.feature.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.feature.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.feature.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.feature.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.feature.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `legend.region` | `NULL` | 可空 R 对象；具体类型/继承见说明 | NULL 显示并继承公共设置，ggplot2::element_blank() 隐藏该角色的全部图例。 |
| `legend.region.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例位置。 |
| `legend.region.direction` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例方向。 |
| `legend.region.background` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例背景。 |
| `legend.region.key` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键元素。 |
| `legend.region.key.size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.width` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.key.height` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例键的尺寸。自动 Identity 标尺独立测量边界文字，不把通用色块尺寸当成整条标尺长度。需要手动指定标尺长度时，使用 Ribbon 专用尺寸或显式 guide。 |
| `legend.region.text` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字元素。 |
| `legend.region.text.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例文字位置。 |
| `legend.region.title` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题元素。 |
| `legend.region.title.position` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例标题位置。 |
| `legend.region.margin` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例边距。 |
| `legend.region.spacing` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 各角色图例间距。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 标准 ggplot2 主题参数，包括公共文字、图及图例设置；公共图例设置也提供角色继承。已删除的下划线风格参数会给出迁移错误。 |

</details>


## 输入、转换、诊断与输出 {#workflows}

<span id="api-as_ggchord_alignment"></span>
<details><summary><code>as_ggchord_alignment()</code></summary>


[完整参考](reference/as_ggchord_alignment.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `sequences` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 已对齐字符串的命名字符向量，或含唯一 accver 和字符列 sequence 的数据框。所有字符串的正宽度必须相同，每条记录须至少含一个非 gap 字符。 |
| `alphabet` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 显式指定 "DNA"、"RNA" 或 "AA"。接受 IUPAC 核酸含糊码及蛋白含糊/稀有残基码；* 仅用于 AA。不猜测字母表。 |
| `source_data` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数据框，每个对齐 accver 对应一行，包含 start（首个非 gap 源位置）和 strand（"+" 或 "-"）。可选 source_accver 声明不同的源 ID，可选字符列 sequence 提供完整正向源字符串以供核对。坐标必须为正且精确的整数。输入元数据保留在 ggchord_conversion 属性中。 |

</details>

<span id="api-as_ggchord_feature"></span>
<details><summary><code>as_ggchord_feature()</code></summary>


[完整参考](reference/as_ggchord_feature.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `x` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 数据框或 GRanges 对象。 |
| `mapping` | `NULL` | aes 映射 | 可选的命名字符向量，将输出角色对应到源列，例如 as_ggchord_seq()。 |

</details>

<span id="api-as_ggchord_link"></span>
<details><summary><code>as_ggchord_link()</code></summary>


[完整参考](reference/as_ggchord_link.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `x` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 数据框或 GInteractions 对象。 |
| `mapping` | `NULL` | aes 映射 | 可选的命名字符向量，将标准连接角色对应到源列，例如 as_ggchord_seq()。 |

</details>

<span id="api-as_ggchord_seq"></span>
<details><summary><code>as_ggchord_seq()</code></summary>


[完整参考](reference/as_ggchord_seq.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `x` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 数据框、GRanges、Seqinfo、GInteractions，或 Biostrings 的 XStringSet/XString 对象。 |
| `mapping` | `NULL` | aes 映射 | 数据框输入可使用命名字符向量 output_role = "input_column"，例如 c(accver = "chrom", length = "size")。 |
| `lengths` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的已知序列长度命名数值向量，用于补足未知长度；必须与源中已知长度一致。 |
| `ids` | `NULL` | 字符路径、标识或向量；见约束 | 可选的显式标识字符向量，仅用于未命名的 XStringSet/XString 输入；必须与已有源名称一致。 |

</details>

<span id="api-build_ggchord_cds"></span>
<details><summary><code>build_ggchord_cds()</code></summary>


[完整参考](reference/build_ggchord_cds.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `cds_data` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 数据框，含 feature_id、accver、有序一基闭区间 start/end、strand（+/-）及 GFF3 phase，可为 0–2 的数值或 read_gff3() 返回的原文。. 保持未知；嵌套 CDS 初始 phase 未知时须显式提供 codon_start。存在 type 时应先明确选择 CDS 行。也可每个 CDS 一行，提供 segments 列表列及显式初始 phase 或 codon_start（1–3）。嵌套基因组片段保留显式 segment_order，不会被父级外包区间替代。 |
| `genetic_code` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选正整数遗传密码 ID，可为单值或按 CDS feature_id 命名，须与已有 genetic_code 列一致。缺失 ID 保持未知，但不妨碍只做坐标投影。 |
| `partial` | `c("error", "trim")` | 向量；候选/多选含义见说明 | 默认 "error"，或显式选择 "trim"。 |

</details>

<span id="api-bundle_ggchord_ribbons"></span>
<details><summary><code>bundle_ggchord_ribbons()</code></summary>


[完整参考](reference/bundle_ggchord_ribbons.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `ribbon_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | ggchord 连接格式的比对数据。 |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 含 accver、length 的序列数据。 |
| `bins` | `80L` | 数值/整数；范围见说明 | 每条序列归一化中点分箱的正整数数量，默认 80。 |
| `min_bundle` | `2L` | 数值/整数；范围见说明 | 聚合至少需要的行数，默认 2；较小组仍保留独立连接。 |
| `weight` | `c("length", "count", "pident")` | 向量；候选/多选含义见说明 | 聚束权重："length"、"count" 或 "pident"；相似度权重为比对长度乘以相似度比例。 |
| `group_by` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选分组列，值须一致才能聚束。 |

</details>

<span id="api-clean_ggchord_data"></span>
<details><summary><code>clean_ggchord_data()</code></summary>


[完整参考](reference/clean_ggchord_data.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 必需 data.frame/tibble，含 accver、length，作为裁剪坐标参照。不论 unknown_id 如何，参照 ID 均须非缺失、非空、唯一。 |
| `ribbon_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选比对 data.frame/tibble 或 NULL；其他类型会被拒绝，不视作缺省。 |
| `gene_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选基因注释 data.frame/tibble 或 NULL；其他类型会被拒绝，不视作缺省。 |
| `unknown_id` | `c("drop", "error", "keep")` | 向量；候选/多选含义见说明 | 未知 ID 策略，默认 "drop"。ribbon_data 或 gene_data 中 ID 缺失、为空或未知时，"drop" 删除，"error" 报错停止，"keep" 保持原样。保留但无参照序列的行会报告，不清理其坐标/注释，也不保证能绘制。不会猜测或替换 ID；有限数值坐标和有效基因链等结构检查仍执行。 |
| `out_of_range` | `c("clip", "drop", "error", "keep")` | 向量；候选/多选含义见说明 | 越界策略，默认 "clip"。坐标超出 [1, sequence length] 时，"clip" 限制到边界，"drop" 删除，"error" 停止，"keep" 保持。采用一基闭区间，start == end 是有效单碱基区间，会保留。"clip" 时，完全越界的区间删除并报告；部分重叠保留交集，即使仅一个碱基。 |
| `reversed_interval` | `c("sort", "drop", "error", "keep")` | 向量；候选/多选含义见说明 | 反向区间策略，默认 "sort"。对 start > end，"sort" 交换端点以稳定绘制，并在报告记录原方向；"drop" 删除，"error" 停止，"keep" 保持。 |
| `invalid_pident` | `c("clip", "drop", "error", "keep")` | 向量；候选/多选含义见说明 | 相似度策略，默认 "clip"。pident 超出 [0, 100] 时，"clip" 限制范围，"drop" 删除，"error" 停止，"keep" 保持。 |
| `empty_annotation` | `c("keep", "drop", "replace")` | 向量；候选/多选含义见说明 | 空注释策略，默认 "keep"。基因 anno 缺失或为空时，"replace" 用 replacement_annotation 填充，"drop" 删除，"keep" 保持。 |
| `replacement_annotation` | `"unannotated"` | 字符/枚举；取值见说明 | 替换文字，默认 "unannotated"，用于 empty_annotation = "replace"。 |
| `feature_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 统一注释输入，与 gene_data 互斥。 |

</details>

<span id="api-deduplicate_ggchord_ribbons"></span>
<details><summary><code>deduplicate_ggchord_ribbons()</code></summary>


[完整参考](reference/deduplicate_ggchord_ribbons.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `ribbon_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 连接格式的数据框。 |
| `tolerance` | `0` | 数值/整数；范围见说明 | 数值，默认 0。两行在 qstart/qend/sstart/send 各端点允许的最大绝对差，单位 bp，用于 by = "coordinates"。 |
| `by` | `c("exact", "coordinates", "overlap")` | 向量；候选/多选含义见说明 | 去重方式，默认 "exact"。"exact" 去除同序列对、同坐标；"coordinates" 另将 tolerance bp 内视为重复；"overlap" 要求两条序列上的互惠重叠均至少 min_reciprocal_overlap，每个区间都须达标，不仅较短者。所有方式保留相反相对方向；单碱基端点使用提供的链/方向证据，否则方向未知。 |
| `keep` | `c("best_pident", "longest", "first")` | 向量；候选/多选含义见说明 | 保留代表的规则，默认 "best_pident"；"best_pident" 取最高 pident，"longest" 取最长 length，"first" 取输入首次出现。 |
| `min_reciprocal_overlap` | `0.9` | 数值/整数；范围见说明 | 0–1 数值，默认 0.9，用于 by = "overlap" 的互惠重叠阈值。 |

</details>

<span id="api-diagnose_ggchord"></span>
<details><summary><code>diagnose_ggchord()</code></summary>


[完整参考](reference/diagnose_ggchord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `plot` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 一个 ggchord 图。 |
| `width` | `7` | 数值/整数；范围见说明 | 正目标画布尺寸。 |
| `height` | `5` | 数值/整数；范围见说明 | 正目标画布尺寸。 |
| `units` | `c("in", "cm", "mm", "px")` | 向量；候选/多选含义见说明 | 尺寸单位 "in"、"cm"、"mm" 或 "px"。 |
| `device` | `c("png", "pdf", "svg")` | 向量；候选/多选含义见说明 | 目标设备 "png"、"pdf" 或 "svg"，必须在当前 R 环境可用；仅检查此设备。 |
| `dpi` | `144` | 数值/整数；范围见说明 | 栅格分辨率，同时确定像素尺寸。 |
| `min_size` | `2` | 数值/整数；范围见说明 | 建议的最小字号及线宽，均为 mm，不修改输入样式。 |
| `min_linewidth` | `0.2` | 数值/整数；范围见说明 | 建议的最小字号及线宽，均为 mm，不修改输入样式。 |
| `colour_vision` | `FALSE` | 逻辑值 | 是否请求可选的色觉模拟提示。 |

</details>

<span id="api-export_ggchord_fusion"></span>
<details><summary><code>export_ggchord_fusion()</code></summary>


[完整参考](reference/export_ggchord_fusion.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `plot` | `必需/互斥入口，见说明` | R 对象；约束见说明 | Fusion对象。 |
| `width` | `10` | 数值/整数；范围见说明 | 以英寸指定的正画布宽高。 |
| `height` | `8` | 数值/整数；范围见说明 | 以英寸指定的正画布宽高。 |

</details>

<span id="api-export_ggchord_layout"></span>
<details><summary><code>export_ggchord_layout()</code></summary>


[完整参考](reference/export_ggchord_layout.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `plot` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 一个 ggchord 图。 |
| `include` | `c("seq", "ribbon", "gene", "feature", "labels", "axis", "track")` | 向量；候选/多选含义见说明 | 选择组件的字符向量："seq"、"ribbon"、"gene"、"feature"、"labels"、"alignment"、"link"、"restriction"、"track"及"axis"。点连接和酶切位点须显式选择。Dotplot默认选择alignment及axis；axis_layout保留可逆轴偏移、方向和源原点。比对导出保留四端点、图层身份和源行，省略原因写入connection_report。 |
| `original_data` | `FALSE` | 逻辑值 | 逻辑值，是否在 original_data 包含逐层解析输入表，默认 FALSE。 |

</details>

<span id="api-filter_ggchord_ribbons"></span>
<details><summary><code>filter_ggchord_ribbons()</code></summary>


[完整参考](reference/filter_ggchord_ribbons.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `ribbon_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 数据框，至少含 qaccver、saccver，通常为完整连接格式。 |
| `seq_ids` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选字符向量，只保留 query、subject 都位于 seq_ids 的行。 |
| `min_pident` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值，pident 的下限/上限。 |
| `max_pident` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值，pident 的下限/上限。 |
| `min_length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值，length 下限。 |
| `max_evalue` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值，evalue 上限，要求 evalue 列。 |
| `min_bitscore` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数值，bitscore 下限，要求 bitscore 列。 |
| `min_query_coverage` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选 0–100 数值，查询覆盖率下限。优先使用 qcovs，否则由 qlen 与查询区间计算，要求 qlen。 |
| `min_subject_coverage` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选 0–100 数值，目标覆盖率下限，由 slen 与目标区间计算，要求 slen 列。 |
| `keep_pairs` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数据框、列表或矩阵，指定无向序列对。数据框/矩阵前两列为 query/subject ID，列表各项为长度 2 的字符向量；任一方向匹配即保留。 |
| `drop_self_links` | `TRUE` | 逻辑值 | 逻辑值，默认 TRUE，是否删除 qaccver == saccver 的行。 |
| `sort_by` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选排序列名字符向量，用 "-" 或 "desc:" 前缀表示降序，例如 c("pident", "-evalue")。 |

</details>

<span id="api-filter_restriction_sites"></span>
<details><summary><code>filter_restriction_sites()</code></summary>


[完整参考](reference/filter_restriction_sites.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `sites` | `必需/互斥入口，见说明` | R 对象；约束见说明 | find_restriction_sites() 的结果。 |
| `set` | `c("all", "unique", "unique_dual", "six_plus", "unique_6plus",      "commercial", "reference")` | 向量；候选/多选含义见说明 | 显示预设。"reference" 复用精确识别参考序列保存的酶集合，没有该配置会报错；已保存的 "None" 配置返回空位点。 |
| `enzymes` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 额外筛选条件。 |
| `min_site_length` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 额外筛选条件。 |
| `cuts` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 额外筛选条件。 |
| `window` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 额外筛选条件。 |
| `commercial_only` | `FALSE` | 逻辑值 | 额外筛选条件。 |
| `parent_set` | `c("all", "commercial_all", "commercial_nonredundant")` | 向量；候选/多选含义见说明 | 酶目录子集。"commercial_nonredundant" 为等价识别组和位点保留稳定的单个商业代表，优先使用数据库首选酶。 |

</details>

<span id="api-find_common_features"></span>
<details><summary><code>find_common_features()</code></summary>


[完整参考](reference/find_common_features.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `sequence` | `必需/互斥入口，见说明` | 序列对象/字符串/完整序列表，见说明 | DNA 文本、命名字符向量，或含 accver、sequence 的数据框。 |
| `database` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选规范化数据库列表或自定义数据框。 |
| `types` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的注释类型及名称/ID 筛选。 |
| `features` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的注释类型及名称/ID 筛选。 |
| `circular` | `TRUE` | 逻辑值 | 是否跨序列原点搜索。 |
| `mode` | `c("auto", "dna", "protein")` | 向量；候选/多选含义见说明 | 匹配模式。"auto" 先为完全已知序列使用整理注释，再按数据库各记录的检测模式处理。标为 exactProteinMatch 的记录先做精确蛋白匹配，再做近精确 DNA 匹配；其余使用精确 DNA。"protein" 显式启用精确及近似蛋白匹配。 |
| `min_protein_identity` | `0.9` | 数值/整数；范围见说明 | 近似蛋白匹配阈值，各为 [0, 1] 内单个有限数。 |
| `min_protein_coverage` | `0.8` | 数值/整数；范围见说明 | 近似蛋白匹配阈值，各为 [0, 1] 内单个有限数。 |
| `resolve` | `c("best", "all")` | 向量；候选/多选含义见说明 | 返回确定性的最佳注释或所有候选。 |

</details>

<span id="api-find_ggchord_loci"></span>
<details><summary><code>find_ggchord_loci()</code></summary>


[完整参考](reference/find_ggchord_loci.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 序列表，含唯一 accver 及正整数 length。可选 genome、genome_id 保留在返回窗口中。 |
| `feature_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 注释表，含 feature_id、accver、start、end，可选 strand；缺失链表示未知。 |
| `targets` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 非空、无重复的 feature_id 字符向量。精确选择，不猜测注释文字或 accession version。 |
| `upstream` | `0` | 数值/整数；范围见说明 | 单个非负整数，表示沿生物学方向在每个目标前后扩展的碱基数。 |
| `downstream` | `0` | 数值/整数；范围见说明 | 单个非负整数，表示沿生物学方向在每个目标前后扩展的碱基数。 |

</details>

<span id="api-find_primer_bindings"></span>
<details><summary><code>find_primer_bindings()</code></summary>


[完整参考](reference/find_primer_bindings.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `sequence` | `必需/互斥入口，见说明` | 序列对象/字符串/完整序列表，见说明 | DNA 文本、命名字符向量，或含 accver、sequence 的数据框。 |
| `primers` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选自定义引物，含 sequence，可选 name、primer_id（每个不同序列一个 ID），或命名字符向量。 |
| `set` | `c("universal", "all", "reference")` | 向量；候选/多选含义见说明 | 随包目录子集。reference 为识别出的随包参考质粒返回权威引物注释。 |
| `match` | `c("exact", "three_prime")` | 向量；候选/多选含义见说明 | 完整精确匹配，或最长精确 3′ 端退火。 |
| `min_annealed_bases` | `15L` | 数值/整数；范围见说明 | 至少 8 的有限整数，three_prime 匹配的最小后缀长度；拒绝小数。 |
| `circular` | `TRUE` | 逻辑值 | 允许跨序列原点匹配。 |
| `unique_only` | `FALSE` | 逻辑值 | 只保留每个目标上仅一个结合位点的引物序列。 |

</details>

<span id="api-find_restriction_sites"></span>
<details><summary><code>find_restriction_sites()</code></summary>


[完整参考](reference/find_restriction_sites.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `sequence` | `必需/互斥入口，见说明` | 序列对象/字符串/完整序列表，见说明 | DNA 文本、命名字符向量，或含 accver、sequence 的数据框。 |
| `enzymes` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选酶名称；数据框是 patterns 的兼容写法。 |
| `patterns` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选自定义定义。 |
| `circular` | `TRUE` | 逻辑值 | 允许识别和切割坐标跨原点回绕。 |
| `database` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选已解析的模式数据库。 |
| `methylation` | `c("auto", "none", "dam_dcm")` | 向量；候选/多选含义见说明 | 甲基化背景。"auto" 为识别出的参考质粒使用已审计背景，其他情况不声称位点被阻断；"none" 关闭阻断；"dam_dcm" 按源支持的 Dam/Dcm 敏感性检查输入序列。 |

</details>

<span id="api-focus_ggchord_data"></span>
<details><summary><code>focus_ggchord_data()</code></summary>


[完整参考](reference/focus_ggchord_data.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 含 accver、length 的序列表。 |
| `ribbon_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选 ggchord 比对表。 |
| `gene_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选基因或注释表。 |
| `loci` | `必需/互斥入口，见说明` | 显式窗口/锚点表或已记录的等价输入 | 数据框，含 accver、start、end，可选 locus_id。 |
| `expand` | `0` | 数值/整数；范围见说明 | 两侧扩展的非负碱基数，可为单值或每个区域一个值。 |
| `boundary` | `c("trim", "drop")` | 向量；候选/多选含义见说明 | 边界策略："trim" 裁剪部分重叠的基因/比对，"drop" 只保留完全位于区域内的行。 |
| `feature_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 统一注释输入，与 gene_data 互斥。 |

</details>

<span id="api-focus_ggchord_track"></span>
<details><summary><code>focus_ggchord_track()</code></summary>


[完整参考](reference/focus_ggchord_track.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 定量输入表，保留元数据和映射值。 |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 序列参照表。 |
| `loci` | `必需/互斥入口，见说明` | 显式窗口/锚点表或已记录的等价输入 | 显式 accver、start、end，可选 locus_id，与 focus_ggchord_data() 相同。 |
| `type` | `c("point", "interval")` | 向量；候选/多选含义见说明 | "point" 要求 pos；"interval" 要求 start/end。 |
| `boundary` | `c("trim", "drop")` | 向量；候选/多选含义见说明 | 部分保留区间使用 "trim"（默认）或 "drop"。 |

</details>

<span id="api-fusion_panel"></span>
<details><summary><code>fusion_panel()</code></summary>


[完整参考](reference/fusion_panel.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `plot` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 横向或纵向Collinear ggchord图形。此处不支持Sankey、Radial、Dotplot、分面及外部ggplot。 |
| `id` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 唯一的非空面板名称；"center"为保留名称。 |
| `side` | `c("right", "top", "left", "bottom")` | 向量；候选/多选含义见说明 | 可选right、top、left或bottom；每个方向最多一个面板。 |
| `size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的相对宽度与高度。指定区域重叠或超出画布时报错，不擅自移动用户明确配置的面板。 |

</details>

<span id="api-fusion_ribbon"></span>
<details><summary><code>fusion_ribbon()</code></summary>


[完整参考](reference/fusion_ribbon.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 保留原始端点的关系数据框。 |
| `mapping` | `NULL` | aes 映射 | 标准ggplot2美学映射，可映射fill、colour、alpha、linewidth和linetype。存在pident且未明确指定fill时，默认使用Identity映射及0/50/80/90/95/100刻度。组合层使用标准scale_fill_*等标度；各面板的角色标度只在其内部生效。 |
| `qpanel` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的固定面板ID，不得与同名数据列同时指定。 |
| `spanel` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的固定面板ID，不得与同名数据列同时指定。 |
| `gradient` | `FALSE` | 逻辑值 | 显式启用query到subject的序列颜色插值，不得与固定fill或fill映射同时使用。 |
| `tip` | `c("none", "subject")` | 向量；候选/多选含义见说明 | "none"或装饰性的"subject"尖头。即使绘制尖头，仍导出完整来源前缘。 |
| `allow_repeated` | `FALSE` | 逻辑值 | 显式允许重复展示source_hsp_id。 |
| `curvature` | `0.3` | 数值/整数；范围见说明 | 控制距离相对于端点间距的比例，范围为[0,1]。 |
| `fill` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 固定样式优先于对应映射。没有pident时，默认填充为steelblue。 |
| `colour` | `NA` | R 对象；约束见说明 | 固定样式优先于对应映射。没有pident时，默认填充为steelblue。 |
| `alpha` | `0.35` | 数值/整数；范围见说明 | 固定样式优先于对应映射。没有pident时，默认填充为steelblue。 |
| `linewidth` | `0.25` | 数值/整数；范围见说明 | 固定样式优先于对应映射。没有pident时，默认填充为steelblue。 |
| `linetype` | `1` | 数值/整数；范围见说明 | 固定样式优先于对应映射。没有pident时，默认填充为steelblue。 |
| `show.legend` | `NA` | R 对象；约束见说明 | 标准ggplot2图例选择：NA按映射自动决定；TRUE/FALSE或命名逻辑向量控制参与图例的美学。 |

</details>

<span id="api-get_chord_layout"></span>
<details><summary><code>get_chord_layout()</code></summary>


[完整参考](reference/get_chord_layout.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `plot` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 一个 ggchord 图。 |
| `build` | `TRUE` | 逻辑值 | 逻辑值，未缓存时是否构建布局，默认 TRUE。 |

</details>

<span id="api-ggchord"></span>
<details><summary><code>ggchord()</code></summary>


[完整参考](reference/ggchord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 含 accver 和 length 的 data.frame/tibble，或已知完整长度的 GRanges、Seqinfo、GInteractions，或已命名的核酸或蛋白 Biostrings XStringSet。未命名序列先通过 as_ggchord_seq() 显式提供 ids。可选 S4 输入经 as_ggchord_seq() 适配。sequence_unit 为 "bp"（省略时的默认值）或 "aa"，同一图使用一个单位。蛋白长度来自完整序列。aa 图接受连接和通用特征；基因、引物及酶切图层需要 bp。任意字符宽度不等于生物序列长度。 |
| `ribbon_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选比对表，包含 qaccver、saccver、qstart、qend、sstart 和 send。只有映射、筛选或统计需要时，才要求 length 和 pident。GInteractions 通过 as_ggchord_link() 转换，不编造相似度或比对长度。 |
| `gene_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 兼容的旧注释输入，第三位置参数含义不变。新代码使用 feature_data；非 NULL 旧输入在构造时提示一次，渲染时不重复。无类别旧表保留有方向基因语义，未安排删除版本。两种注释输入不能同时提供。 |
| `debug` | `FALSE` | 逻辑值 | 逻辑值，是否输出调试信息，默认 FALSE。 |
| `validate` | `c("warn", "error", "none")` | 向量；候选/多选含义见说明 | 字符值，默认 "warn"。控制结构化输入校验，见 validate_ggchord_data："warn" 在发现问题时发出一条汇总警告，并在图对象的 p$ggchord$validation 缓存完整报告；"error" 在严重问题时停止；"none" 跳过诊断性校验，但仍执行防止崩溃的基础结构检查。 |
| `feature_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选统一注释表或 GRanges，要求 accver、start、end，缺少链时为未知（"."）。通用特征层使用全表；默认基因层按 feature_class、type、role 的优先顺序选择 gene/CDS，并要求已知链。tRNA/rRNA 保留为通用特征。图层显式 data 优先；组合 gene 与 feature 层时使用互斥子集，各图层独立保留自己的数据。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 保留用于为已删除的构造函数参数提供明确迁移错误。外观设置使用 labs()、coord_chord() 和 theme()。 |

</details>

<span id="api-ggchord_fusion"></span>
<details><summary><code>ggchord_fusion()</code></summary>


[完整参考](reference/ggchord_fusion.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `center` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 使用Circular或Chord坐标的ggchord图形。 |
| `center_size` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的宽度与高度，按扣除跨面板图例后的组合区域比例指定。NULL自动分配可用空间。 |
| `layout` | `c("fixed", "adaptive")` | 向量；候选/多选含义见说明 | "fixed"使用固定面板区域；"adaptive"还根据真实物理宽高比调整中央区域。两者均按实际面板尺寸重新计算标签和图例，不缩放用户明确设置的字号。 |

</details>

<span id="api-link_ggchord_clusters"></span>
<details><summary><code>link_ggchord_clusters()</code></summary>


[完整参考](reference/link_ggchord_clusters.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `feature_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 注释数据框，含 feature_id、accver、start、end，可选 strand，应符合 as_ggchord_feature()。 |
| `membership` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 含 feature_id 和 cluster_id 的数据框；每个注释只能属于一个组，其余列保留在报告中。 |
| `pairs` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 以 qaccver 和 saccver 明确指定比较的数据框。 |
| `paralogs` | `c("drop", "all")` | 向量；候选/多选含义见说明 | 同一组、同一序列有多个成员时，默认 "drop"，或显式 "all"。 |
| `segment_policy` | `c("error", "span")` | 向量；候选/多选含义见说明 | 默认 "error" 拒绝重复注释 ID；"span" 显式使用已校验不连续注释的外边界。 |
| `max_links` | `10000` | 数值/整数；范围见说明 | 生成行数的正整数上限，默认 10000；超限会在分配笛卡尔积前停止。 |

</details>

<span id="api-map_alignment_position"></span>
<details><summary><code>map_alignment_position()</code></summary>


[完整参考](reference/map_alignment_position.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | as_ggchord_alignment() 返回的对齐长表。 |
| `accver` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 单个对齐行 ID，即使其 source_accver 使用不同名称。 |
| `position` | `必需/互斥入口，见说明` | Position 对象或已记录的字符串 | 数值位置，必须为正整数或 NA。 |
| `inverse` | `FALSE` | 逻辑值 | FALSE 将对齐列映射到源位置；TRUE 将源位置映射到对齐列。 |

</details>

<span id="api-map_ggchord_position"></span>
<details><summary><code>map_ggchord_position()</code></summary>


[完整参考](reference/map_ggchord_position.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `plot` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 使用coord_collinear()、coord_radial_synteny()或coord_dotplot()的ggchord图。 |
| `accver` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 序列标识，单值或每个位置一个。 |
| `position` | `必需/互斥入口，见说明` | Position 对象或已记录的字符串 | 数值基因组位置，inverse = TRUE时为布局x。Radial反查须改为提供包含x及y列的矩阵或数据框，位置必须在指定骨架上，骨架外位置会被拒绝。 |
| `inverse` | `FALSE` | 逻辑值 | 是否将布局 x 反向转换为原生序列位置。 |
| `source` | `FALSE` | 逻辑值 | 聚焦后是否使用原始源坐标。 |
| `axis` | `NULL` | 可空 R 对象；具体类型/继承见说明 | Dotplot映射轴为"x"（query）或"y"（subject）。Dotplot必须指定，线性与径向布局省略。Dotplot间隙的反查结果为NA，返回值记录所选轴。 |

</details>

<span id="api-merge_ggchord_ribbons"></span>
<details><summary><code>merge_ggchord_ribbons()</code></summary>


[完整参考](reference/merge_ggchord_ribbons.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `ribbon_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 连接格式的数据框。 |
| `max_gap` | `0` | 数值/整数；范围见说明 | 数值，默认 0；两条序列上允许合并的最大块间距，单位 bp。 |
| `min_pident_difference` | `0` | 数值/整数；范围见说明 | 数值，默认 0；大于 0 时，两块 pident 差值须不超过该值。 |
| `require_same_orientation` | `TRUE` | 逻辑值 | 逻辑值，默认 TRUE；仅合并相对方向相同且边界推进相容的块。两端同时反转保持相对方向；FALSE 明确允许更粗略的混向汇总。 |
| `group_by` | `c("qaccver", "saccver")` | 向量；候选/多选含义见说明 | 额外分组列，默认 c("qaccver", "saccver")；query/subject ID 始终是分组键。 |
| `extra_columns` | `c("consistent", "first")` | 向量；候选/多选含义见说明 | 多行合并时额外列的处理：默认 "consistent" 仅在所有源行一致时保留，否则记同类型缺失值；"first" 保留首行值。 |

</details>

<span id="api-optimize_ggchord_layout"></span>
<details><summary><code>optimize_ggchord_layout()</code></summary>


[完整参考](reference/optimize_ggchord_layout.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 含 accver、length 的序列数据。 |
| `ribbon_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | ggchord 连接格式的比对数据。 |
| `seq_order` | `NULL` | 几何数值/向量/已记录格式 | 可选完整初始序列顺序，默认使用 seq_data 行顺序。 |
| `seq_orientation` | `NULL` | 几何数值/向量/已记录格式 | 可选初始方向 1 或 -1，支持单值、完整向量或命名向量。 |
| `optimize` | `c("order", "orientation")` | 向量；候选/多选含义见说明 | 选择 "order" 和/或 "orientation"。 |
| `weight` | `c("length", "count", "pident")` | 向量；候选/多选含义见说明 | 目标权重："length"、"count" 或 "pident"。 |
| `max_iter` | `50L` | 数值/整数；范围见说明 | 改善搜索轮数的正整数上限。 |

</details>

<span id="api-prepare_ggchord_data"></span>
<details><summary><code>prepare_ggchord_data()</code></summary>


[完整参考](reference/prepare_ggchord_data.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | as_ggchord_seq() 接受的序列输入。 |
| `ribbon_data` | `NULL` | 表格/受支持对象；data 还可为函数 | as_ggchord_link() 接受的可选成对区间。 |
| `gene_data` | `NULL` | 表格/受支持对象；data 还可为函数 | as_ggchord_feature() 接受的可选有方向基因输入；无方向区间放入 features。 |
| `features` | `list()` | 命名列表 | 注释或成对区间输入的命名列表。 |
| `links` | `list()` | 命名列表 | 注释或成对区间输入的命名列表。 |
| `mapping` | `list()` | aes 映射 | 命名列表，可含 seq、ribbon、gene 列映射，以及与集合同名的 features/links 映射列表；每项列映射为命名字符向量。 |
| `steps` | `list()` | 命名列表 | 明确辅助步骤参数的命名列表，支持 clean、filter、deduplicate、merge、bundle、optimize。数据由准备函数提供，不可覆盖。 |
| `feature_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 统一注释表，与 gene_data 互斥。 |

</details>

<span id="api-primer_catalog"></span>
<details><summary><code>primer_catalog()</code></summary>


[完整参考](reference/primer_catalog.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `set` | `c("all", "universal")` | 向量；候选/多选含义见说明 | 返回全部记录或常用 universal 子集。 |
| `aliases` | `FALSE` | 逻辑值 | 是否返回源别名记录，而非序列目录。 |

</details>

<span id="api-project_ggchord_features"></span>
<details><summary><code>project_ggchord_features()</code></summary>


[完整参考](reference/project_ggchord_features.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `feature_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 结构域/注释数据框，含 feature_id、parent_id（或每项单值的 parent_ids 列表列）、start、end。 |
| `cds` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 由 build_ggchord_cds() 创建的映射。 |
| `from` | `c("protein", "genome")` | 向量；候选/多选含义见说明 | 默认 "protein" 或 "genome"。 |
| `codon_boundary` | `c("error", "expand")` | 向量；候选/多选含义见说明 | 反向投影策略："error" 或 "expand"。 |

</details>

<span id="api-project_ggchord_links"></span>
<details><summary><code>project_ggchord_links()</code></summary>


[完整参考](reference/project_ggchord_links.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `link_data` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 数据框，以 qfeature_id、sfeature_id 标识 CDS 模型，并提供有序闭区间蛋白坐标 qstart、qend、sstart、send。两端残基跨度须相等；若存在 gapopen 则必须为零，源 length 须等于残基跨度。 |
| `cds` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 由 build_ggchord_cds() 创建的映射。 |
| `ungapped` | `FALSE` | 逻辑值 | 须显式设为 TRUE，确认无缺口残基对应，不能仅由端点跨度相等推断。 |

</details>

<span id="api-read_aligned_fasta"></span>
<details><summary><code>read_aligned_fasta()</code></summary>


[完整参考](reference/read_aligned_fasta.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `file` | `必需/互斥入口，见说明` | 字符路径、标识或向量；见约束 | 单个本地已对齐 FASTA 文件路径，不拼接多个独立对齐。 |
| `alphabet` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 显式指定 "DNA"、"RNA" 或 "AA"。接受 IUPAC 核酸含糊码及蛋白含糊/稀有残基码；* 仅用于 AA。不猜测字母表。 |
| `source_data` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选数据框，每个对齐 accver 对应一行，包含 start（首个非 gap 源位置）和 strand（"+" 或 "-"）。可选 source_accver 声明不同的源 ID，可选字符列 sequence 提供完整正向源字符串以供核对。坐标必须为正且精确的整数。输入元数据保留在 ggchord_conversion 属性中。 |
| `source_file` | `FALSE` | 逻辑值 | 逻辑值，是否在每个字符行额外保留 .source_file。 |

</details>

<span id="api-read_bed"></span>
<details><summary><code>read_bed()</code></summary>


[完整参考](reference/read_bed.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `file` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径向量或通配模式，在 files 之前读取。 |
| `files` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径或通配模式，用于多个 BED 文件。 |
| `source_file` | `FALSE` | 逻辑值 | 逻辑值，是否添加 .source_file 列，默认 FALSE。 |

</details>

<span id="api-read_bedpe"></span>
<details><summary><code>read_bedpe()</code></summary>


[完整参考](reference/read_bedpe.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `file` | `NULL` | 字符路径、标识或向量；见约束 | 可选的本地路径向量或通配模式，在files之前读取。 |
| `files` | `NULL` | 字符路径、标识或向量；见约束 | 可选的其他本地路径或通配模式。 |
| `source_file` | `FALSE` | 逻辑值 | 逻辑值，是否添加.source_file列，默认FALSE。 |
| `unmapped` | `c("error", "keep")` | 向量；候选/多选含义见说明 | 不完整锚点的处理方式：默认"error"，或显式"keep"。 |

</details>

<span id="api-read_blast"></span>
<details><summary><code>read_blast()</code></summary>


[完整参考](reference/read_blast.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `file` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径向量或通配模式，在 files 之前读取。 |
| `files` | `NULL` | 字符路径、标识或向量；见约束 | 可选 BLAST 表格文件字符向量，支持路径和通配模式；读取所有匹配文件并合并为一个数据框。 |
| `format` | `c("auto", "outfmt6", "outfmt7", "custom")` | 向量；候选/多选含义见说明 | 格式选项："auto"（默认）自动检测列布局；"outfmt6" 要求标准 12 列；"outfmt7" 解析 # Fields: 声明；"custom" 要求提供 col_names。 |
| `col_names` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的文件列名字符向量，用于 format = "custom" 或覆盖自动检测。 |
| `comment` | `"#"` | 字符/枚举；取值见说明 | 注释字符，默认 "#"；BLAST outfmt 7 的头部注释行以 # 开始。 |
| `source_file` | `FALSE` | 逻辑值 | 逻辑值，读取一个或多个文件时是否添加 .source_file 列，默认 FALSE。 |
| `...` | `必需/互斥入口，见说明` | 上游转发参数 | 传给 utils::read.table() 的其他参数，例如 na.strings。 |

</details>

<span id="api-read_fasta_lengths"></span>
<details><summary><code>read_fasta_lengths()</code></summary>


[完整参考](reference/read_fasta_lengths.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `file` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径向量或通配模式，在 files 之前读取。 |
| `files` | `NULL` | 字符路径、标识或向量；见约束 | 可选 FASTA 文件字符向量，支持路径和通配模式；读取并合并所有匹配文件。 |
| `header_delim` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选分隔字符串；提供后在每个出现位置拆分头部，只保留第一段。 |
| `source_file` | `FALSE` | 逻辑值 | 逻辑值，读取一个或多个文件时是否添加 .source_file 列，默认 FALSE。 |
| `sequence_unit` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选单位声明 "bp" 或 "aa"，提供后在结果中添加该列。函数不推断字母表。省略时保留既有输出列，ggchord 按 bp 解释。 |
| `keep_sequence` | `FALSE` | 逻辑值 | 逻辑值，默认 FALSE。在 sequence 中保留真实 FASTA 字符，仅去除空白，保留大小写。不推断字母表、反向互补或比对列。若有 gap 字符，也计入 length；将有 gap 的输入用于基因组坐标前，必须提供显式的比对列到源位置映射。 |

</details>

<span id="api-read_genbank"></span>
<details><summary><code>read_genbank()</code></summary>


[完整参考](reference/read_genbank.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `file` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径向量或通配模式，在 files 之前读取。 |
| `files` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径或通配模式，用于多个文件。 |
| `feature_types` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的注释类型，NULL 保留全部。 |
| `lengths` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选的已知完整长度命名数值向量，遵循 as_ggchord_seq()。 |
| `python` | `"python3"` | 字符路径、标识或向量；见约束 | 解释器名称或路径，默认 "python3"，必须事先安装可选 Biopython；本函数不修改环境。 |
| `unsupported` | `c("error", "report")` | 向量；候选/多选含义见说明 | 对无法安全表示为本地区间的位置，选择 "error" 或 "report"。 |
| `source_file` | `FALSE` | 逻辑值 | 是否在输出表添加 .source_file 列，默认 FALSE；无论此选项如何，报告均保留文件来源。 |

</details>

<span id="api-read_gff3"></span>
<details><summary><code>read_gff3()</code></summary>


[完整参考](reference/read_gff3.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `file` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径向量或通配模式，在 files 之前读取。 |
| `files` | `NULL` | 字符路径、标识或向量；见约束 | 可选 GFF3 文件字符向量，支持路径和通配模式；读取并合并所有匹配文件。 |
| `feature_types` | `"CDS"` | 字符/枚举；取值见说明 | 保留的注释类型字符向量，默认 "CDS"。 |
| `anno_from` | `c("product", "Name", "gene", "ID")` | 向量；候选/多选含义见说明 | GFF3 属性键字符向量，按顺序尝试填充 anno 列。 |
| `unstranded` | `c("plus", "drop", "keep")` | 向量；候选/多选含义见说明 | "plus"（默认值）将未知链设为 +，"drop" 删除未知链。使用 "keep" 并通过 geom_feature() 可保留为无方向的 . 区间；source_strand 始终记录原始链。 |
| `source_file` | `FALSE` | 逻辑值 | 逻辑值，读取一个或多个文件时是否添加 .source_file 列，默认 FALSE。 |

</details>

<span id="api-read_paf"></span>
<details><summary><code>read_paf()</code></summary>


[完整参考](reference/read_paf.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `file` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径向量或通配模式，在 files 之前读取。 |
| `files` | `NULL` | 字符路径、标识或向量；见约束 | 可选路径或通配模式，用于多个 PAF 文件。 |
| `source_file` | `FALSE` | 逻辑值 | 逻辑值，是否添加 .source_file 列，默认 FALSE。 |

</details>

<span id="api-segment_ggchord_data"></span>
<details><summary><code>segment_ggchord_data()</code></summary>


[完整参考](reference/segment_ggchord_data.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 完整序列表，含唯一 accver 和整数长度。 |
| `ribbon_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选成对区间表。 |
| `gene_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选基因或通用注释表，支持分段。 |
| `windows` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 窗口表，含 accver、start、end、唯一 locus_id，可选 row_id（默认源 accver）与 orientation（1 或 -1，默认 1）。边界须有序、整数且位于序列内。 |
| `point_data` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 可选点连接，含 qaccver、saccver、qpos、spos。 |
| `cross_gap` | `c("split", "drop")` | 向量；候选/多选含义见说明 | 被选定窗口截断的区间连接处理方式。 |
| `max_links` | `10000` | 数值/整数；范围见说明 | 区间/点连接候选副本数的正整数上限；有意使用更大重叠视图时显式提高。 |
| `feature_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 统一注释输入，与 gene_data 互斥。 |

</details>

<span id="api-validate_ggchord_data"></span>
<details><summary><code>validate_ggchord_data()</code></summary>


[完整参考](reference/validate_ggchord_data.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 必需 data.frame/tibble，包含基本序列列 accver、length。 |
| `ribbon_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选 data.frame/tibble，比对列为 qaccver、saccver、qstart、qend、sstart、send。 |
| `gene_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 可选 data.frame/tibble，基因列为 accver、start、end、strand；anno 可选。 |
| `strict` | `FALSE` | 逻辑值 | 逻辑值。TRUE 时收集问题后遇严重错误停止；默认 FALSE 返回完整报告而不停止。 |
| `check_coordinates` | `TRUE` | 逻辑值 | 逻辑值，默认 TRUE，是否检查连接/基因坐标位于 [1, sequence length] 内。 |
| `check_duplicates` | `TRUE` | 逻辑值 | 逻辑值，默认 TRUE，是否检查完全重复、近似重复及高度重叠记录。互惠重叠要求交集覆盖两区间各至少 90%，包括端点递减的连接。反向注释边界可能跨圆形原点，不参与线性重叠诊断，但仍报告完全重复。 |
| `check_self_links` | `TRUE` | 逻辑值 | 逻辑值，默认 TRUE，是否标记满足 qaccver == saccver 的自身连接。 |
| `feature_data` | `NULL` | 表格/受支持对象；data 还可为函数 | 统一注释表，与 gene_data 互斥，要求 accver、start、end。缺少链时为未知，接受 +、-、.、+/-。报告使用表名 "feature" 和 data_summary$n_feature，不修改输入。 |

</details>

<span id="api-validate_ggchord_hierarchy"></span>
<details><summary><code>validate_ggchord_hierarchy()</code></summary>


[完整参考](reference/validate_ggchord_hierarchy.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `seq_data` | `必需/互斥入口，见说明` | 表格/受支持对象；data 还可为函数 | 包含 accver 和 genome 外键的数据框。 |
| `genomes` | `必需/互斥入口，见说明` | R 对象；约束见说明 | 包含唯一 genome 主键的数据框。 |
| `strict` | `FALSE` | 逻辑值 | 逻辑值，关系错误时停止；genome 没有 contig 的警告不会导致失败。 |

</details>

<span id="api-view_ggchord"></span>
<details><summary><code>view_ggchord()</code></summary>


[完整参考](reference/view_ggchord.html)


| 参数 | 签名默认 | 类型/取值规则 | 含义、范围与限制 |
| --- | --- | --- | --- |
| `plot` | `ggplot2::last_plot()` | R 对象；约束见说明 | ggchord 或 ggplot 对象，默认 last_plot()。 |
| `width` | `12.39` | 数值/整数；范围见说明 | 正输出宽度。未显式提供尺寸时，coord_chord() 使用 12.39 英寸，coord_circular() 从注释外包区域推导宽度。 |
| `height` | `9.71` | 数值/整数；范围见说明 | 正输出高度。未提供尺寸时，coord_chord() 使用 9.71 英寸，coord_circular() 推导宽高。NULL 在保持坐标单位等比例时拟合序列、标签和图例。显式宽高始终保留。 |
| `units` | `c("in", "cm", "mm", "px")` | 向量；候选/多选含义见说明 | 输出单位："in"、"cm"、"mm" 或 "px"。 |
| `device` | `c("png", "svg")` | 向量；候选/多选含义见说明 | 预览设备 "png" 或 "svg"；SVG 需要可选 svglite 或可用的基础 Cairo SVG 设备。 |
| `dpi` | `144` | 数值/整数；范围见说明 | 正栅格分辨率，默认 144。 |
| `bg` | `NULL` | 可空 R 对象；具体类型/继承见说明 | 传给 ggsave() 的可选背景色。 |
| `viewer` | `c("auto", "ide", "browser", "none")` | 向量；候选/多选含义见说明 | 打开方式："auto" 优先 getOption("viewer")，交互会话中可回退浏览器；"ide"、"browser"、"none" 选择明确目标。"none" 适合脚本和测试。 |

</details>


## 上游重导出与对象 {#upstream}

这些入口保留上游 ggplot2/grid 语义；它们不是新增的 ggchord 图层。after_stat/after_scale 的阶段不同，不能由 identical 结果推断可互换。

<span id="api-aes"></span>
<details><summary><code>aes</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-after_scale"></span>
<details><summary><code>after_scale</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-after_stat"></span>
<details><summary><code>after_stat</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-annotate"></span>
<details><summary><code>annotate</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `geom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `x` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `y` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `xmin` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `xmax` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `ymin` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `ymax` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `xend` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `yend` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `na.rm` | `FALSE` | 逻辑值; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-arrow"></span>
<details><summary><code>arrow</code> · grid</summary>

上游入口参考: [Rd](reference/grid-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `angle` | `30` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `length` | `unit(0.25, "inches")` | grid 物理单位对象; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `ends` | `"last"` | 字符/枚举；取值见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `type` | `"open"` | 字符/枚举；取值见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-element_blank"></span>
<details><summary><code>element_blank</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |

</details>

<span id="api-element_line"></span>
<details><summary><code>element_line</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `colour` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `linewidth` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `linetype` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `lineend` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `color` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `linejoin` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `arrow` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `arrow.fill` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `inherit.blank` | `FALSE` | 逻辑值; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `size` | `deprecated()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-element_rect"></span>
<details><summary><code>element_rect</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `fill` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `colour` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `linewidth` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `linetype` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `color` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `linejoin` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `inherit.blank` | `FALSE` | 逻辑值; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `size` | `deprecated()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-element_text"></span>
<details><summary><code>element_text</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `family` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `face` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `colour` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `size` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `hjust` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `vjust` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `angle` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `lineheight` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `color` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `margin` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `debug` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `inherit.blank` | `FALSE` | 逻辑值; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-expansion"></span>
<details><summary><code>expansion</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `mult` | `0` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `add` | `0` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-ggsave"></span>
<details><summary><code>ggsave</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `filename` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot` | `get_last_plot()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `device` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `path` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `scale` | `1` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `width` | `NA` | 几何数值/向量/已记录格式; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `height` | `NA` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `units` | `c("in", "cm", "mm", "px")` | 向量；候选/多选含义见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `dpi` | `300` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `limitsize` | `TRUE` | 逻辑值; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `bg` | `NULL` | 可空 R 对象；具体类型/继承见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `create.dir` | `FALSE` | 逻辑值; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-ggtitle"></span>
<details><summary><code>ggtitle</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `label` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `subtitle` | `waiver()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-guides"></span>
<details><summary><code>guides</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-labs"></span>
<details><summary><code>labs</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `title` | `waiver()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `subtitle` | `waiver()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `caption` | `waiver()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `tag` | `waiver()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `dictionary` | `waiver()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `alt` | `waiver()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `alt_insight` | `waiver()` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-last_plot"></span>
<details><summary><code>last_plot</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |

</details>

<span id="api-margin"></span>
<details><summary><code>margin</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `t` | `0` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `r` | `0` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `b` | `0` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `l` | `0` | 数值/整数；范围见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `unit` | `"pt"` | 字符/枚举；取值见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-position_identity"></span>
<details><summary><code>position_identity</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |

</details>

<span id="api-PositionIdentity"></span>
<details><summary><code>PositionIdentity</code> · object</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |

</details>

<span id="api-rel"></span>
<details><summary><code>rel</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-theme"></span>
<details><summary><code>theme</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `...` | `必需/互斥入口` | 上游转发参数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `line` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `rect` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `text` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `title` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `point` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `polygon` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `geom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `spacing` | `必需/互斥入口` | 几何数值/向量/已记录格式; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `margins` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `aspect.ratio` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.title` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.title.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.title.x.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.title.x.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.title.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.title.y.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.title.y.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text.x.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text.x.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text.y.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text.y.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text.theta` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.text.r` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.x.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.x.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.y.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.y.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.theta` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.r` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.x.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.x.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.y.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.y.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.theta` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.r` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length.x.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length.x.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length.y.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length.y.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length.theta` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.ticks.length.r` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length.x.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length.x.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length.y.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length.y.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length.theta` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.minor.ticks.length.r` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line.x.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line.x.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line.y.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line.y.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line.theta` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `axis.line.r` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.background` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.margin` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.spacing` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.spacing.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.spacing.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.key` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.key.size` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.key.height` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.key.width` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.key.spacing` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.key.spacing.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.key.spacing.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.key.justification` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.frame` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.ticks` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.ticks.length` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.axis.line` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.text` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.text.position` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.title` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.title.position` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.position` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.position.inside` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.direction` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.byrow` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.justification` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.justification.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.justification.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.justification.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.justification.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.justification.inside` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.location` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.box` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.box.just` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.box.margin` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.box.background` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `legend.box.spacing` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.background` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.border` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.spacing` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.spacing.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.spacing.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.grid` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.grid.major` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.grid.minor` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.grid.major.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.grid.major.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.grid.minor.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.grid.minor.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.ontop` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.widths` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `panel.heights` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.background` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.title` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.title.position` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.subtitle` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.caption` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.caption.position` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.tag` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.tag.position` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.tag.location` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `plot.margin` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.background` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.background.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.background.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.clip` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.placement` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.text` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.text.x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.text.x.bottom` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.text.x.top` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.text.y` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.text.y.left` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.text.y.right` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.switch.pad.grid` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `strip.switch.pad.wrap` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `complete` | `FALSE` | 逻辑值; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `validate` | `TRUE` | 逻辑值; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-unit"></span>
<details><summary><code>unit</code> · grid</summary>

上游入口参考: [Rd](reference/grid-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |
| `x` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `units` | `必需/互斥入口` | R 对象；约束见说明; 遵循上游接口；具体类型、单位和允许值见其帮助。 |
| `data` | `NULL` | 表格/受支持对象；data 还可为函数; 遵循上游接口；具体类型、单位和允许值见其帮助。 |

</details>

<span id="api-waiver"></span>
<details><summary><code>waiver</code> · ggplot2</summary>

上游入口参考: [Rd](reference/ggplot2-helpers.html)

| 参数 | 签名默认 | 类型/约束 |
| --- | --- | --- |

</details>
