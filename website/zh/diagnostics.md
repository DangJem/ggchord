{{code:courseblockbj}}

{{prepare:frontlinearbase}}

## 按字段名消费几何 {#geometry}
{{code:frontgeometry}}

{{output:frontgeometry}}

`export_ggchord_layout()` 是推荐的分层表格接口；`get_chord_layout()` 保留完整的
已记录布局，两者都未弃用。选择字段名，不依赖固定列位置。
`source_row` 指已解析图层输入行，ID 是构建期资源，不是永久生物学标识。
详见完整[布局结果说明](reference/ggchord_layout_contract.html)。

## 理解诊断证据 {#diagnostics}
`diagnose_ggchord()` 区分面板已测量与几何收敛。面板已测量不表示文字无碰撞、字体
已嵌入或跨应用表现已认证。应检查警告、文字/图例测量及放置原因。
有界物理重建未收敛时，绘制会发出警告，几何度量仍为估计值；诊断会捕获该警告，
并保留 `physical_layout` 检查项。可扩大画布或调整标签/图例设置后重新测量。
`optimize_ggchord_layout()` 按文档执行主动选择的调整，应保留副本并查看结果。
预览和诊断使用独立渲染副本。直接改动共享 ggproto 图层环境仍遵循普通 ggplot2
引用语义，不会自动隔离为一次独立编辑。

## 完整布局、表格与主动选择的排序 {#consumption}

前面的结果表已足够追踪骨架图层。需要附上原始输入行时，给 `export_ggchord_layout()` 加 `original_data = TRUE`。需要查看完整布局或改变顺序时，分别查阅 [布局契约](reference/ggchord_layout_contract.html) 与 [排序和方向](reference/optimize_ggchord_layout.html)。

## 测量标签接触 {#label-quality}

`diagnose_ggchord(p, width = 10, height = 8)$label_quality` 报告实测面板文字对象、注释引线、文字框接触及物理引线长度。`summary` 区分图层策略保留/隐藏的源标签、文字对象和字形框；`conflicts`、`text` 和 `leaders` 定位每个接触。`measurement` 条目分别记录面板测量、几何收敛和物理重建次数。圆图可在两次测量后试一次有界的尺度预测，`refinement_acceleration` 记录结果；只有重新构建并通过原实测收敛门槛才接受。拒绝后保留通常的重建序列，不能计为通过。普通重建在各视图使用同一有界预算。弧形字符串内部字形共享文字对象身份，折线共享端点不算交叉。检查表格时也应查看实际 PNG。

已登记文字在目标设备开启时使用该设备的物理度量，包含普通标签及圆弧字形步距；无活动设备的布局仍保留测量后备方式。字体替换和嵌入仍须另行检查。酶切标签的实绘、外围排版和引线路由共享粗体酶名、普通坐标及物理间隙的目标设备包络；有效主题或显式文字设置在定位之前进入布局。

这些是旋转文字包围框，不是精确墨迹轮廓，引线与框接触可能需要人工复核；空格及标签背景也不等同已认证的字形形状。`feature_contacts` 另检查已登记 gene/feature 多边形包围范围，包括透明多边形，分别记录有同一输入来源证据的自身标签贴合（`own_source`）、其他构造输入记录接触（`other_source`）及证据不足的关联（`unresolved`）；不按内容相等猜来源。`feature_polygons` 给出实测顶点，`feature_parts_unchecked` 列出跳过的非多边形线条/标记。线宽、可见性/遮挡、外层装饰及未登记 grob 仍未检查，多边形接触为零不能认证这些项目。诊断不会靠隐藏标签或改变源数据改善布局。面板已测量，或已检查类别中的零接触，都不能证明全部布局约束通过。

SnapGene 作为清晰注释层级与规整引线的显示效果参考，不要求取得该应用的配对导出。使用 ggchord 内部控制变量对照，评价内容保留、物理尺寸、拥挤和可读性的变化。

自动序列轴的轴线和长短刻度绘制在第一个 gene、feature 或 primer 图层之前；不透明特征填充会覆盖相交刻度，透明特征仍保持透明。坐标文字按实测留白放置并绘制在前景。用户各图层仍按原有添加顺序排列，包括显式添加的 `geom_seq_scale()` 比例尺；刻度坐标、长度和来源特征均不变。

启用通常的重叠策略时，若自动定位的圆图平行标签（含 `1 / 2856` 这类接缝字符串）文字框接触 feature 多边形或固定的注释引线根，可在有界范围内向内或向外调整径向间距，同时保护其他文字和骨架；随后内部及外部注释引线绕开已确定的轴文字。会随路由改变的中间折点不反过来影响轴定位。自动端点字符串还可在同一刻度的两个可读侧选择位置，以不超过 4 mm 的切向间距避开固定的 feature 出口射线。`axis_clearance` 表分别记录状态、径向有符号位移、`label_side` 和 `tangent_gap_mm`。刻度坐标和标签内容保持；显式提供 `axis.text.offset`（包括与预设相等的值）以及显式 scale breaks 会禁止此调整。这是局部留白规则，不是任意轴布局的质量保证。 圆图圈外特征和限制位点引线随后可联合尝试直线或单拐点出口，源起点和标签接入点不变；只有整组避开文字且不新增引线交叉时才接受。失败组保持原路线及覆盖显示策略。 显式 straight 和 elbow 限制位点样式保留其已选路线。

圆图外侧 callout 离开 feature band 时也会保护实测原点轴与内侧特征文字。只对与旋转文字框接触的路线绕行，真实起点、最终标签连接点及文字内容保持不变。有界搜索对全部已登记静态文字复核候选，`leader_text_status` 报告已清晰、已绕行或未解决；未解决路线保留证据供复核。无用的零长度占位段会移除，不当作真实引线或 marker 观测。

圆形小画布中，已经交叉的外侧引线可在实测文字与其他引线约束下有限绕行，真实起点、文字边缘连接点和标签位置保持不变。未解决路线仍在诊断中保留；受检交叉为零不表示任意输入或精确墨迹均获保证。

`scale` 与 `adaptive` 策略也会调整原生序列、点链接、ribbon 和 alignment 图层中省略的固定线宽。显式线宽（包括与默认值相等的值）保持固定，映射线宽保留其标度规则；`fixed` 和未添加策略的图保持既有样式。这不表示所有几何常量或原生主题派生的美学表达式都统一相乘。

{{prepare:comparescaledata}}

## 目标尺寸诊断与明确顺序优化 {#target-size}

diagnose_ggchord 检查一个指定 PNG 尺寸，返回测量、文字与 guide 提示，不修改原图。export_ggchord_layout 按层导出数据；original_data=TRUE 保留源行。get_chord_layout 返回完整公共布局。optimize_ggchord_layout 返回建议顺序、方向和前后分数，不改输入也不自动把建议应用到图上。优化的目标是连接可读性，不是生物学上更正确的基因组顺序。

{{code:comparediagnostics}}
{{console:comparediagnostics}}
