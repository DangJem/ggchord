<div class="task-intro"><p>按真实导出尺寸检查图形，再用标准 ggsave() 保存，并了解布局表与诊断可以证明什么。</p></div>

## 本页准备 {#prepare}
用四条随包序列画出共线性骨架，后面用同一张图演示预览和保存。

{{code:frontlinearbase}}
{{figure:frontlinearbase}}

## 按目标尺寸保存 {#size}
下面以本页创建的 `p_linear_base` 为例。有自己的图时，把它换成你的绘图对象。

{{code:frontpreview}}

圆图省略宽高可自动估计尺寸。弦图和共线性图可指定宽度并设 `height = NULL`，按图形估计高度。网页缩略图还会按正文栏缩放；仅增大 dpi 不会让网页中的字变大。

确认字号、标签和图例后，使用相同尺寸保存。

{{code:frontsave}}

`view_ggchord()` 返回临时图片路径，不保存最终导出文件。
网站、README 和 vignette 共用由 `view_ggchord()` 生成的预览 PNG，尺寸见[图集](gallery.html)，
不截取临时调整后的交互面板。最终矢量输出应显式选择目标 PDF/SVG 设备；字体可用性
和嵌入仍取决于设备。

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
`optimize_ggchord_layout()` 按文档执行主动选择的调整，应保留副本并查看结果。
预览和诊断使用独立渲染副本。直接改动共享 ggproto 图层环境仍遵循普通 ggplot2
引用语义，不会自动隔离为一次独立编辑。

## 完整布局、表格与主动选择的排序 {#consumption}

前面的结果表已足够追踪骨架图层。需要附上原始输入行时，给 `export_ggchord_layout()` 加 `original_data = TRUE`。需要查看完整布局或改变顺序时，分别查阅 [布局契约](reference/ggchord_layout_contract.html) 与 [排序和方向](reference/optimize_ggchord_layout.html)。

## 完成后检查 {#check-result}
预览与保存使用相同宽高。检查文字、图例与警告；一次测量成功不代表跨设备字体和发布质量已认证。

<div class="next-step">继续：<a href="faq.html">常见问题 →</a></div>

## 测量标签接触 {#label-quality}

`diagnose_ggchord(p, width = 10, height = 8)$label_quality` 报告实测面板文字对象、注释引线、文字框接触及物理引线长度。`summary` 区分图层策略保留/隐藏的源标签、文字对象和字形框；`conflicts`、`text` 和 `leaders` 定位每个接触。弧形字符串内部字形共享文字对象身份，折线共享端点不算交叉。检查表格时也应查看实际 PNG。

酶切标签的实绘、外围排版和引线路由共享粗体酶名、普通坐标及物理间隙的目标设备包络；有效主题或显式文字设置在定位之前进入布局。

这些是旋转文字包围框，不是精确墨迹轮廓，引线与框接触可能需要人工复核；空格及标签背景也不等同已认证的字形形状。`feature_contacts` 另检查已登记 gene/feature 多边形包围范围，包括透明多边形，分别记录有同源证据的自身标签贴合（`own_source`）、其他构造输入记录接触（`other_source`）及证据不足的关联（`unresolved`）；不按内容相等猜来源。`feature_polygons` 给出实测顶点，`feature_parts_unchecked` 列出跳过的非多边形线条/标记。线宽、可见性/遮挡、外层装饰及未登记 grob 仍未检查，多边形接触为零不能认证这些项目。诊断不会靠隐藏标签或改变源数据改善布局。面板已测量，或已检查类别中的零接触，都不能证明全部布局约束通过。

SnapGene 作为清晰注释层级与规整引线的显示效果参考，不要求取得该应用的配对导出。使用 ggchord 内部控制变量对照，评价内容保留、物理尺寸、拥挤和可读性的变化。

启用通常的重叠策略时，若自动定位的圆图平行标签（含 `1 / 2856` 这类接缝字符串）文字框接触 feature 多边形，可在有界范围内向内或向外调整径向间距，同时保护其他文字和骨架；随后注释引线绕开已确定的轴文字。`axis_clearance` 表记录实际状态和以毫米计的有符号位移。刻度坐标和标签内容保持；显式提供 `axis.text.offset`（包括与预设相等的值）以及显式 scale breaks 会禁止此调整。这是局部留白规则，不是任意轴布局的质量保证。

圆图外侧 callout 离开 feature band 时也会保护实测原点轴与内侧特征文字。只对与旋转文字框接触的路线绕行，真实起点、最终标签连接点及文字内容保持不变。有界搜索对全部已登记静态文字复核候选，`leader_text_status` 报告已清晰、已绕行或未解决；未解决路线保留证据供复核。无用的零长度占位段会移除，不当作真实引线或 marker 观测。

圆形小画布中，已经交叉的外侧引线可在实测文字与其他引线约束下有限绕行，真实起点、文字边缘连接点和标签位置保持不变。未解决路线仍在诊断中保留；受检交叉为零不表示任意输入或精确墨迹均获保证。

`scale` 与 `adaptive` 策略也会调整原生序列、点链接、ribbon 和 alignment 图层中省略的固定线宽。显式线宽（包括与默认值相等的值）保持固定，映射线宽保留其标度规则；`fixed` 和未添加策略的图保持既有样式。这不表示所有几何常量或原生主题派生的美学表达式都统一相乘。
