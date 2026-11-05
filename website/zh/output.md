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
