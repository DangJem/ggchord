<div class="task-intro"><p>按真实导出尺寸检查图形，再用标准 ggsave() 保存，并了解布局表与诊断可以证明什么。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。

{{prepare:setup,docstyle,collinear}}

## 按目标尺寸保存 {#size}
先在[三种模式](modes.html)中创建 `p_collinear`：

{{code:save}}

RStudio Plots 使用当前面板尺寸，标签与图例可能随面板调整而重新排列。
最终 `ggsave()` 前，按相同显式尺寸预览：

{{code:preview}}

`view_ggchord()` 返回临时图片路径，不保存最终导出文件。
网站、README 和 vignette 共用标准 `ggsave()` 生成的 PNG，尺寸见[图集](gallery.html)，
不截取临时调整后的交互面板。最终矢量输出应显式选择目标 PDF/SVG 设备；字体可用性
和嵌入仍取决于设备。

## 按字段名消费几何 {#geometry}
{{code:provenance}}

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
{{code:layoutresults}}

先在[三种模式](modes.html#collinear)创建 `p_collinear`。
`full_layout` 保留完整的已记录 getter 结果，`layer_tables` 选择序列/基因表与已解析
原始数据，不是原始文件档案。layer ID 和源行角色区分独立同类图层。
空组件可没有固定列原型，应按字段名读取。

`layout_choice` 以计数权重选择顺序/方向，因为本例的最小绘图表有意未保留比对长度统计列。
返回的选择应显式用于序列设置，不修改原图，也不确立生物学方向。
采用前应保存该选择并对照结果。文档精确尺寸见[资产清单](../assets/plots.csv)。

## 完成后检查 {#check-result}
预览与保存使用相同宽高。检查文字、图例与警告；一次测量成功不代表跨设备字体和发布质量已认证。

<div class="next-step">继续：<a href="faq.html">常见问题 →</a></div>
