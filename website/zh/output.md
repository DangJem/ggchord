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
详见完整[布局契约](reference/ggchord_layout_contract.html)。

## 理解诊断证据 {#diagnostics}
`diagnose_ggchord()` 区分面板已测量与几何收敛。面板已测量不表示文字无碰撞、字体
已嵌入或跨应用表现已认证。应检查警告、文字/图例测量及放置原因。
`optimize_ggchord_layout()` 按文档执行主动选择的调整，应保留副本并查看结果。
预览和诊断使用独立渲染副本。直接改动共享 ggproto 图层环境仍遵循普通 ggplot2
引用语义，不会自动隔离为一次独立编辑。
