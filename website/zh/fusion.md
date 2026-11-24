Fusion 组合使用各自坐标系的独立图。这里的圆图展示完整的 MT108731.1 基因组，两条竖直骨架并列展示该基因组和 MT118296.1。一条实测 HSP 将圆图中的源区间连到第二条骨架。

<a id="draw"></a>

## 排列两幅独立图 {#draw}

先构建环形总览和直线详图。将具名列表 `list(overview = overview, detail = detail)` 传给 `ggchord_fusion()`。添加 `layout_fusion_grid(ncol = 2, gap = 3)`，将两个视图排在一行，中间留出 3 mm。列表名称用于标识连接图层中的视图路径。

{{example:fusion-panels}}

<a id="adapt"></a>

## 连接指定视图中的 HSP {#adapt}

继续使用这个查询方向的第一条 HSP。在 `geom_link_ribbon()` 中将 `qpanel_path` 和 `spanel_path` 映射到两个视图的名称。实测点对点关系可用 `geom_link_line()`。查询端连到 `overview` 中 MT108731.1 的 26,298–51,139 bp，目标端连到 `detail` 中 MT118296.1 的 7,121–31,959 bp。参考序列标识指定来源，面板路径指定视图。

{{example:first-fusion}}

<span id="shared-guides"></span>

每个视图须包含指定的源位置。用 `fusion_magnify_auto()` 从序列总览派生窗口，或用 `fusion_magnify_manual()` 在两幅已准备好的视图中标记同一源区间。有真实字符时，Auto 会在面板容得下实测字形时显示残基详图。设 `level = "feature"` 保持注释窗口，或设 `level = "residue"` 指定字符详图。字符详图可用 `read_fasta_lengths(keep_sequence = TRUE)` 保留参考字符，也可从已登记的字符图开始。在 `ggchord_fusion()` 中设置 `guides = "collect"` 可共同排布外置图例。给 `scale_for()` 传入含 `panel_path` 和 `layer_id` 的表，可以对所选视图共同训练视觉标度。下一步可沿[六面板案例](case-fusion.md)从基因组总览读到局部比对，或[按最终宽度导出图形](output.md)。
