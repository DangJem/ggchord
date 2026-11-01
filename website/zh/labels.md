## 选择布局 {#layout}
基因标签提供 `"auto"`、`"radial"`、`"arc"` 布局。通用 feature 标签还支持优先
内部放置的默认 `"feature"` 和外部 `"callout"`。应按实际图层选择已记录的取值，
这些布局并不是可在所有入口互换的预设。

先运行[共享准备代码](index.html#quick-start)：

{{code:labels}}
{{figure:labels}}

## 固定与自动放置 {#fixed}
手动 `label_rotation`、径向/周向位移和方向控制使用 `geom_gene_label()` 或
`geom_feature_label()`。自动层拒绝这些手动控制，应选择 `label_layout` 或删除手动参数。
固定文字样式使用受支持的 `size`、`colour`、`family`。

在支持的入口中，`label_fit`、`label_wrap`、`label_max_lines`、`label_min_size`
控制适配。缩小字号须显式选择且有下限。无法放置的标签可能移位或隐藏，feature
几何仍保留。应在最终导出尺寸查看 `label_fit_diagnostics`、`label_diagnostics`。

## 跨层匹配锚点 {#anchors}
共线性模式中，省略的 feature 标签 Position 可继承无歧义的匹配 feature 层。
唯一 `feature_id` 允许重排或取子集，但标识、区间、链仍须一致。
多个匹配层需要明确 Position 或无歧义 ID。先选最终设备尺寸再预览；调整 Plots
面板会改变可用空间。
