<div class="task-intro"><p>在保持注释几何的前提下安排文字：先选输出尺寸和标签子集，再选择自动或固定位置。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。

{{prepare:setup,sequencepalette,docstyle}}

## 选择布局 {#layout}
基因标签提供 `"auto"`、`"radial"`、`"arc"` 布局。通用注释标签还支持优先
内部放置的默认 `"feature"` 和外部 `"callout"`。应按实际图层选择已记录的取值，
这些布局并不是可在所有入口互换的预设。



{{code:labels}}
{{figure:labels}}

## 固定与自动放置 {#fixed}
手动 `label_rotation`、径向/周向位移和方向控制使用 `geom_gene_label()` 或
`geom_feature_label()`。自动层拒绝这些手动控制，应选择 `label_layout` 或删除手动参数。
固定文字样式使用受支持的 `size`、`colour`、`family`。

在支持的入口中，`label_fit`、`label_wrap`、`label_max_lines`、`label_min_size`
控制适配。缩小字号须显式选择且有下限。无法放置的标签可能移位或隐藏，注释
几何仍保留。应在最终导出尺寸查看 `label_fit_diagnostics`、`label_diagnostics`。

## 跨层匹配锚点 {#anchors}
共线性模式中，省略的注释标签 Position 可继承无歧义的匹配注释层。
唯一 `feature_id` 允许重排或取子集，但标识、区间、链仍须一致。
多个匹配层需要明确 Position 或无歧义 ID。先选最终设备尺寸再预览；调整 Plots
面板会改变可用空间。

## 对照三种基因标签布局 {#comparison}
{{code:labellayouts}}
{{figure:labelauto}}
{{figure:labelradial}}
{{figure:labelarc}}

每条序列最多选最前和最后两个有信息的标签，保留全部基因箭头，便于比较文字方向。
打印 `label_plots$auto`、`$radial` 或 `$arc` 查看候选图。三图输入和请求字号相同，
仅布局选择不同。auto 在适用处使用共享侧列，radial 跟随局部轮廓，arc 跟随序列切线。
`p_manual_labels` 展示独立固定标签入口，显式旋转不属于自动 repel 图层。
通用注释标签还可选内部优先/callout，密集混合注释见[完整圆图](plasmids.html)，
继承或显式 Position 锚点见[堆叠](annotations.html#stack)。

## 完成后检查 {#check-result}
三张布局对照使用同一输入。判断可读性时同时看文字方向、连线和隐藏原因；没有标签不等于注释丢失。

<div class="next-step">继续：<a href="output.html">预览和保存图形 →</a></div>
