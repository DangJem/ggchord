<div class="task-intro"><p>在保持注释几何的前提下安排文字：先选输出尺寸和标签子集，再选择自动或固定位置。</p></div>

<span id="prepare"></span>

## 选择弦图标签布局 {#layout}
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

{{code:manualgenelabels}}
{{figure:manualgenelabels}}

这里把固定标签旋转为 30 度。固定位置不会自动分散相邻文字，空间不足的标签可能隐藏；本图只显示了一条标签，四个 CDS 箭头仍保留。需要显示更多名称时，对照下面的自动布局。

## 让标签与箭头匹配 {#anchors}
本页的箭头与标签使用同一份 CDS 表，并保留一致的序列名称、区间与真实链。筛选标签时只筛选这份表，不修改源坐标。固定布局可分别设置径向位移与文字旋转。

## 对照三种基因标签布局 {#comparison}
{{code:labellayouts}}
{{figure:labelauto}}
{{figure:labelradial}}
{{figure:labelarc}}

只取四条有名称的真实 CDS，用完全相同的输入比较三种布局。
在 `view_ggchord()` 中换用 `label_plots$radial` 或 `label_plots$arc`，可按相同尺寸查看候选图。三图输入和请求字号相同，
仅布局选择不同。auto 在适用处使用共享侧列，radial 跟随局部轮廓，arc 跟随序列切线。
`p_manual_labels` 展示独立固定标签入口，显式旋转不属于自动 repel 图层。
通用注释标签还可选内部优先/callout，密集混合注释见[完整圆图](plasmids.html)，
继承或显式 Position 锚点见[堆叠](annotations.html#stack)。

## 检查要点 {#check-result}
三张布局对照使用同一输入。判断可读性时同时看文字方向、连线和隐藏原因；没有标签不等于注释丢失。

<div class="next-step">继续：<a href="output.html">预览和保存图形 →</a></div>
