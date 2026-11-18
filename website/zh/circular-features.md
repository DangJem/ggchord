{{code:courseblockak}}

{{prepare:comparefinishcircle}}

<a id="feature-shapes"></a>

## 通用形状与质粒预设 {#feature-shapes}

先用同六条真实注释只改变 arrow/block。block 不画方向箭头，但原始 strand 没有改变。第三幅展示 geom_feature(preset = 'plasmid') 的组合预设，默认 width=0.058、arrow_head_length=0.068、arrow_head_width=1.42，并使用质粒 Position；这会同时改变多个默认，不能当作只变一个参数。其他形状用途见[形状表](annotations.md#shapes)，箭头细节与短特征见[几何课](feature-geometry.md)。

{{code:comparefeaturearrow}}
{{figure:comparefeaturearrow}}
{{code:comparefeatureblock}}
{{figure:comparefeatureblock}}
{{code:comparefeaturepreset}}
{{figure:comparefeaturepreset}}

<a id="stacking"></a>

## 让重叠特征与文字共享位置 {#stacking}

下面用同一六条实际注释比较 identity 与共享分道。Position 同时传给几何和文字，使锚点和图形采用同一道位。`position_feature_stack(spacing = 0.1, base_position = position_strand(), direction = "auto")` 先确定起始道，再按重叠分道。两链同侧时，用 `position_offset()` 指定共同起点，并独立设置增道方向。两链同侧或分置两侧均保留来源链。

{{code:comparefinishstack}}
{{code:comparepositionidentity}}
{{figure:comparepositionidentity}}
{{code:comparepositionstack}}
{{figure:comparepositionstack}}

<a id="adapt"></a>

## 操作示例、替换输入与下一步 {#adapt}

只把双线间距改为 0.04；再把中心文字的 length_style 改为不同字号。保留完整长度，不用注释最大 end 猜长度。用自己的 FASTA/GenBank/GFF3 时核对 ID 和坐标，见[格式导入](formats.md)。接下来了解[圆图标签](circular-labels.md)、[引物](primers.md)、[酶切](restriction-sites.md)和[轨道](tracks.md)。所有预览使用图中尺寸；用同样宽高的 ggsave 保存。
