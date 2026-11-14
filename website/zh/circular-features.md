{{code:courseblockak}}

{{prepare:comparefinishcircle}}

## 通用形状与质粒预设 {#feature-shapes}

先用同六条真实注释只改变 arrow/block。block 不画方向箭头，但原始 strand 没有改变。第三幅展示 geom_feature_plasmid 的组合预设，默认 width=0.058、arrow_head_length=0.068、arrow_head_width=1.42，并使用质粒 Position；这会同时改变多个默认，不能当作只变一个参数。其他形状用途见[形状表](annotations.html#shapes)，箭头细节与短特征见[几何课](feature-geometry.html)。

{{code:comparefeaturearrow}}
{{figure:comparefeaturearrow}}
{{code:comparefeatureblock}}
{{figure:comparefeatureblock}}
{{code:comparefeaturepreset}}
{{figure:comparefeaturepreset}}

## 让重叠特征与文字共享位置 {#stacking}

下面用同一六条实际注释比较 identity 与共享的 feature stack。Position 同时传给几何和文字，不单独挪动文字假装注释分层。`position_feature_stack(spacing = 0.1, side = "strand", base_position = NULL)` 先应用可选基础 Position，再分配重叠轨道。`position_plasmid()` 是局部预设；`position_strand()` 按真实链应用带符号位移。检查方向时不要用内外侧代替生物学链。

{{code:comparefinishstack}}
{{code:comparepositionidentity}}
{{figure:comparepositionidentity}}
{{code:comparepositionstack}}
{{figure:comparepositionstack}}

## 操作示例、替换输入与下一步 {#adapt}

只把双线间距改为 0.04；再把中心文字的 length_style 改为不同字号。保留完整长度，不用注释最大 end 猜长度。用自己的 FASTA/GenBank/GFF3 时核对 ID 和坐标，见[格式导入](formats.html)。接下来了解[圆图标签](circular-labels.html)、[引物](primers.html)、[酶切](restriction-sites.html)和[轨道](tracks.html)。所有预览使用图中尺寸；用同样宽高的 ggsave 保存。
