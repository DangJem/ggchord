共线性坐标把每条序列放在独立的一行，便于比较比对区间的顺序和方向，也能观察注释与比对在同一骨架上的相对位置。

<a id="draw"></a>

## 比较两条全长序列 {#draw}

序列表提供以 bp 为单位的全长。示例保留输入中指定 query→subject 方向的全部比对，共得到三条带。可沿两条骨架分别读取带子两端的宽度。

{{example:first-collinear}}

<a id="adapt"></a>

## 选择行顺序和显示方向 {#adapt}

有多行序列时，连接层默认显示相邻行间的比对。给图层设置 `adjacent_only = FALSE` 可加入选定的跨行连接。用 `layout_sequences(order = ...)` 设定行序，用具名方向向量翻转指定行，同时保留其来源坐标和链。

[带 CDS 标签的四基因组比较](case-collinear.md) · [线性坐标与锚定](../../man/coord_collinear.Rd) · [序列排序](../../man/layout_sequences.Rd)
