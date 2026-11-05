序列名用来识别整条骨架，基因名用来识别一个区间。两个标签层可以分别调整。

## 把序列名放到内侧 {#sequence-labels}

{{code:frontseqlabels}}

{{figure:frontseqlabels}}

`seq_label_radius = 0.7` 把文字放到对应骨架半径的 70% 处。和[默认图](chord-start.html#draw)比，标签离刻度更远。默认文字沿弧转动，`seq_label_orientation = "horizontal"` 让四条序列名都横排。若想分辨两项设置的作用，每次只改一项。

## 标注一个 CDS {#gene-labels}

{{code:frontgenes,frontgenelabels}}

{{figure:frontgenelabels}}

在绿色的 `OQ646790.1` 骨架上找到红色 integrase 箭头。新标签取自 `anno`。这里用 `label_radial_offset = 0.10` 给文字多留一点径向距离，方便把它与箭头区分开。多个注释相邻时，可用 `geom_gene_label_repel()` 分散文字，参见[密集标签](labels.html#layout)。

继续看[四序列综合案例](chord.html)，把这些调整合在一起。
