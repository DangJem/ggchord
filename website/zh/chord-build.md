先完成每个小步骤，再查看综合图。四步复用同一组数据，直接比较每次新增图层的作用。

## 本页准备 {#prepare}
先解压[教学输入](../assets/tutorial-inputs.zip)。
{{prepare:sequencepalette,docstyle,showcasedata}}

## 第一步：只有序列骨架 {#backbones}
{{code:learnseq}}
{{figure:learnseq}}
`geom_seq()` 画骨架，`geom_seq_label()` 标识序列，`coord_chord()` 选择整体视图；此时不需要连接或基因输入。

## 第二步：加入已有比对 {#links}
{{code:learnlinks}}
{{figure:learnlinks}}
相似度由aes和scale表达；显式保留非相邻连接，让后面替换坐标时仍使用相同的连接数据。

## 第三步：加入精选基因 {#genes}
{{code:learngenes}}
{{figure:learngenes}}
`position_strand()` 将正负链分到骨架两侧。仅展示真实筛选的注释，原表仍保留全部CDS。

## 第四步：给每个展示特征加标签 {#labels}
{{code:learnlabels}}
{{figure:learnlabels}}
同一特征表交给箭头层与repel标签层；标签层不另行筛选。与单独添加图层相比，改名、删行或排序会改变业务含义，应保留明确中间对象。

继续：[1.4 方向、半径、曲率和避障](chord-tuning.html)。
