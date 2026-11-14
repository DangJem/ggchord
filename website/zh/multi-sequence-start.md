一课完成骨架、连接、特征和标签。输入来自随包真实记录，按下面顺序运行；需要 ggchord 和可选 dplyr。

## 绘制骨架 {#draw}

按[安装页](installation.html)完成安装后，在全新 R 会话运行本例。

{{code:courseblockac}}

{{code:frontbackbone}}

{{figure:frontbackbone}}

## 读图 {#read}

每条弧代表序列表的一行。用标签找到 `MT108731.1`，再与 `MT118296.1` 比较。它们分别为 64,323 bp 和 32,090 bp，所以第一条弧约为第二条的两倍长。

`geom_seq()` 画出骨架，`geom_seq_label()` 标出序列 ID，`coord_chord()` 把四条骨架排成弧形。

## 试着去掉标签 {#try}

删除 `geom_seq_label()` 那一行并重新运行。四条弧会保留，序列名会消失。

## 加入一条比对带 {#links}

取随包比对表的第一行。它连接四条骨架中的两条。

{{code:frontlinks,frontribbon}}

{{figure:frontribbon}}

沿着同一条带看两端。它在 `MT108731.1` 上覆盖 26,298–51,139 bp，在 `MT118296.1` 上覆盖 7,121–31,959 bp。两端都有跨度，因此这里使用 `geom_link_ribbon()`。

默认图用颜色显示该行的比对一致性（约 98.6%）。带的交叉形状来自两个区间端点的连接方式，仅凭它不能判定倒位。

## 加入一个 CDS {#add-genes}

<span id="prepare"></span><span id="inspect-genes"></span><span id="genes"></span>选出 `OQ646790.1` 上从 43,765 bp 开始的原始注释，并取回它的真实链方向。

{{code:frontgenes,frontfeature}}

{{figure:frontfeature}}

绿色的 `OQ646790.1` 上出现红色箭头。它是 integrase 的 CDS，范围为 43,765–46,140 bp，原始注释为正链。连接带和这个 CDS 使用各自的数据。

如果需要显示 CDS 名称，可以加上 `geom_gene_label(data = chord_genes_example)`，标签来自 `anno` 列。下一节用相同骨架比较[旋转与间隙](chord-geometry.html#orientation)。

## 显示特征名称 {#labels}

{{code:frontgenelabels}}

{{figure:frontgenelabels}}

## 切换线性布局 {#collinear}

{{code:courseblockad}}

这一步保持输入与图层，但线性视图默认只连接相邻行。完整四视图对照见[下一课](multi-sequence-views.html)。

## 换成自己的序列 {#adapt}

把 `seq_data_example` 换成含 `accver` 和 `length` 的表即可。列的含义见[基本数据表](data.html#prepare)。接下来[添加连接和特征](multi-sequence-start.html#draw)。
