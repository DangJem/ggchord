在同一张弦图上先加入一条比对，再加入一个 CDS。代码按本页顺序运行。

## 从四条骨架开始 {#bundled}

<span id="backbones"></span>用[上一节](chord-start.html)的四条序列。运行下方代码即可重现骨架。

{{code:frontbackbone}}
{{figure:frontbackbone}}

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

<span id="labels"></span>如果需要显示 CDS 名称，可以加上 `geom_gene_label(data = chord_genes_example)`，标签来自 `anno` 列。下一节用相同骨架比较[旋转与间隙](chord-tuning.html)。
