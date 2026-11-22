同样的五条 pUC19 特征可以回答两个问题。统一填色可显示区间的位置。再把填色映射到记录中的注释类型，就能同时看出每个区间属于哪一类。

<a id="styles"></a>

## 用固定填色观察区间位置 {#styles}

选取两个 CDS、两个引物结合区和一个多克隆位点，直接在特征图层设置 `fill`。长短区间都呈相同的蓝绿色，靠近起点、已经分道的区间也是如此。

{{example:fixed-features}}

<a id="categories"></a>

## 把类型列映射到具名调色板 {#categories}

沿用同样的筛选和位置设置。把 `fill = type` 写进 `aes()`，在 `scale_fill_manual()` 中用三个原始类别值给颜色命名。两个 CDS 使用同一种颜色，两个引物结合区也各自同色，第三种颜色单独标出了多克隆位点。

{{example:category-fill}}

调色板中的名称使重新排序或筛选后的类别仍对应原来的颜色。要调整字体或移动图例，可使用 `theme()`。调整不同序列的半径、偏移或宽度，请继续看[灵活参数](parameters.md)。要用高度表示实测值，可试画 [GC 轨道](tracks.md)。

序列标识可用 `scale_colour_sequence()` 着色，链方向用 `scale_fill_strand()`，注释类型用 `scale_fill_feature_type()`，比对一致性百分比用 `scale_fill_pident()`。当基因和连带都映射 `fill` 时，给基因层指定 `id = "genes"`，再加上 `scale_for("genes", scale_fill_strand())`。基因图例显示链方向，连带保留自己的百分比色标。

默认画布为暖白色。需要纯白绘图区时可使用 `theme_ggchord_publication()`。连带默认带有细绿色描边，在连带图层设置 `colour = NA` 即可隐藏描边。
