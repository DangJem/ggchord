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

用 `scale_colour_sequence()` 为参考 ID 配色，用 `scale_fill_strand()` 为链方向配色，用 `scale_fill_feature_type()` 为注释类型配色，用 `scale_fill_pident()` 为百分比一致性配色。基因与连带都映射填色时，为基因层设置 `id = "genes"`，再添加 `scale_for(scale_fill_strand(), layers = "genes")`。使用 `scales_ggchord()` 可为识别到的映射启用传统配色。

默认画布为暖白色。需要纯白绘图区时可使用 `theme_ggchord_publication()`。连带默认带有细绿色描边，在连带图层设置 `colour = NA` 即可隐藏描边。
