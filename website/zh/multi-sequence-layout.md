序列布局决定参考序列如何显示。先从 MT108731.1 和 MT118296.1 的两行对照图开始，再尝试调换行序、重复显示一条参考序列，或将骨架画成等长。

<a id="direction"></a>

## 翻转一条序列的显示方向 {#direction}

在[起始布局](collinear.md#adapt)中，MT108731.1 位于上方，两个箭头都指向右侧。把 MT118296.1 排在首位，再用具名向量将 MT108731.1 的方向设为 −1。现在较长的骨架在下方，箭头朝左，三条 HSP 仍沿用原始坐标。

{{example:order-orientation}}

<a id="length"></a>

## 用等长骨架比较相对位置 {#length}

只看第一条 HSP，将 `length_mode = "equal"`。两行现在有相同的显示宽度，但末端刻度分别仍为 64,323 bp 和 32,090 bp。这条 HSP 在目标序列上覆盖的相对比例更大。

{{example:layout-equal}}

<a id="instances"></a>

## 指定哪一份骨架接收连接带 {#instances}

在顺序向量中给三行分别取名为 `top`、`middle` 和 `bottom`。上下两行都是 MT108731.1。第一层把最长的 HSP 连到上方骨架，第二层把另外两条连到下方骨架。每一层都用 `between` 指定该组记录连接的两行。

{{example:layout-instances}}

下一步可[比较区间连接带和中点连线](multi-sequence-links.md)。要调整间距、方向或逐序列几何参数，请查阅[灵活参数](parameters.md)和[序列布局参考](../../man/layout_sequences.Rd)。
