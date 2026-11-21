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

声明三个全长显示窗口，分别将 `instance_id` 设为 `top`、`middle` 和 `bottom`。最上和最下的窗口来自同一条序列 MT108731.1。连接表中的 `qinstance` 和 `sinstance` 列分别指定两端的显示实例。最长的 HSP 从上方副本连出，另外两条从下方副本连出。

{{example:layout-instances}}

下一步可[比较区间连接带和中点连线](multi-sequence-links.md)。要调整间距、方向或逐序列几何参数，请查阅[灵活参数](parameters.md)和[序列布局参考](../../man/layout_sequences.Rd)。
