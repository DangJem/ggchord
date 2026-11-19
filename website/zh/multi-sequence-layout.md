<a id="direction"></a>

## 翻转一条序列的显示方向 {#direction}

把 MT118296.1 放在第一行，并翻转 MT108731.1。用 accession 给方向向量命名，之后重排时，翻转仍作用于同一条序列。可以对照[原排列](multi-sequence-links.md)观察三条连接带的两端。

{{example:order-orientation}}

<a id="instances"></a>

## 在图中两次显示同一条参考序列 {#instances}

要比较 A–B–A，可为两个 A 命名，例如 `layout_sequences(order = c(left = "A", middle = "B", right = "A"))`。在每个连接层筛选相应的原始记录，再用 `between = c("left", "middle")` 或 `between = c("middle", "right")` 选择显示副本。来源坐标与生物学链仍属于原序列 A，两个副本均保留其完整长度。

<a id="length"></a>

## 选择按比例或等长显示 {#length}

要比较不同序列上的实际距离，保留默认长度比例。要用等长骨架比较特征在各序列上的相对位置，可用 `layout_sequences(length_mode = "equal")`。碱基数或残基数从各自的局部序列轴读取。
