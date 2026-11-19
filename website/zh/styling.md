<a id="categories"></a>

## 按特征类别填色 {#categories}

来源列 `type` 的每种类别对应一种填色。用类别名给调色板命名，排序或筛选后各类别仍能找到原来的颜色。标度同时给出图例标题和图例键。

{{example:category-fill}}

<a id="styles"></a>

## 设置图层或整张图的样式 {#styles}

整层同色时直接设置 `fill = "#2B7580"`。要表示数据差异，将来源列映射到 `aes()` 中的 `fill`，再选择标度。各个起始映射都读取图层源表。如果用 `aes(accver = ref)` 决定特征位置，想按这个显示参考着色，可同时映射 `fill = ref`。

图层设 `show.legend = FALSE` 可隐藏它对所有图例的贡献。点轨道设 `show.legend = c(fill = FALSE)` 则只隐藏这一层对填色图例的贡献。全图字体和图例位置用 `ggplot2::theme()` 调整。数值高度使用[轨道数值标度](tracks.md)，调色板使用填色或颜色标度。
