<a id="intervals"></a>

## 每条实测比对绘制一条连接带 {#intervals}

这三条局部比对连接 MT108731.1 和 MT118296.1。每行有四个端点，带的两端分别覆盖对应骨架上的区间。坐标递减可以表示反向匹配，准备连接表时应保留原始端点顺序。

{{example:paired-ribbons}}

<a id="points"></a>

## 选择点连接、区间连接和序列对 {#points}

连接两个位点时使用 `geom_link_line()`，两端字段为 `qaccver, qpos, saccver, spos`。连接两段区间时使用 `geom_link_ribbon()`。示例固定填色，使三条记录同色。要表示一致性百分数，可保留输入的 `pident`，使用连接带默认的一致性映射。

共线性图默认比较相邻行。把关心的序列对表传给 `pairs`，或用 `adjacent_only = FALSE` 纳入非相邻行。解读缺失的连接前，可以先核对[来源行](diagnostics.md)。
