Sankey 坐标把参考序列排成左右两列。带子左端取 query 区间，右端取 subject 区间，两端宽度分别由记录的跨度决定。

<a id="draw"></a>

## 从左列的 query 追踪到右列 {#draw}

输入的前三行均连接 MT108731.1 和 MT118296.1。在第一条序列的左列副本找到三个区间，沿带子可到达第二条序列的右列副本。每条参考序列均在两列出现。

{{example:first-sankey}}

<a id="adapt"></a>

## 理解宽度与方向 {#adapt}

带宽按来源端点的坐标单位读取。一条比对在两端可以占据不同数量的位置。若输入同时提供互换 query 和 subject 的记录，这些行各自表示一次观测。可用 `column_gap` 调整两列之间的距离。

[四基因组双列图](case-sankey.md) · [Sankey 坐标](../../man/coord_sankey.Rd)
