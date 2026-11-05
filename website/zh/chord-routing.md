先用控制点改变带的走向，再用同一条比对带比较 CDS 附近的端点间距。这两种调整分别作用于路径中段和端点附近。

## 改变连接带经过的位置 {#control}

{{code:frontbackbone,frontlinks,frontribbon}}

{{figure:frontribbon}}

默认控制点在布局原点。将 `link_ctrl_point` 设为 `c(0, 0.8)`，让控制点在布局坐标中上移。对比带在图中央的走向，两端仍对应原来的区间。

{{code:frontroute}}

{{figure:frontroute}}

控制点用旋转前的布局坐标。默认旋转为 0°，所以 (0, 0.8) 沿图的正上方牵引连接带。改变整图旋转后，控制点也随图一起转动。贝塞尔曲线由控制点牵引，曲线一般不经过该点。

## 用一个真实 CDS 比较避障 {#clearance}

取 `MT108731.1` 上的 48,053–49,111 bp CDS。它在这条比对的查询区间内。这里明确把正链箭头向内移 `0.05` 个布局单位，使箭头落在连接带一侧。

{{code:frontobstacle,frontnoavoid}}

{{figure:frontnoavoid}}

`link_avoid = "none"` 使用固定端点间距。沿`MT108731.1` 骨架靠近 49 kb 的位置看，比对带会盖住部分箭头。

{{code:frontsmooth}}

{{figure:frontsmooth}}

`link_avoid = "smooth"` 在 CDS 附近平滑增大间距。沿同一段`MT108731.1` 骨架看，带的边缘局部缩进去，给箭头腾出空间。默认的 `"uniform"` 会让同一个端点区间采用一致的避让间距。如果指定数值 `link_gap`，就直接使用该间距。

## 用线连两个位置，用色带标区间 {#point-markers}

{{code:frontmarkers}}

{{figure:frontmarkers}}

外侧色带标出查询序列的完整比对区间。 `region_side = "outside"` 将它放在骨架外侧。直线只连接两个区间的中点，`arrow` 指向目标端。需要看两端各覆盖多长时，使用比对带。

下一节[放置序列和基因标签](chord-labels.html)。

## 分支、筛选与固定样式的边界 {#branch-contract}

| 参数 | 当前默认 | 什么时候改变与注意事项 |
| --- | --- | --- |
| `link_branch` | `"none"` | query/subject 让相同端点共享主干；端点不相同不能凭空组成同源分支 |
| `link_branch_fraction` | `0.15` | 主干相对最短成员的长度比例，范围 (0, 0.5]；仅启用分支时有意义 |
| `pairs` | `NULL` | 显式选择序列对；不创建新比对记录 |
| `link_gap` | `NULL` | 交给避障策略；显式数值覆盖自动间距 |
| `link_avoid` | `"uniform"` | none 固定间距，smooth 局部避让；不避让文字与引线 |
| `link_ctrl_point` | `NULL` | 弦图默认原点；共线性依据端点推导；不是 bp 坐标 |
| `nonadjacent_linewidth` | `0.4` | 只控制共线性非相邻显示行边框，0 沿用普通 linewidth |

固定 `fill/alpha/colour/linewidth` 放在 `aes()` 外；数据差异映射到 ribbon 专用属性并配匹配标度。显式固定值优先于相应映射，所以映射透明度时省略固定 alpha。改变颜色不能证明不同生物学关系。完整签名见 [`geom_link_ribbon()`](reference/geom_link_ribbon.html)。

接着用[共享端点分支与点连接](chord-branches.html)比较分支主干、主干比例及曲线/直线。

## 区间高亮的侧向、宽度与偏移 {#region-parameters}

从第一条真实 HSP 的 query 端保留原始 ID 与区间。默认宽度 0.08、偏移 0、侧向 inside；后图分别只改变 outside、宽度或偏移。它们是局部几何量，不是 bp。固定填色写在 geom 外观参数，数据类别映射使用 region_fill 和其专用标度。不能把高亮当作新的功能注释。

{{code:comparefinishchord}}
{{code:compareregioninside}}
{{figure:compareregioninside}}
{{code:compareregionoutside}}
{{figure:compareregionoutside}}
{{code:compareregionwidth}}
{{figure:compareregionwidth}}
{{code:compareregionoffset}}
{{figure:compareregionoffset}}
