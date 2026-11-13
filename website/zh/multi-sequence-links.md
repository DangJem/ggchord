{{code:courseblockag}}

{{prepare:frontbackbone,frontlinks}}

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

下一节[放置序列和基因标签](multi-sequence-labels.html#sequence-labels)。
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

接着用[共享端点分支与点连接](multi-sequence-links.html#input)比较分支主干、主干比例及曲线/直线。

## 跨行连接 {#cross-row}

{{prepare:frontlinearbase,frontlinks,frontlinearlink}}

## 允许跨行连接 {#routes}

{{code:frontlinearcross}}

{{figure:frontlinearcross}}

相邻带默认使用细实线，非相邻行之间的带默认使用较粗虚线。普通线宽由 `linewidth` 控制；跨行线宽由 `nonadjacent_linewidth` 单独控制，默认 0.4。

第二个连接层只接收一条第一行到第三行的比对。`adjacent_only = FALSE` 让这条带可见，默认虚线轮廓提示它跨过中间一行。原来第一、二行的连接仍由第一个层绘制。
## 单独设置跨行边框 {#width}

{{code:frontcrosswidth}}
{{figure:frontcrosswidth}}

`nonadjacent_linewidth = 0.8` 加粗第二个图层中的跨行带，不改变第一个图层的相邻带，也不改变线型。

{{code:frontcrosssame}}
{{figure:frontcrosssame}}

设为 0 时，跨行带使用该图层的普通 `linewidth`，这里是 0.2。它仍是虚线。若使用 `aes(linewidth = ...)` 映射，实际线宽由映射及标度控制。

## 共享分支与点连接 {#branches-and-points}

## 准备共享端点的示意输入 {#input}

{{code:comparebranchdata}}

三条示意序列各为 1000 bp。两条连接共用 query 的 200–600 bp，分别到两个 target 的不同区间。所有坐标为一基闭区间。每张图固定 9 × 7 英寸；不为改变分支重写端点。
## 独立带与共享 query 主干 {#branches}

{{code:comparebranchnone}}
{{figure:comparebranchnone}}

默认 `link_branch = "none"` 时两条带各自走向目标，query 端附近可以重叠。

{{code:comparebranchquery}}
{{figure:comparebranchquery}}

只改成 `"query"`，完全相同的 query 端点可形成共享主干。仅仅区间有重叠不等于端点完全相同；真实数据若没有合适组，启用分支也可能不改变图。`"subject"` 则以相同 subject 端点分组。
## 改变主干比例 {#fraction}

{{code:comparebranchfraction}}
{{figure:comparebranchfraction}}

与 query 分支图相比，仅把 `link_branch_fraction` 从默认 0.15 改为 0.35。主干相对最短成员变长；它不是 identity、相似性概率或新的比对长度。

| 参数 | 默认 | 作用与限制 |
| --- | --- | --- |
| `link_branch` | `"none"` | none/query/subject；选择共享哪一端 |
| `link_branch_fraction` | `0.15` | 范围 (0, 0.5]；none 时不产生共享主干 |
| `link_ctrl_point` | `NULL` | 分支时影响远端形状，接合点由主干确定 |
| `pairs` | `NULL` | qaccver/saccver 表选择已有序列对，不生成比对 |
## 用点代替区间 {#point-links}

点连接要显式提供 `qpos/spos`，不能把区间表直接当作点表。本例选择两个区间的中点；这项选择会丢失带宽所表达的区间长度，因此先说明选点依据。

{{code:comparepointdata}}
{{code:comparelinkcurve}}
{{figure:comparelinkcurve}}

默认 `link_type = "curve"` 用曲线连接位置。

{{code:comparelinkstraight}}
{{figure:comparelinkstraight}}

只改为 `"straight"`；位置、输入和尺寸相同。需要有向箭头时另加 `arrow = grid::arrow()`，其方向表达 query 到 subject 的绘图约定，不能单独作为基因链或转录证据。
## 真实数据、模式与练习 {#adapt}

真实比对必须保留自己的 q/s 端点。若没有完全相同的端点，不要为了形成主干而改写坐标。筛选与聚束见[连接处理](ribbons.html)，几何分支不会自动完成这些统计操作。

弦图与共线性图可用这些连接层；单序列圆图拒绝 ribbon/link。共线性图的邻接筛选仍适用，要显示跨行带先检查 `adjacent_only` 或显式 `pairs`，见[跨行连接](multi-sequence-links.html#routes)。

练习：只将主干比例改为 0.25；再回到默认比例，比较曲线与直线的选点。用自己的数据时先决定需要展示区间长度还是单个位置，随后查看 [`geom_link_ribbon()`](reference/geom_link_ribbon.html) 与 [`geom_link_line()`](reference/geom_link_line.html)。

{{prepare:ribboninputs}}

## 选择区间区间连接或点连接 {#interval-points}

{{code:pointlinks}}

{{figure:pointlinks}}

点连接只画前三条真实比对的中点，箭头指向目标端。需要显示整个比对区间时，继续用连接带。

