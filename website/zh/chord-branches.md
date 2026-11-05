连接带表达区间，点连接表达选定位置。本页用最小简化数据比较分支与路径；这些记录仅用于几何示意，不是实际序列比对或同源证据。

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

弦图与共线性图可用这些连接层；单序列圆图拒绝 ribbon/link。共线性图的邻接筛选仍适用，要显示跨行带先检查 `adjacent_only` 或显式 `pairs`，见[跨行连接](collinear-routing.html)。

练习：只将主干比例改为 0.25；再回到默认比例，比较曲线与直线的选点。用自己的数据时先决定需要展示区间长度还是单个位置，随后查看 [`geom_link_ribbon()`](reference/geom_link_ribbon.html) 与 [`geom_link_line()`](reference/geom_link_line.html)。
