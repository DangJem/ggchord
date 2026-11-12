## 连接真实来源区间 {#intervals}

Fusion组合独立的基因组图形，并连接其真实序列区间。中央可用Circular或Chord，外围为横向或纵向Collinear。每张子图仍可单独使用；组合层先在实际物理尺寸内重新构建各面板，再绘制连接。

从[输入下载与来源页](formats.html)取得原始EMALE PAF。本课的数据整理需要可选依赖`dplyr`。PAF包含11条已观测HSP；下面显式选择长度≥900、Identity≥85%的两条用于跨面板展示。计数轨道按明确的2,000 bp分箱统计全部原始query起点观测，不解释为覆盖深度。最后一个显示分箱明确裁切到真实序列边界，来源观测及计数保留。

{{code:fusiondata}}

保留`read_paf()`规范化后的反向端点。构造器对降序端点的提示不意味着应将它们排序。先准备两张独立图形；明确指定序列配色，使相同accession具有一致颜色。本例主动隐藏局部图例，因为序列颜色已在骨架显示，只需保留跨面板Identity图例。

{{code:fusionpanels}}

{{code:fusionlinear}}

固定的`qpanel`与`spanel`选择准确的query和subject面板。Ribbon落在多序列面板中匹配序列的具体区间，不是面板外框或中心。query/subject角色不表示遗传物质的真实生物流动。

{{code:fusionview}}

{{figure:fusionview}}

## 注册面板与显示实例 {#panels}

每个方向最多一个面板：`right`、`top`、`left`、`bottom`。ID必须唯一，`center`是保留名称。`size=c(width,height)`与`center_size`使用扣除跨面板图例后组合区域的相对比例。明确指定的区域发生重叠时会报错，不擅自移动面板。
方位只决定摆放位置，不旋转文字；需要纵向序列轴时，在子图上设置`coord_collinear(seq_axis="vertical")`。

四方向示例复用真实PAF第3、2、9、1行，明确连接到上、左、下、右面板。重复视图是不同显示实例，不是新增观测。不同面板无需共享相同的物理bp/mm比例。

{{code:fusionfour}}

{{figure:fusionfour}}

选择信息也可写入`qpanel`、`spanel`、`q_instance`、`s_instance`列。同一panel选择不得同时写固定参数和同名列。如果accession出现在多个面板或重复实例中，只有完整可见区间的候选唯一时才允许省略选择；有歧义会报错，必须指定目标面板和实例。对于明确改名的focus窗口，仍可通过保留的来源原点使用原始accession和源坐标。
重复展示`source_hsp_id`时，每个参与图层均须设置`allow_repeated=TRUE`。没有ID时按各层源行记录，不因值相同就猜测为同一生物学身份。

## 样式与物理布局 {#styles}

跨面板层接受`aes(fill=..., colour=..., alpha=..., linewidth=..., linetype=...)`及标准ggplot2标度。子图的角色标度只在各自面板生效；明确指定的固定样式优先于对应映射。默认Identity范围0–100，按0/50/80/90/95/100分段并保持数值比例。端点渐变与subject尖头须显式开启，复用单图Ribbon绘制器。尖头是装饰性方向标记，仍导出完整来源前缘。

{{code:fusionstyle}}

默认`layout="fixed"`分配面板区域，并按输出尺寸重新测量标签和图例；`layout="adaptive"`还按物理宽高比调整中央区域。用户明确指定的字号保持固定。这不是通用字号缩放或全局碰撞求解器。
每张子图保留自己的图例，不自动合并含义不同的标度；需要简化时，主动设置`show.legend=FALSE`或局部`theme()`。组合层的`theme()`/`guides()`只作用于跨面板图例。必须保留的图例空间会实际测量；画布过小则给出错误诊断。

## 导出与核对锚点 {#export}

{{code:fusionexport}}

`panels`记录实际面板与容器的毫米边界；`connections`保留原始关系列、来源行和最终路由，方向未知的行会明确省略并报告。`fronts`包含局部坐标、来源区间边界、全局毫米坐标及局部切向/法向；`paths`补充连接几何。
`source_layouts`、`source_data`与`relationship_data`保留各自来源。诊断包含局部标签可见性、字号拟合及图例尺寸，不将未测量的全局碰撞写成零。

区间为一基闭区间，完整边界使用`[min(start,end)-1,max(start,end)]`，并按原始端点方向排列。单位置区间仍有非零宽度。这个新增组合契约不重写旧单图Ribbon几何。源位置、区间边界和毫米坐标是不同的量。

`export_ggchord_fusion()`在临时设备上测量并恢复此前设备；保存时会在实际输出设备上重建，因此字体度量可能略有不同，但源坐标不变。使用标准`ggsave()`或共享预览；不要对组合对象调用`ggplot_build()`。

{{code:fusionsave}}

要求区间完整可见；需要局部窗口或拆分时，先在来源工作流中明确准备。不支持跨原点环绕解释、直接bp↔aa连接、外围Sankey/Radial/Dotplot、分面、猜测外部ggplot坐标、自由方位/八方向布局或全局自动路由。
连接可以交叉或穿过面板内部；前景注释保留，但不保证任意组合都无碰撞。


## 源区间与独立局部细节 {#magnify}

`magnify_ggchord()` 组合一张 Circular/Chord 源图和显式构建的 Collinear 细节图。源区间是生物学数据，不是图片像素中的矩形；同一区间须在两个面板完整可见。两个面板都按最终设备尺寸重建，不拉伸截图。

下例选择 pETDuet-1 的原始位置 101–1900，共 1800 bp。5404→2 的跨原点注释与该窗口没有交集；进入线性聚焦前显式排除，不能用 min/max 把它误画成长特征。细节保留 10 条注释，并按明确规则给 CDS/复制起点子集加标签。反转细节的显示方向不改变源坐标或生物学链。

```r
plasmid <- plasmid_example_pETDuet_1
annotations <- find_common_features(plasmid)
selected_region <- data.frame(accver = plasmid$accver, start = 101, end = 1900)
local_source <- focus_ggchord_data(
  plasmid,
  feature_data = annotations |> filter(start <= end),
  loci = selected_region
)
whole_plasmid <- ggchord(plasmid, feature_data = annotations) +
  geom_seq(show.legend = FALSE) +
  geom_feature(position = position_plasmid(), show.legend = FALSE) +
  geom_seq_center_label() +
  coord_circular(rotation = -30)
local_labels <- local_source$feature_data |>
  filter(type %in% c("CDS", "rep_origin"), anno != "ATG")
local_detail <- ggchord(
  local_source$seq_data, feature_data = local_source$feature_data
) +
  geom_seq(seq_orientation = -1, show.legend = FALSE) +
  geom_feature(show.legend = FALSE) +
  geom_feature_label(
    data = local_labels, label_orientation = "horizontal",
    label_overlap = "allow", colour = "#202020", size = 3
  ) +
  coord_collinear(seq_axis = "vertical") +
  scale_seq_position_continuous(coordinates = "source")
focused_view <- magnify_ggchord(whole_plasmid, selected_region, local_detail)
focused_view
focus_geometry <- export_ggchord_fusion(focused_view, width = 12, height = 8)
focus_geometry$focus_regions
focus_geometry$leaders
```

用普通 `ggsave()` 指定物理尺寸保存，例如 12×8 英寸；`view_ggchord()` 使用同一渲染器。`center_size`、`detail_size` 是容器比例，`highlight_width` 则以物理 mm 计。细节方位和显示方向会影响引线是否交叉，这里不承诺全局避让。`connector = "none"` 只保留高亮。

虚线是装饰性的聚焦引线，与生物学 `connections`、`relationship_data` 分开导出；添加真实 `fusion_ribbon()` 仍须提供真实关系数据。`focus_region_data` 保留原始选区及任意元数据，`focus_regions` 记录解析出的实例、源原点和采样物理路径长度比。细节容器较小时，该比值可能小于 1，这不改变生物学长度。

细节图拥有自己的数据、标度、图层和统计。若使用 count/bin/coverage 轨道，请自行准备并说明局部或全局统计总体；组合器只记录 `explicit_detail_plot`，不会静默选择口径。它不复制任意源图层、不推断比对，也不接收普通 Cartesian 字符面板。当前使用一个有序区间；同一 accession 重复时显式指定 `instance`/`detail_instance`，跨原点区间须先拆分。多个或嵌套 inset 盒及其他中央 Coord 不属于此接口。
