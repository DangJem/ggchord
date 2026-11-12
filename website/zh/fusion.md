## 连接真实来源区间 {#intervals}

Fusion组合独立的基因组图形，并连接其真实序列区间。中央和外围可使用已注册基因组Coord，Sankey/Radial保留明确显示实例，也支持嵌套Fusion和已注册字符细节。每张子图仍可单独使用；组合层先在实际物理尺寸内重新构建各面板，再绘制连接。

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

要求区间完整可见；需要局部窗口或拆分时，先在来源工作流中明确准备。跨原点仍需明确拆分，不支持直接bp↔aa连接、Dotplot端点、分面、猜测外部ggplot坐标。实测障碍路由见后文。其他已注册视图及显式区域分配见后文。
连接可以交叉或穿过面板内部；前景注释保留，但不保证任意组合都无碰撞。


## 源区间与独立局部细节 {#magnify}

`magnify_ggchord()` 组合一张 Circular/Chord/Collinear 源图和显式构建的 Collinear 细节图。源区间是生物学数据，不是图片像素中的矩形；同一区间须在两个面板完整可见。两个面板都按最终设备尺寸重建，不拉伸截图。

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

细节图拥有自己的数据、标度、图层和统计。若使用 count/bin/coverage 轨道，请自行准备并说明局部或全局统计总体；组合器只记录 `explicit_detail_plot`，不会静默选择口径。它不复制任意源图层、不推断比对，也不接收未注册的 Cartesian 字符面板。当前使用一个有序区间；同一 accession 重复时显式指定 `instance`/`detail_instance`，跨原点选区按后文显式声明 `wrap = TRUE`。下节提供独立细节；后文进一步介绍显式嵌套与自由区域分配，实测障碍路由见后文。

## 保留真实 HSP 的成对源窗口 {#paired-detail}

本例的 `region` 每行选择一个原始源区间，所有行共享一个显式细节面板。本例从四序列数据中选取一条真实 HSP 的两个端点窗口，基因方向和局部 HSP 均来自聚焦后的输入表。派生区域、前缘和装饰引线通过 `region_row` 对应输入行；任意原始元数据保留在 `focus_region_data`。

```r
library(dplyr)
seqs <- seq_data_example |> select(accver, length)
features <- gene_data_example |> mutate(type = 'gene')
hsp <- ribbon_data_example |> slice(1)
regions <- bind_rows(
  transmute(hsp, accver = qaccver, start = pmin(qstart, qend), end = pmax(qstart, qend)),
  transmute(hsp, accver = saccver, start = pmin(sstart, send), end = pmax(sstart, send))
)
local <- focus_ggchord_data(seqs, feature_data = features,
  ribbon_data = hsp, loci = regions)
source <- ggchord(seqs, feature_data = features, ribbon_data = hsp) +
  geom_seq(show.legend = FALSE) + geom_gene(show.legend = FALSE) +
  geom_link_ribbon(show.legend = FALSE)
detail <- ggchord(local$seq_data, feature_data = local$feature_data,
  ribbon_data = local$ribbon_data) + geom_seq(show.legend = FALSE) +
  geom_gene(show.legend = FALSE) + geom_link_ribbon(show.legend = FALSE) + coord_collinear()
paired_detail <- magnify_ggchord(
  source + coord_collinear(), regions, detail, side = "bottom"
)
paired_detail
export_ggchord_fusion(paired_detail, width = 12, height = 8)$focus_regions
```

这仍然只有一个细节面板，不推断逐碱基比对。每个区间须在两个面板完整可见；实例选择按行独立解析。

## 独立细节面板 {#multiple-details}

将一条真实 HSP 的两端分别绘为独立局部图。命名列表提供细节图，`detail_panel`
为每个源行选择目标，命名 `side` 向量指定不同侧边。每个面板拥有自己的数据和标度。

```r
library(ggchord)
library(dplyr)
seqs <- seq_data_example |> select(accver, length)
hsp <- ribbon_data_example |> slice(1)
regions <- bind_rows(
  transmute(hsp, accver = qaccver, start = pmin(qstart, qend),
    end = pmax(qstart, qend), detail_panel = "query"),
  transmute(hsp, accver = saccver, start = pmin(sstart, send),
    end = pmax(sstart, send), detail_panel = "subject")
)
make_detail <- function(region) {
  local <- focus_ggchord_data(seqs, loci = region)
  ggchord(local$seq_data) +
    geom_seq(show.legend = FALSE) +
    coord_collinear(seq_axis = "vertical")
}
source <- ggchord(seqs, hsp) +
  geom_seq(show.legend = FALSE) +
  geom_link_ribbon(show.legend = FALSE)
details <- list(
  query = make_detail(regions[1, ]),
  subject = make_detail(regions[2, ])
)
multiple <- magnify_ggchord(
  source, regions, details,
  side = c(query = "left", subject = "right")
)
export_ggchord_fusion(multiple, width = 12, height = 8)$focus_regions
```

本例使用四侧命名细节，每侧一个，每个面板须有选区；缺失或未知目标会报错。
源/目标实例选择相互独立。`detail_size` 接受按面板命名的宽高比例列表，NULL 使用
Fusion 默认分配；区域重叠会报错。标准保存按目标尺寸重建面板和引线。小画布需要
检查标签并适当增大尺寸。引线可能穿过主图或相互交叉，全局避障是独立策略；显式嵌套与矩形区域分配见后文。

## 精确关联源坐标的字符细节 {#character-detail}

普通 Cartesian 字符图通过 `register_ggchord_alignment()` 可成为独立细节面板。注册实际显示的字符行：x为原始对齐列/编码列，y为对齐行ID，label为真实字符。源映射与普通 ggplot 标度保持独立。以下真实60列片段保留上游 gap 和 mismatch；宏观长度表示**去gap后的对齐片段**，不是完整基因组。

```r
library(dplyr)
cells <- read_aligned_fasta("examples/alignment/opuntia-aligned.fasta", "DNA")
rows <- unique(cells$accver)[c(1, 3)]
window <- cells |>
  filter(accver %in% rows, between(alignment_column, 321, 380))
sequences <- cells |>
  filter(accver %in% rows, !is_gap) |>
  group_by(accver) |>
  summarise(length = max(source_position), .groups = "drop")
regions <- window |>
  filter(!is_gap) |>
  group_by(source_accver) |>
  summarise(start = min(source_position), end = max(source_position), .groups = "drop") |>
  rename(accver = source_accver)
letters <- ggplot2::ggplot(window, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), size = 2.7, overflow = "error") +
  ggplot2::scale_fill_hue(l = 85, c = 40, guide = "none") +
  ggplot2::theme_minimal()
letters <- register_ggchord_alignment(letters, window)
macro <- ggchord(sequences) + geom_seq() + coord_collinear()
character_detail <- magnify_ggchord(
  macro, regions, letters,
  center_size = c(.22, .7), detail_size = c(.5, .8)
)
export_ggchord_fusion(character_detail, width = 16, height = 5)$focus_regions
# ggsave("character-detail.png", character_detail, width = 16, height = 5)
```

普通图仍可编辑，各自颜色、图层与注释保持独立。支持 identity/reverse x 标度和离散 y 行。Fusion 每次构建都核对修改后的字符位置/标签，并重新测量最终面板毫米边界；半格边界对应真实端点字符。gap 无源位置，不会被插值为碱基。选区必须按生物学顺序包含全部已选源碱基，内含子需拆成独立源区间，跨原点区域使用后文的显式环绕声明。同一来源有多个对齐行时用 `detail_instance` 明确选择。不注册的 ggplot、分面、非线性/翻转坐标和移动已注册字符均报错。

`translate_ggchord_cds()` 的碱基行可以同样注册。按 `examples/m2-detail.R` 将翻译标在 codon 中心，完整跨段 codon 来源查 `source_positions`。AA残基行使用独立aa坐标域，基因组bp宏观图直接连接aa位置仍需显式CDS投影。magnify装饰引线与生物学Fusion连接分开；字符端点没有骨架颜色可用于渐变，应明确提供跨面板fill。窄画布仍执行 `geom_residue()` 的tile/error策略，应增加物理尺寸，不能只增加DPI。全局引线避障是独立布局问题；下文提供显式嵌套和自由区域分配。

## 八方位与独立物理区域 {#advanced-placement}

`fusion_panel()` 在四个基本方向之外支持 `topright`、`topleft`、`bottomright`、`bottomleft`。默认角落布局为八个方向预留空间，每个侧边名称只使用一次。其他排布可用 `angle`/`distance` 或显式 `box`，显式区域/角度没有固定面板数上限，但每个区域仍须可容纳。位置参数只改变容器，不旋转子图Coord或文字。角度在矩形画布上也表示真实物理方向；distance是分配区域较短边一半的比例，以中心面板的分配区域中心为起点。

```r
data(seq_data_example)
sequences <- seq_data_example |>
  dplyr::select(accver, length)
base_panels <- ggchord(sequences) + geom_seq(show.legend = FALSE)
central <- base_panels + coord_chord()
mirrored <- base_panels + coord_sankey()
free_angle <- ggchord_fusion(central) +
  fusion_panel(mirrored, "mirror", angle = 31, distance = .8, size = c(.16, .22))
free_box <- ggchord_fusion(central, center_box = c(.02, .35, .2, .8)) +
  fusion_panel(mirrored, "mirror", box = c(.45, .98, .1, .9))
physical_box <- ggchord_fusion(central, center_box = unit(c(5, 85, 10, 100), "mm")) +
  fusion_panel(mirrored, "mirror", box = unit(c(100, 210, 10, 100), "mm"))
export_ggchord_fusion(physical_box, width = 12, height = 6)$panels
```

数值box按 xmin、xmax、ymin、ymax给出比例；绝对单位box保持真实物理尺寸，坐标相对于扣除外层连接图例后的组合区域。box不与side/size/angle/distance同时使用，center_box不与center_size同时使用。每次目标尺寸均重新检查物理区域，不暗中夹紧。`ggchord_fusion()` 或 `magnify_ggchord()` 的显式 `allow_overlap = TRUE` 允许有意叠放，后加入的面板盖住先加入者。使用矢量grid裁剪，不截取PNG；透明连接不在同一区域重复叠画。导出的 `overlaps` 表记录覆盖顺序。这允许有意遮罩，不代表自动避障。

Sankey来源同时出现在左右两列，需从 `get_chord_layout(plot)$display_instances` 明确选择 `q_instance`/`s_instance`，歧义时不能省略。Radial、Chord使用各自已注册骨架frame。q/s区间、单位和来源记录保持不变。同一显式 `source_hsp_id` 重复显示时，所有参与的Fusion图层（含嵌套内部）均须设置 `allow_repeated = TRUE`；`source_observations` 区分明确来源ID和显示次数，不通过比较表值猜测观测身份。Dotplot有独立双轴，当前尚未注册为Fusion端点。

## 嵌套源区域保留独立命名空间 {#nested-detail}

细节图可以本身是Fusion。本例选择未改动的pUC19来源上的两个真实坐标窗口，内部主图是较大窗口，内部细节是较小窗口；外层选区通过 `detail/center` 指向内部主图，并非指向图片中的任意矩形。

```r
data(plasmid_example_pUC19)
source_genome <- plasmid_example_pUC19 |>
  dplyr::select(accver, length)
large_region <- data.frame(accver = "pUC19", start = 101, end = 1900)
small_region <- data.frame(accver = "pUC19", start = 301, end = 600)
large_window <- focus_ggchord_data(source_genome, loci = large_region)
small_window <- focus_ggchord_data(source_genome, loci = small_region)
large_plot <- ggchord(large_window$seq_data) +
  geom_seq(show.legend = FALSE) + coord_collinear()
small_plot <- ggchord(small_window$seq_data) +
  geom_seq(show.legend = FALSE) + coord_collinear()
inner <- magnify_ggchord(large_plot, small_region, small_plot)
whole <- ggchord(source_genome) + geom_seq() + coord_circular()
large_region$detail_panel <- "detail/center"
outer <- magnify_ggchord(
  whole, large_region, inner,
  center_box = c(.02, .22, .15, .85), detail_box = c(.27, .99, .1, .9)
)
report <- export_ggchord_fusion(outer, width = 14, height = 7)
report$focus_regions
report$focus_region_data_by_composition
```

`composition_id` 与 `region_row` 将派生选区、front和引线关联到 `focus_region_data_by_composition` 的各层原始选区。顶层 `focus_region_data` 仍是原始外层输入。叶面板用 `detail/center`、`detail/detail` 等路径命名，组合容器不会被当作单一基因组轴。来源在多个叶面板重复时，应指定完整路径与准确实例；`source_panel` 可选择 `center` 下的源叶面板。导出保留实际毫米边界和原始box/angle请求。

各子图保留独立图层、Scale和统计总体，组合器不会暗中复制图层或重算统计。跨原点区间采用后文的显式完整源长度声明。完整负链CDS→基因组focus→真实字符示例见 `examples/m3-nested.R`。字符保留按物理尺寸的tile/error策略，较小嵌套区域不能仅凭增加DPI容纳60个字母。完整多设备发布验收与更广泛的联合布局优化属于独立工作。

## 显式跨原点选区 {#origin-detail}

圆形排版不是环状来源的证据。在 `start > end` 的选区上设置 `wrap = TRUE`。
两个选定叶面板必须保留相同的完整源长度。字符适配器还需显式 `circular = TRUE`，
并在 `source_data` 提供完整字符串或精确 `source_length`；缺少完整字符串时，字符
仍标为未核对。最多覆盖完整一圈，保证源位置到字符列的反查唯一。负链按自身方向
环绕并保留 gap。不会暗中拼接不同的 focus 实例。

完整真实流程见 `examples/origin-detail.R`：从未改动的随包 pUC19 FASTA 提取
82 个实际碱基（2640–2686，再到 1–35）。`focus_region_data` 保留原来的一行；
派生表以 `region_piece` 区分两段，共享同一 `region_row`。只有外侧两个边界有引线，
接缝不是另一条选区或生物学关系。

```r
# Requires dplyr and the source checkout or tutorial input bundle.
source("examples/origin-detail.R")
origin_detail
```

## 实测路由与显示通道 {#routing}

`fusion_ribbon(route = "avoid")` 保留两端完整区间及物理离开法线，在实测的其他
叶面板容器、原生图层文字和最终 gtable 图例/装饰单元之间寻找较窄的多边形路线。
`route_width` 单位为毫米，不表示比对长度或数量。默认 `"direct"` 保留已有
Bezier 绘制并报告碰撞，不会暗中筛掉记录。

`route = "bundle"` 配合显式字符列 `route_group`，在同一图层的每组及有序面板对内
共享通道。每个源行仍有独立端面、颜色/透明度及多边形。共享通道的叠画不代表
read support 或新增观测。分组由用户声明，不推断同源簇或关系。运行真实三条 HSP 示例：

```r
# Requires dplyr and the source checkout or tutorial input bundle.
source("examples/fusion-routing.R")
routed
```

`magnify_ggchord(route = "avoid")` 将相同实测空间契约用于装饰引线。导出的
`leader_paths` 含最终折线顶点，与生物学 `paths` 分开。`routing_corridors` 记录
共享通道及成员源行，`routing_obstacles` 记录物理包围框；连接/引线报告路由状态
及碰撞。每个目标尺寸均重新构建，嵌套分配后还会核对祖先新增障碍。

找不到可行路线时明确报错：可移动或扩大面板、缩小选区、调整声明分组，或使用
`route = "direct"` 查看诊断。祖先覆盖阻断已测子路线时也明确失败。采用确定性的
可见路径，不声称全局最优布局。旋转文字使用保守矩形包围框，曲线字符串采用更大
包络，原生装饰使用最终 gtable 单元边界。不声称任意自定义 grob 和所有端点数据
形状均已进入逐字形障碍模型。主题改变测量值，不改变源坐标或证据。

## 外层 Guide 与面板空间预算 {#outer-guides}

各子图保留独立 Scale 和 Guide。组合层的 Scale 与 `guides()` 描述跨面板连接带；
原生 Guide 可以分别放在不同外侧，无需全部共用同一侧：

```r
library(ggchord)
guide_source <- ggchord(seq_data_example) + geom_seq(show.legend = FALSE)
guide_fusion <- ggchord_fusion(guide_source, layout = "adaptive") +
  fusion_panel(guide_source + coord_collinear(), "detail") +
  fusion_ribbon(ribbon_data_example, aes(fill = qaccver, alpha = pident),
    qpanel = "center", spanel = "detail") +
  ggplot2::guides(
    fill = ggplot2::guide_legend(position = "bottom", nrow = 1),
    alpha = ggplot2::guide_legend(position = "right")
  )
guide_fusion
```

`layout = "fixed"` 保留已有位置；`"adaptive"` 复用原生 Guide 候选/测量引擎，
联合考虑各侧预算和面板区域，再在剩余物理空间重建各子图。可以移动或重排未锁定
Guide；显式位置、方向、行列、字体、key 尺寸和绝对面板框仍是约束。例如宽矮
画布上的默认 Identity 尺可以转为横向，按原来六个标签重新测量长度，不丢掉95或100。

`metadata$guide_layout` 报告选择及锁定项；`guide_bounds_mm`、`guide_sizes_mm`
记录完成图例组合后、测量面板前的最终外层框。嵌套组合的 Guide 报告位于
`diagnostics`。无法满足约束时明确失败。`theme(legend.position = "none")` 沿用
ggplot2 语义隐藏组合图例。组合内部 Guide 暂不支持；各独立子图仍可拥有局部
Guide。这是有界物理空间搜索，不声称全局最优。
