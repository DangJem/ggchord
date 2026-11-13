{{code:courseblockbr}}

{{prepare:fusiondata,fusionpanels,fusionlinear,fusionview}}

## 源区间与独立局部细节 {#magnify}

`magnify_ggchord()` 组合一张 Circular/Chord/Collinear 源图和显式构建的 Collinear 细节图。源区间是生物学数据，不是图片像素中的矩形；同一区间须在两个面板完整可见。两个面板都按最终设备尺寸重建，不拉伸截图。

下例选择 pETDuet-1 的原始位置 101–1900，共 1800 bp。5404→2 的跨原点注释与该窗口没有交集；进入线性聚焦前显式排除，不能用 min/max 把它误画成长特征。细节保留 10 条注释，并按明确规则给 CDS/复制起点子集加标签。反转细节的显示方向不改变源坐标或生物学链。

{{code:courseblockbs}}

![单一来源局部细节](../assets/coursefocusedview.png)

用普通 `ggsave()` 指定物理尺寸保存，例如 12×8 英寸；`view_ggchord()` 使用同一渲染器。`center_size`、`detail_size` 是容器比例，`highlight_width` 则以物理 mm 计。细节方位和显示方向会影响引线是否交叉，这里不承诺全局避让。`connector = "none"` 只保留高亮。

虚线是装饰性的聚焦引线，与生物学 `connections`、`relationship_data` 分开导出；添加真实 `fusion_ribbon()` 仍须提供真实关系数据。`focus_region_data` 保留原始选区及任意元数据，`focus_regions` 记录解析出的实例、源原点和采样物理路径长度比。细节容器较小时，该比值可能小于 1，这不改变生物学长度。

细节图拥有自己的数据、标度、图层和统计。若使用 count/bin/coverage 轨道，请自行准备并说明局部或全局统计总体；组合器只记录 `explicit_detail_plot`，不会静默选择口径。它不复制任意源图层、不推断比对，也不接收未注册的 Cartesian 字符面板。当前使用一个有序区间；同一 accession 重复时显式指定 `instance`/`detail_instance`，跨原点选区按后文显式声明 `wrap = TRUE`。下节提供独立细节；后文进一步介绍显式嵌套与自由区域分配，实测障碍路由见后文。
## 保留真实 HSP 的成对源窗口 {#paired-detail}

本例的 `region` 每行选择一个原始源区间，所有行共享一个显式细节面板。本例从四序列数据中选取一条真实 HSP 的两个端点窗口，基因方向和局部 HSP 均来自聚焦后的输入表。派生区域、前缘和装饰引线通过 `region_row` 对应输入行；任意原始元数据保留在 `focus_region_data`。

{{code:courseblockbt}}

![HSP 两端共享细节图](../assets/coursepaireddetail.png)

这仍然只有一个细节面板，不推断逐碱基比对。每个区间须在两个面板完整可见；实例选择按行独立解析。
## 独立细节面板 {#multiple-details}

将一条真实 HSP 的两端分别绘为独立局部图。命名列表提供细节图，`detail_panel`
为每个源行选择目标，命名 `side` 向量指定不同侧边。每个面板拥有自己的数据和标度。

{{code:courseblockbu}}

![两个独立来源细节面板](../assets/coursemultiple.png)

本例使用四侧命名细节，每侧一个，每个面板须有选区；缺失或未知目标会报错。
源/目标实例选择相互独立。`detail_size` 接受按面板命名的宽高比例列表，NULL 使用
Fusion 默认分配；区域重叠会报错。标准保存按目标尺寸重建面板和引线。小画布需要
检查标签并适当增大尺寸。引线可能穿过主图或相互交叉，全局避障是独立策略；显式嵌套与矩形区域分配见后文。
## 精确关联源坐标的字符细节 {#character-detail}

普通 Cartesian 字符图通过 `register_ggchord_alignment()` 可成为独立细节面板。注册实际显示的字符行：x为原始对齐列/编码列，y为对齐行ID，label为真实字符。源映射与普通 ggplot 标度保持独立。以下真实60列片段保留上游 gap 和 mismatch；宏观长度表示**去gap后的对齐片段**，不是完整基因组。

{{code:courseblockbv}}

![已有真实字符的注册细节](../assets/coursecharacterdetail.png)

普通图仍可编辑，各自颜色、图层与注释保持独立。支持 identity/reverse x 标度和离散 y 行。Fusion 每次构建都核对修改后的字符位置/标签，并重新测量最终面板毫米边界；半格边界对应真实端点字符。gap 无源位置，不会被插值为碱基。选区必须按生物学顺序包含全部已选源碱基，内含子需拆成独立源区间，跨原点区域使用后文的显式环绕声明。同一来源有多个对齐行时用 `detail_instance` 明确选择。不注册的 ggplot、分面、非线性/翻转坐标和移动已注册字符均报错。

`translate_ggchord_cds()` 的碱基行可以同样注册。按 `examples/m2-detail.R` 将翻译标在 codon 中心，完整跨段 codon 来源查 `source_positions`。AA残基行使用独立aa坐标域，基因组bp宏观图直接连接aa位置仍需显式CDS投影。magnify装饰引线与生物学Fusion连接分开；字符端点没有骨架颜色可用于渐变，应明确提供跨面板fill。窄画布仍执行 `geom_residue()` 的tile/error策略，应增加物理尺寸，不能只增加DPI。全局引线避障是独立布局问题；下文提供显式嵌套和自由区域分配。
