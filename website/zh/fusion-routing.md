{{code:courseblockbk}}

{{prepare:fusiondata,fusionpanels,fusionlinear,fusionview}}

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

首例的每个命名方向最多一个面板：`right`、`top`、`left`、`bottom`；斜向、自由区域与自动排布见后文。ID必须唯一，`center`是保留名称。`size=c(width,height)`与`center_size`使用扣除原生外层装饰及跨面板图例后内容区的相对比例。明确指定的区域发生重叠时会报错，不擅自移动面板。
方位只决定摆放位置，不旋转文字；需要纵向序列轴时，在子图上设置`coord_collinear(seq_axis="vertical")`。

四方向示例复用真实PAF第3、2、9、1行，明确连接到上、左、下、右面板。重复视图是不同显示实例，不是新增观测。不同面板无需共享相同的物理bp/mm比例。

{{code:fusionfour}}

{{figure:fusionfour}}

选择信息也可写入`qpanel`、`spanel`、`q_instance`、`s_instance`列。同一panel选择不得同时写固定参数和同名列。如果accession出现在多个面板或重复实例中，只有完整可见区间的候选唯一时才允许省略选择；有歧义会报错，必须指定目标面板和实例。对于明确改名的focus窗口，仍可通过保留的来源原点使用原始accession和源坐标。
重复展示`source_hsp_id`时，每个参与图层均须设置`allow_repeated=TRUE`。没有ID时按各层源行记录，不因值相同就猜测为同一生物学身份。

## 导出与核对锚点 {#export}

{{code:fusionexport}}

`panels`记录实际面板与容器的毫米边界；`connections`保留原始关系列、来源行和最终路由，方向未知的行会明确省略并报告。`fronts`包含局部坐标、来源区间边界、全局毫米坐标及局部切向/法向；`paths`补充连接几何。
`source_layouts`、`source_data`与`relationship_data`保留各自来源。诊断包含局部标签可见性、字号拟合及图例尺寸，不将未测量的全局碰撞写成零。

区间为一基闭区间，完整边界使用`[min(start,end)-1,max(start,end)]`，并按原始端点方向排列。单位置区间仍有非零宽度。这个新增组合契约不重写旧单图Ribbon几何。源位置、区间边界和毫米坐标是不同的量。

`export_ggchord_fusion()`在临时设备上测量并恢复此前设备；保存时会在实际输出设备上重建，因此字体度量可能略有不同，但源坐标不变。使用标准`ggsave()`或共享预览；不要对组合对象调用`ggplot_build()`。

{{code:fusionsave}}

要求区间完整可见；需要局部窗口或拆分时，先在来源工作流中明确准备。跨原点仍需明确拆分，不支持直接bp↔aa连接、Dotplot端点、分面、猜测外部ggplot坐标。实测障碍路由见后文。其他已注册视图及显式区域分配见后文。
连接可以交叉或穿过面板内部；前景注释保留，但不保证任意组合都无碰撞。

## 实测路由与显示通道 {#routing}

`fusion_ribbon(route = "avoid")` 保留两端完整区间及物理离开法线，在实测的其他
叶面板容器、原生图层文字和最终 gtable 图例/装饰单元之间寻找较窄的多边形路线。
`route_width` 单位为毫米，不表示比对长度或数量。默认 `"direct"` 保留已有
Bezier 绘制并报告碰撞，不会暗中筛掉记录。

`route = "bundle"` 配合显式字符列 `route_group`，在同一图层的每组及有序面板对内
共享通道。每个源行仍有独立端面、颜色/透明度及多边形。共享通道的叠画不代表
read support 或新增观测。分组由用户声明，不推断同源簇或关系。运行真实三条 HSP 示例：

{{code:courseblockbl}}

![显式显示通道与障碍](../assets/courserouted.png)

`magnify_ggchord(route = "avoid")` 将相同实测空间契约用于装饰引线。导出的
`leader_paths` 含最终折线顶点，与生物学 `paths` 分开。`routing_corridors` 记录
共享通道及成员源行，`routing_obstacles` 记录物理包围框；连接/引线报告路由状态
及碰撞。每个目标尺寸均重新构建，嵌套分配后还会核对祖先新增障碍。

找不到可行路线时明确报错：可移动或扩大面板、缩小选区、调整声明分组，或使用
`route = "direct"` 查看诊断。祖先覆盖阻断已测子路线时也明确失败。采用确定性的
可见路径，不声称全局最优布局。旋转文字使用保守矩形包围框，曲线字符串采用更大
包络，原生装饰使用最终 gtable 单元边界。不声称任意自定义 grob 和所有端点数据
形状均已进入逐字形障碍模型。主题改变测量值，不改变源坐标或证据。
