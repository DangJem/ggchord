先确定数据含义，再选角色和标度。这里把同一批真实 HSP、注释与 GC 窗口分别用于连续值、类别、identity、透明度、线宽、线型与主题；无需为每个辅助函数创建一页。

## 独立输入与标度职责 {#input}

本页独立加载 ggchord 和可选 dplyr。四条完整序列、HSP identity 与 query/subject ID 来自随包真实记录。命名向量用于把已有类别对应到显示值，不创造新的样本或统计量。映射写 aes，固定样式写 geom，用户 scale 优先于自动默认。显式固定属性会覆盖对应映射，想映射透明度时不要再传固定 alpha。

{{code:comparescaledata}}


## 连续、分段与直接颜色值 {#numeric}

两图都将真实 pident 映射到 ribbon_fill，范围固定 0–100、配色和尺寸固定。gradientn 连续插值；stepsn 按给定 breaks 分段。颜色相同不证明新的同源关系。values 是渐变中的颜色锚点，不是原始 pident；breaks 控制刻度/分段，limits 控制数据域。NA 与越界值应检查标度策略，不当作零。

{{code:comparescalegradient}}
{{figure:comparescalegradient}}
{{code:comparescalesteps}}
{{figure:comparescalesteps}}


## identity 不再替数据选色 {#identity}

下面先从真实 identity 按显式 90% 阈值生成显示颜色，再用 scale_ribbon_fill_identity 直接读取颜色字符串。这是事先声明的显示分类，不是恢复连续 identity 的方法。若需要连续图例，用上一节的渐变。identity 通常不显示图例；要显示时必须明确 breaks、labels 和 guide，使图例解释颜色的业务含义。

{{code:comparescaleidentity}}


## 类别填色、透明度、边框与线型 {#categories}

同一 qaccver 映射到同一显示值。manual 系列的 values 推荐按实际类别命名，避免类别顺序变化后错配。此例同时训练多个独立角色，是组合示范，不是单参数视觉对照。连续透明度用 `scale_ribbon_alpha_continuous(range=...)`；离散透明度用 `scale_ribbon_alpha_manual`。colour 与 linetype 只改变边框；不能用边框外观代替区间方向。

{{code:comparescalemanual}}


## 点连接的独立标度家族 {#point-links}

先从真实 HSP 显式计算中点；不把点连接当作覆盖整个比对区间。manual 系列分别给 link_colour、link_alpha、link_linewidth、link_linetype 映射类别；gradientn、alpha_continuous、linewidth_continuous 对应数值。下面先构造类别组合，再用 pident 训练连续组合。

自动托管的数值 pident 映射使用0–100域和0、50、80、90、95、100六个Identity刻度；显式用户标度继续优先。此例为三个通道设置相同域、刻度与标题，使图例键可以合并。

{{code:comparescalelinkdata}}
{{code:comparescalelinkmanual}}
{{code:comparescalelinkcontinuous}}
{{figure:comparescalelinkcontinuous}}


## 注释、基因与区域角色 {#features}

feature_fill 和 feature_shape 可按同一个 type 映射，形状与填色图例可合并；独立 shape guide 需要显式请求。固定 shape 会优先于形状映射。这里所有 type 先用 arrow，是为了演示映射规则而不赋予未验证的形状功能。质粒 fill/shape 预设使用其记录的类别规则。

{{code:comparescalefeature}}
{{code:comparescalepreset}}


## 基因链、序列 ring 与局部区域 {#rings}

基因链采用原始 source_strand。region_fill 与 gene_fill、ribbon_fill 互不覆盖。seq_ring 必须有显式正半径表，不推断轨道数量；不能再同时传 seq_radius，也不能用于共线性图。它是序列半径映射，不等于定量 position_track 的命名轨道。

{{code:comparescalegenering}}
{{code:comparescaleregion}}


## 命名轨道的值域与颜色 {#tracks}

先用真实 pUC19 DNA 计算窗口。track_value 控制数值高度，track_fill/track_colour 训练各自颜色；track ID 必须和 position_track 一致。下面的分类阈值 0.5 是显示选择，完整比例值仍保留。tile 不自动求均值，不把 NA 高度变成零。多个 ID 的颜色和数值标度分别声明。

{{code:trackdata}}
{{code:comparescaletrack}}
{{code:comparescaletrackgradient}}


## guide、别名与 staged mapping {#guides}

guide_ggchord_legend 管理离散角色图例；guide_ggchord_colourbar 管理连续颜色条。show.legend 控制图层是否参与，guides 选择图例，theme 控制位置、方向与外观。colour/color 拼写别名有相同语义：guide_ggchord_colorbar、scale_seq_color_manual、scale_ribbon_color_manual、scale_link_color_manual、scale_link_color_gradientn 分别对应 colour 版本，不能同时添加两种同义标度来叠加控制。

after_stat 读取统计变量，例如下一课的 bundle_n；after_scale 在标度之后计算受支持的视觉属性。两者来自 ggplot2，即使 R 对象看起来相同，staged mapping 的语义也不同。所有 ggplot2/grid 重导出遵循其上游帮助，参数图谱明确区分归属。



## 主题只表达非数据外观 {#themes}

theme_ggchord、minimal、publication、dark 是共享非数据样式的不同起点；plasmid 适合单序列圆图。下面两个图只有主题不同，骨架与数据映射不变。theme 中的轴间距、刻度长、文字偏移使用文档声明的物理单位；自定义主题字段不是按 accver 命名的几何列表。只关闭文字会隐藏文字，不会删除序列。

{{code:comparethemes}}
{{figure:comparethememinimal}}
{{figure:comparethemedark}}


## 练习与选择规则 {#adapt}

先把 pident 映射从连续渐变改为显式分段，再改回渐变，只改变 alpha 输出范围；核对图例是否仍表示相同数据单位。新增图层时给每层独立 data/mapping，保留已存在的用户 scale。完整参数、类型、默认与转发边界见[参数图谱](parameter-atlas.html)；聚束和密度见[连接统计](ribbons.html#statistics)。
