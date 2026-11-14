按 API 查找用途、主讲说明、模式、综合案例和完整参数。精确默认值、格式、单位和限制见函数参考；覆盖等级保留在页末维护记录。上游重导出及纯拼写别名按家族说明，不重复制造说明。

<a id="apis"></a>

## API 功能索引 {#apis}

| API 与参考 | 用途 | 用法说明 | 模式/输入边界 | 参数图谱 | 综合案例/工作流 |
| --- | --- | --- | --- | --- | --- |
| [`adapt_ggchord`](../../man/adapt_ggchord.Rd) | 按目标尺寸启用有界主题缩放与图例联合排布 | [用法](styling.md#adaptive) | fixed/scale/adaptive；显式设置优先，不保证全局零碰撞 | [↗](parameter-atlas.md#api-adapt_ggchord) | [↗](styling.md) |
| [`aes`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-aes) | [↗](chord.md) |
| [`after_scale`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-after_scale) | [↗](chord.md) |
| [`after_stat`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-after_stat) | [↗](chord.md) |
| [`annotate`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-annotate) | [↗](chord.md) |
| [`arrow`](../../man/grid-helpers.Rd) | 物理单位与方向箭头 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-arrow) | [↗](chord.md) |
| [`as_ggchord_alignment`](../../man/as_ggchord_alignment.Rd) | 适配已有对齐字符串，保留精确且可缺失的源位置 | [用法](aligned-sequences.md#aligned-source) | 已有DNA/RNA/AA对齐列；Cartesian局部单元与显式源映射 | [↗](parameter-atlas.md#api-as_ggchord_alignment) | [↗](formats.md) |
| [`as_ggchord_feature`](../../man/as_ggchord_feature.Rd) | 将区间转换为 ggchord 通用注释 | [用法](import-sequences.md#adapters) | 表格/受支持对象适配；完整长度与 ID 不得猜测 | [↗](parameter-atlas.md#api-as_ggchord_feature) | [↗](formats.md) |
| [`as_ggchord_link`](../../man/as_ggchord_link.Rd) | 将成对基因组区间转换为 ggchord 连接 | [用法](import-sequences.md#adapters) | 表格/受支持对象适配；完整长度与 ID 不得猜测 | [↗](parameter-atlas.md#api-as_ggchord_link) | [↗](formats.md) |
| [`as_ggchord_seq`](../../man/as_ggchord_seq.Rd) | 将序列信息转换为 ggchord 输入 | [用法](import-sequences.md#adapters) | 表格/受支持对象适配；完整长度与 ID 不得猜测 | [↗](parameter-atlas.md#api-as_ggchord_seq) | [↗](formats.md) |
| [`build_ggchord_cds`](../../man/build_ggchord_cds.Rd) | 建立显式的剪接 CDS 坐标映射 | [用法](projection.md#boundaries) | 显式 CDS/相位/单位的编码坐标工作流 | [↗](parameter-atlas.md#api-build_ggchord_cds) | [↗](homology.md) |
| [`bundle_ggchord_ribbons`](../../man/bundle_ggchord_ribbons.Rd) | 显式聚束密集比对连接 | [用法](ribbons.md#statistics-bundle) | 处理已有比对；不生成新的序列比对 | [↗](parameter-atlas.md#api-bundle_ggchord_ribbons) | [↗](chord.md) |
| [`clean_ggchord_data`](../../man/clean_ggchord_data.Rd) | 按明确策略清理输入并返回报告 | [用法](data.md#repairs) | 显式校验/修复/成员工作流；保留源表及报告 | [↗](parameter-atlas.md#api-clean_ggchord_data) | [↗](real-data.md) |
| [`coord_chord`](../../man/coord_chord.Rd) | 弦图坐标系 | [用法](chord-geometry.md#shape) | 弦图；整图旋转与逐序列几何分工 | [↗](parameter-atlas.md#api-coord_chord) | [↗](chord.md) |
| [`coord_circular`](../../man/coord_circular.Rd) | 单序列圆图坐标系 | [用法](circular-layout.md#coordinates) | 单条序列；bp/aa 布局；不接受 ribbon/link | [↗](parameter-atlas.md#api-coord_circular) | [↗](circular-start.md) |
| [`coord_collinear`](../../man/coord_collinear.Rd) | 直线型多序列共线性坐标系 | [用法](collinear-layout.md#parameters) | 共线性；bp/aa；拒绝圆图专属参数 | [↗](parameter-atlas.md#api-coord_collinear) | [↗](collinear.md) |
| [`coord_dotplot`](../../man/coord_dotplot.Rd) | HSP双轴投影 | [用法](dotplot-start.md#dotplot) | 独立原生轴上的identity HSP线段 | [↗](parameter-atlas.md#api-coord_dotplot) | [↗](dotplot-start.md) |
| [`coord_radial_synteny`](../../man/coord_radial_synteny.Rd) | 将线性序列投影到径向辐条 | [用法](radial-synteny.md#radial) | 比例或equal辐条、真实环形邻接HSP、共享Feature与Track框架 | [↗](parameter-atlas.md#api-coord_radial_synteny) | [↗](radial-synteny.md) |
| [`coord_sankey`](../../man/coord_sankey.Rd) | query 左列与 subject 右列镜像区间关系 | [用法](sankey.md#sankey) | 首批基础架构范围；详见本章支持矩阵 | [↗](parameter-atlas.md#api-coord_sankey) | [↗](multi-sequence.md) |
| [`deduplicate_ggchord_ribbons`](../../man/deduplicate_ggchord_ribbons.Rd) | 去除重复比对区间 | [用法](ribbons.md#transformations) | 处理已有比对；不生成新的序列比对 | [↗](parameter-atlas.md#api-deduplicate_ggchord_ribbons) | [↗](chord.md) |
| [`diagnose_ggchord`](../../man/diagnose_ggchord.Rd) | 检查指定尺寸的一次渲染而不改变原图 | [用法](diagnostics.md#target-size) | 单一目标尺寸和设备；提示不等于全面质量保证 | [↗](parameter-atlas.md#api-diagnose_ggchord) | [↗](chord.md) |
| [`element_blank`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-element_blank) | [↗](chord.md) |
| [`element_line`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-element_line) | [↗](chord.md) |
| [`element_rect`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-element_rect) | [↗](chord.md) |
| [`element_text`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-element_text) | [↗](chord.md) |
| [`expansion`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-expansion) | [↗](chord.md) |
| [`export_ggchord_fusion`](../../man/export_ggchord_fusion.Rd) | 导出实测组合前缘与来源 | [用法](fusion-routing.md#export) | 已注册基因组图、嵌套 Fusion 与已注册字符细节；显式面板/实例路由；不含 Dotplot | [↗](parameter-atlas.md#api-export_ggchord_fusion) | [↗](fusion.md) |
| [`export_ggchord_layout`](../../man/export_ggchord_layout.Rd) | 导出 ggchord 计算几何 | [用法](diagnostics.md#geometry) | 所有受支持图；按图层导出公共表 | [↗](parameter-atlas.md#api-export_ggchord_layout) | [↗](chord.md) |
| [`filter_ggchord_ribbons`](../../man/filter_ggchord_ribbons.Rd) | 绘图前筛选比对区间 | [用法](ribbons.md#input) | 处理已有比对；不生成新的序列比对 | [↗](parameter-atlas.md#api-filter_ggchord_ribbons) | [↗](chord.md) |
| [`filter_restriction_sites`](../../man/filter_restriction_sites.Rd) | 筛选要显示的酶切位点 | [用法](restriction-sites.md#parameters) | 显示集合选择；保留完整目标上的位点计数 | [↗](parameter-atlas.md#api-filter_restriction_sites) | [↗](plasmids.md) |
| [`find_common_features`](../../man/find_common_features.Rd) | 在 DNA 中查找常见生物学注释 | [用法](recognition.md#features) | DNA/参考注释识别；保留方法与来源 | [↗](parameter-atlas.md#api-find_common_features) | [↗](plasmids.md) |
| [`find_ggchord_loci`](../../man/find_ggchord_loci.Rd) | 围绕明确指定的注释查找基因组区域 | [用法](feature-geometry.md#locus-selection) | 基因组窗口；显式目标 ID 与真实链 | [↗](parameter-atlas.md#api-find_ggchord_loci) | [↗](views.md) |
| [`find_primer_bindings`](../../man/find_primer_bindings.Rd) | 查找测序引物结合位点 | [用法](primers.md#input) | 实际 DNA 与寡核苷酸；命中不等于实验验证 | [↗](parameter-atlas.md#api-find_primer_bindings) | [↗](plasmids.md) |
| [`find_restriction_sites`](../../man/find_restriction_sites.Rd) | 查找酶切模式匹配及切割坐标 | [用法](restriction-sites.md#input) | 实际 DNA；明确甲基化假设；完整位点结果 | [↗](parameter-atlas.md#api-find_restriction_sites) | [↗](plasmids.md) |
| [`focus_ggchord_data`](../../man/focus_ggchord_data.Rd) | 将 ggchord 数据聚焦到一个或多个区域 | [用法](views.md#focus) | 显式输入、窗口与来源的工作流 | [↗](parameter-atlas.md#api-focus_ggchord_data) | [↗](homology.md) |
| [`focus_ggchord_track`](../../man/focus_ggchord_track.Rd) | 聚焦定量区间或采样点而不跨越省略窗口 | [用法](tracks-statistics.md#track-focus) | 点/区间轨道；与显示窗口共享 loci | [↗](parameter-atlas.md#api-focus_ggchord_track) | [↗](tracks.md) |
| [`fusion_panel`](../../man/fusion_panel.Rd) | 注册指定方位的物理面板 | [用法](fusion-layout.md#advanced-placement) | 已注册基因组/字符面板；八方位、显式区域和嵌套路径 | [↗](parameter-atlas.md#api-fusion_panel) | [↗](fusion.md) |
| [`fusion_ribbon`](../../man/fusion_ribbon.Rd) | 连接完整来源区间 | [用法](fusion-routing.md#intervals) | 已注册基因组图、嵌套 Fusion 与已注册字符细节；显式面板/实例路由；不含 Dotplot | [↗](parameter-atlas.md#api-fusion_ribbon) | [↗](fusion.md) |
| [`geom_alignment_segment`](../../man/geom_alignment_segment.Rd) | 绘制HSP完整端点 | [用法](dotplot-start.md#dotplot) | 独立原生轴上的identity HSP线段 | [↗](parameter-atlas.md#api-geom_alignment_segment) | [↗](dotplot-start.md) |
| [`geom_feature`](../../man/geom_feature.Rd) | 绘制通用基因组注释 | [用法](circular-features.md#feature-shapes) | 兼容骨架坐标（见支持表）；通用注释可保留未知链及 aa 区间；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_feature) | [↗](chord.md) |
| [`geom_feature_label`](../../man/geom_feature_label.Rd) | 为通用基因组注释添加标签 | [用法](circular-labels.md#fixed-labels) | 兼容骨架坐标（见支持表）；直线图的 radial/tangent 为竖排/横排；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_feature_label) | [↗](chord.md) |
| [`geom_feature_label_repel`](../../man/geom_feature_label_repel.Rd) | 自动安排通用注释标签 | [用法](circular-labels.md#automatic-labels) | 兼容骨架坐标（见支持表）；直线图拒绝 external=FALSE；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_feature_label_repel) | [↗](chord.md) |
| [`geom_feature_plasmid`](../../man/geom_feature_plasmid.Rd) | 在圆形质粒上绘制 SnapGene 风格注释 | [用法](circular-features.md#feature-shapes) | 通用几何的组合预设；默认面向圆图；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_feature_plasmid) | [↗](plasmids.md) |
| [`geom_gene`](../../man/geom_gene.Rd) | 添加基因箭头图层 | [用法](feature-geometry.md#arrow-heads) | 兼容骨架坐标（见支持表）的 bp 数据；需已知正负链；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_gene) | [↗](chord.md) |
| [`geom_gene_label`](../../man/geom_gene_label.Rd) | 添加固定基因标签图层 | [用法](labels.md#fixed) | bp 基因标签；各坐标的方向含义不同；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_gene_label) | [↗](chord.md) |
| [`geom_gene_label_repel`](../../man/geom_gene_label_repel.Rd) | 添加自动布局的基因标签图层 | [用法](labels.md#comparison) | bp 基因标签；直线图统一使用行布局；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_gene_label_repel) | [↗](chord.md) |
| [`geom_link_line`](../../man/geom_link_line.Rd) | 绘制序列间的点连接 | [用法](multi-sequence-links.md#point-links) | 弦图/共线性；需 qpos/spos；圆图拒绝；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_link_line) | [↗](chord.md) |
| [`geom_link_ribbon`](../../man/geom_link_ribbon.Rd) | 添加比对区间连接图层 | [用法](multi-sequence-links.md#branches) | 弦图/共线性；邻接按显示顺序；圆图拒绝；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_link_ribbon) | [↗](chord.md) |
| [`geom_primer`](../../man/geom_primer.Rd) | 在圆形序列上绘制引物结合范围 | [用法](primers.md#primer-arcs) | bp；推荐单圆图；单序列弦图可组合；共线性拒绝 | [↗](parameter-atlas.md#api-geom_primer) | [↗](plasmids.md) |
| [`geom_primer_label_repel`](../../man/geom_primer_label_repel.Rd) | 排列引物标签和引导线 | [用法](primers.md#primer-labels) | bp；推荐单圆图；共线性拒绝 | [↗](parameter-atlas.md#api-geom_primer_label_repel) | [↗](plasmids.md) |
| [`geom_residue`](../../man/geom_residue.Rd) | 按实际物理空间绘制局部字符，明确选择色格或报错策略 | [用法](aligned-sequences.md#aligned-characters) | 已有DNA/RNA/AA对齐列；Cartesian局部单元与显式源映射 | [↗](parameter-atlas.md#api-geom_residue) | [↗](formats.md) |
| [`geom_restriction_site`](../../man/geom_restriction_site.Rd) | 绘制酶切位点、标签和确定性引导线 | [用法](restriction-sites.md#restriction-labels) | bp；推荐单圆图；共线性拒绝 | [↗](parameter-atlas.md#api-geom_restriction_site) | [↗](plasmids.md) |
| [`geom_seq`](../../man/geom_seq.Rd) | 添加序列骨架图层 | [用法](circular-layout.md#backbones) | 兼容骨架坐标（见支持表）；直线图拒绝半径/曲率/ring；圆图自管方向与开口；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_seq) | [↗](chord.md) |
| [`geom_seq_break`](../../man/geom_seq_break.Rd) | 标记线性序列的裁剪端点 | [用法](views.md#segments) | 仅共线性；需要聚焦/分段端点元数据；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_seq_break) | [↗](collinear.md) |
| [`geom_seq_center_label`](../../man/geom_seq_center_label.Rd) | 添加圆形序列中心标签 | [用法](circular-layout.md#center) | 单条序列；推荐圆图；共线性拒绝 | [↗](parameter-atlas.md#api-geom_seq_center_label) | [↗](chord.md) |
| [`geom_seq_label`](../../man/geom_seq_label.Rd) | 添加序列名称图层 | [用法](multi-sequence-labels.md#sequence-parameters) | 兼容骨架坐标（见支持表）；直线图拒绝显式径向半径；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_seq_label) | [↗](chord.md) |
| [`geom_seq_region`](../../man/geom_seq_region.Rd) | 沿序列骨架突出显示区域 | [用法](annotations.md#region-parameters) | 兼容骨架坐标（见支持表）；区间、ID、单位需合法；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_seq_region) | [↗](chord.md) |
| [`geom_seq_scale`](../../man/geom_seq_scale.Rd) | 在线性图中添加共享的序列位置比例尺 | [用法](collinear-layout.md#gap) | 仅共线性；原生 bp/aa 标尺 | [↗](parameter-atlas.md#api-geom_seq_scale) | [↗](collinear.md) |
| [`geom_track`](../../man/geom_track.Rd) | 源空间统计轨道 | [用法](tracks-statistics.md#tracks) | 首批基础架构范围；详见本章支持矩阵；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_track) | [↗](multi-sequence.md) |
| [`geom_track_bar`](../../man/geom_track.Rd) | 按源位置统计观测数 | [用法](tracks-statistics.md#tracks) | 兼容骨架坐标中的源位置观测计数；基线仅 0；分箱是另行声明的 Stat | [↗](parameter-atlas.md#api-geom_track_bar) | [↗](multi-sequence.md) |
| [`geom_track_col`](../../man/geom_track.Rd) | 显示给定数值高度 | [用法](tracks-statistics.md#tracks) | 首批基础架构范围；详见本章支持矩阵；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_track_col) | [↗](multi-sequence.md) |
| [`geom_track_line`](../../man/geom_track_tile.Rd) | 沿序列骨架绘制定量轨道 | [用法](tracks.md#multiple) | 兼容骨架坐标（见支持表）；每组位置递增；closed 需完整闭圆；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_track_line) | [↗](tracks.md) |
| [`geom_track_tile`](../../man/geom_track_tile.Rd) | 沿序列骨架绘制定量轨道 | [用法](tracks.md#tile) | 兼容骨架坐标（见支持表）；显式轨道；重叠按输入顺序绘制；Radial见新增视图支持矩阵 | [↗](parameter-atlas.md#api-geom_track_tile) | [↗](tracks.md) |
| [`get_chord_layout`](../../man/get_chord_layout.Rd) | 获取 ggchord 图的完整布局 | [用法](diagnostics.md#geometry) | 所有受支持图；完整文档化布局 | [↗](parameter-atlas.md#api-get_chord_layout) | [↗](chord.md) |
| [`ggchord`](../../man/ggchord.Rd) | 用数据表和独立图层创建序列图 | [用法](modes.md#composition) | 兼容坐标（见支持表）；同一图的 bp/aa 单位必须一致 | [↗](parameter-atlas.md#api-ggchord) | [↗](real-data.md) |
| [`ggchord_fusion`](../../man/ggchord_fusion.Rd) | 组合独立基因组面板 | [用法](fusion-routing.md#intervals) | 已注册基因组/字符面板；八方位、显式区域和嵌套路径 | [↗](parameter-atlas.md#api-ggchord_fusion) | [↗](fusion.md) |
| [`ggsave`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-ggsave) | [↗](chord.md) |
| [`ggtitle`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-ggtitle) | [↗](chord.md) |
| [`guide_ggchord_colorbar`](../../man/guide_ggchord_legend.Rd) | ggchord 角色 aesthetic 图例；同义拼写 guide_ggchord_colourbar | [用法](scale-families.md#guides) | 受支持图的相应主题/图例；按角色及设备解释 | [↗](parameter-atlas.md#api-guide_ggchord_colorbar) | [↗](chord.md) |
| [`guide_ggchord_colourbar`](../../man/guide_ggchord_legend.Rd) | ggchord 角色 aesthetic 图例 | [用法](scale-families.md#guides) | 受支持图的相应主题/图例；按角色及设备解释 | [↗](parameter-atlas.md#api-guide_ggchord_colourbar) | [↗](chord.md) |
| [`guide_ggchord_legend`](../../man/guide_ggchord_legend.Rd) | ggchord 角色 aesthetic 图例 | [用法](scale-families.md#guides) | 受支持图的相应主题/图例；按角色及设备解释 | [↗](parameter-atlas.md#api-guide_ggchord_legend) | [↗](chord.md) |
| [`guides`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-guides) | [↗](chord.md) |
| [`labs`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-labs) | [↗](chord.md) |
| [`last_plot`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-last_plot) | [↗](chord.md) |
| [`link_ggchord_clusters`](../../man/link_ggchord_clusters.Rd) | 将显式同源分组成员转换为注释连接 | [用法](homology.md#links) | 消费明确成员表；同名不能证明同源 | [↗](parameter-atlas.md#api-link_ggchord_clusters) | [↗](homology.md) |
| [`magnify_ggchord`](../../man/magnify_ggchord.Rd) | 连接真实源区间与显式独立构建的局部细节图 | [用法](fusion-details.md#magnify) | 已注册基因组图、嵌套 Fusion 与已注册字符细节；显式面板/实例路由；不含 Dotplot | [↗](parameter-atlas.md#api-magnify_ggchord) | [↗](fusion.md) |
| [`map_alignment_position`](../../man/map_alignment_position.Rd) | 在选定行中精确查找对齐列与源位置 | [用法](aligned-sequences.md#aligned-characters) | 已有DNA/RNA/AA对齐列；Cartesian局部单元与显式源映射 | [↗](parameter-atlas.md#api-map_alignment_position) | [↗](formats.md) |
| [`map_ggchord_position`](../../man/map_ggchord_position.Rd) | 映射源位置与显示坐标 | [用法](radial-synteny.md#export) | 原生或源坐标；Dotplot须选轴，Radial反查须给x/y | [↗](parameter-atlas.md#api-map_ggchord_position) | [↗](dotplot-start.md) |
| [`margin`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-margin) | [↗](chord.md) |
| [`merge_ggchord_ribbons`](../../man/merge_ggchord_ribbons.Rd) | 合并同序列对的相邻或重叠比对块 | [用法](ribbons.md#transformations) | 处理已有比对；不生成新的序列比对 | [↗](parameter-atlas.md#api-merge_ggchord_ribbons) | [↗](chord.md) |
| [`optimize_ggchord_layout`](../../man/optimize_ggchord_layout.Rd) | 优化序列顺序和方向以改善连接可读性 | [用法](diagnostics.md#target-size) | 明确请求的顺序/方向启发式；不是生物学推断 | [↗](parameter-atlas.md#api-optimize_ggchord_layout) | [↗](chord.md) |
| [`position_feature_stack`](../../man/position_feature_stack.Rd) | 将重叠基因或注释堆叠到局部轨道 | [用法](circular-features.md#stacking) | 先应用基础 Position，再分配重叠轨道 | [↗](parameter-atlas.md#api-position_feature_stack) | [↗](chord.md) |
| [`position_identity`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-position_identity) | [↗](chord.md) |
| [`position_plasmid`](../../man/position_plasmid.Rd) | 将两条链的注释放在同一带上 | [用法](circular-features.md#stacking) | 共享注释带的局部预设；默认面向圆图 | [↗](parameter-atlas.md#api-position_plasmid) | [↗](circular-features.md) |
| [`position_strand`](../../man/position_strand.Rd) | 按链将注释放在两侧 | [用法](parameters.md#strand) | 合法特征/基因的带符号局部放置 | [↗](parameter-atlas.md#api-position_strand) | [↗](parameters.md) |
| [`position_track`](../../man/position_track.Rd) | 在显式序列局部带中放置定量图层 | [用法](tracks-statistics.md#independent-bands) | 兼容骨架坐标（见支持表）；同 ID 的声明和值域需一致 | [↗](parameter-atlas.md#api-position_track) | [↗](tracks.md) |
| [`PositionIdentity`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-PositionIdentity) | [↗](chord.md) |
| [`prepare_ggchord_data`](../../man/prepare_ggchord_data.Rd) | 准备命名输入表并保留明确审计记录 | [用法](data.md#validation) | 显式校验/修复/成员工作流；保留源表及报告 | [↗](parameter-atlas.md#api-prepare_ggchord_data) | [↗](real-data.md) |
| [`primer_catalog`](../../man/primer_catalog.Rd) | 查看随包测序引物目录 | [用法](recognition.md#catalog) | 目录元数据；不是目标上的结合搜索 | [↗](parameter-atlas.md#api-primer_catalog) | [↗](plasmids.md) |
| [`project_ggchord_features`](../../man/project_ggchord_features.Rd) | 通过剪接 CDS 映射投影蛋白结构域或基因组注释 | [用法](projection.md#cds) | 显式 CDS/相位/单位的编码坐标工作流 | [↗](parameter-atlas.md#api-project_ggchord_features) | [↗](homology.md) |
| [`project_ggchord_links`](../../man/project_ggchord_links.Rd) | 将明确无缺口的蛋白连接投影为基因组编码片段 | [用法](projection.md#links) | 显式 CDS/相位/单位的编码坐标工作流 | [↗](parameter-atlas.md#api-project_ggchord_links) | [↗](homology.md) |
| [`read_aligned_fasta`](../../man/read_aligned_fasta.Rd) | 读取已有DNA/RNA/AA对齐字符，不计算新比对 | [用法](aligned-sequences.md#aligned-characters) | 已有DNA/RNA/AA对齐列；Cartesian局部单元与显式源映射 | [↗](parameter-atlas.md#api-read_aligned_fasta) | [↗](formats.md) |
| [`read_bed`](../../man/read_bed.Rd) | 将 BED 区间读为 ggchord 注释 | [用法](import-intervals.md#zero-based) | 本地文件到标准表；按源格式转换坐标和声明单位 | [↗](parameter-atlas.md#api-read_bed) | [↗](formats.md) |
| [`read_bedgraph`](../../man/read_bedgraph.Rd) | 导入连续区间信号，显式保留坐标及单位边界 | [用法](import-signals.md#signals) | 定量轨道输入；BigWig 为支持平台上的可选依赖 | [↗](parameter-atlas.md#api-read_bedgraph) | [↗](formats.md) |
| [`read_bedpe`](../../man/read_bedpe.Rd) | 读取双端区间，保留来源且不推断缺失比对证据 | [用法](import-intervals.md#bedpe) | 标准q/s区间、配对锚点方向与来源 | [↗](parameter-atlas.md#api-read_bedpe) | [↗](formats.md) |
| [`read_bigwig`](../../man/read_bigwig.Rd) | 导入连续区间信号，显式保留坐标及单位边界 | [用法](import-signals.md#signals) | 定量轨道输入；BigWig 为支持平台上的可选依赖 | [↗](parameter-atlas.md#api-read_bigwig) | [↗](formats.md) |
| [`read_blast`](../../man/read_blast.Rd) | 将一个或多个 BLAST 表格文件读为连接数据 | [用法](alignments.md#standard) | 本地文件到标准表；按源格式转换坐标和声明单位 | [↗](parameter-atlas.md#api-read_blast) | [↗](real-data.md) |
| [`read_chain`](../../man/read_chain.Rd) | 读取已有基因组对齐块及真实源坐标 | [用法](import-relations.md#genome-alignments) | 可选 Biopython；MAF 字符单元与 chain 区间 | [↗](parameter-atlas.md#api-read_chain) | [↗](formats.md) |
| [`read_delta`](../../man/read_delta.Rd) | 通过可选原生解码器读取已有 delta 比对 | [用法](import-relations.md#delta) | 可选 show-coords；bp 坐标与独立声明的比对统计单位 | [↗](parameter-atlas.md#api-read_delta) | [↗](formats.md) |
| [`read_fasta_lengths`](../../man/read_fasta_lengths.Rd) | 将一个或多个 FASTA 文件读为序列数据 | [用法](aligned-sequences.md#characters) | 本地文件到标准表；按源格式转换坐标和声明单位 | [↗](parameter-atlas.md#api-read_fasta_lengths) | [↗](real-data.md) |
| [`read_feature_relations`](../../man/read_feature_relations.Rd) | 读取预计算 feature 对及成员关系 | [用法](import-relations.md#feature-relations) | 精确 feature ID；显式连续锚点；不推断完整长度 | [↗](parameter-atlas.md#api-read_feature_relations) | [↗](formats.md) |
| [`read_genbank`](../../man/read_genbank.Rd) | 使用可选成熟解析器导入本地 GenBank 或 GBFF | [用法](import-sequences.md#genbank) | 本地文件到标准表；按源格式转换坐标和声明单位 | [↗](parameter-atlas.md#api-read_genbank) | [↗](formats.md) |
| [`read_gff3`](../../man/read_gff3.Rd) | 将一个或多个 GFF3 文件读为注释数据 | [用法](real-data.md#read) | 本地文件到标准表；按源格式转换坐标和声明单位 | [↗](parameter-atlas.md#api-read_gff3) | [↗](real-data.md) |
| [`read_maf`](../../man/read_maf.Rd) | 读取已有基因组对齐块及真实源坐标 | [用法](import-relations.md#genome-alignments) | 可选 Biopython；MAF 字符单元与 chain 区间 | [↗](parameter-atlas.md#api-read_maf) | [↗](formats.md) |
| [`read_paf`](../../man/read_paf.Rd) | 将 PAF 比对读为 ggchord 连接 | [用法](import-intervals.md#zero-based) | 本地文件到标准表；按源格式转换坐标和声明单位 | [↗](parameter-atlas.md#api-read_paf) | [↗](formats.md) |
| [`read_signature_matches`](../../man/read_signature_matches.Rd) | 读取已有蛋白 signature 位置及源元数据 | [用法](import-relations.md#protein-signatures) | AA feature；源原生分数；JSON 需可选 jsonlite | [↗](parameter-atlas.md#api-read_signature_matches) | [↗](formats.md) |
| [`read_vcf`](../../man/read_vcf.Rd) | 导入变异记录及显式解码的结构断点 | [用法](import-relations.md#variants) | VCF；BCF 需显式 pysam Python；结构解码为可选依赖 | [↗](parameter-atlas.md#api-read_vcf) | [↗](formats.md) |
| [`register_ggchord_alignment`](../../man/register_ggchord_alignment.Rd) | 注册精确字符/源位置用于独立Fusion细节 | [用法](fusion-details.md#character-detail) | 显式源/对齐列映射；独立Cartesian细节；已观测完整codon | [↗](parameter-atlas.md#api-register_ggchord_alignment) | [↗](fusion.md) |
| [`rel`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-rel) | [↗](chord.md) |
| [`scale_feature_fill_manual`](../../man/scale_gene_fill_manual.Rd) | 基因、注释及区域填色标度 | [用法](scale-families.md#features) | 兼容骨架坐标（见支持表）；匹配通用注释角色或预设 | [↗](parameter-atlas.md#api-scale_feature_fill_manual) | [↗](chord.md) |
| [`scale_feature_fill_plasmid`](../../man/scale_feature_fill_plasmid.Rd) | 质粒注释预设标度 | [用法](scale-families.md#features) | 兼容骨架坐标（见支持表）；匹配通用注释角色或预设 | [↗](parameter-atlas.md#api-scale_feature_fill_plasmid) | [↗](chord.md) |
| [`scale_feature_shape_manual`](../../man/scale_feature_shape_manual.Rd) | 注释形状标度 | [用法](scale-families.md#features) | 兼容骨架坐标（见支持表）；匹配通用注释角色或预设 | [↗](parameter-atlas.md#api-scale_feature_shape_manual) | [↗](chord.md) |
| [`scale_feature_shape_plasmid`](../../man/scale_feature_fill_plasmid.Rd) | 质粒注释预设标度 | [用法](scale-families.md#features) | 兼容骨架坐标（见支持表）；匹配通用注释角色或预设 | [↗](parameter-atlas.md#api-scale_feature_shape_plasmid) | [↗](chord.md) |
| [`scale_gene_fill_manual`](../../man/scale_gene_fill_manual.Rd) | 基因、注释及区域填色标度 | [用法](scale-families.md#rings) | 兼容骨架坐标（见支持表）的 bp 基因图层 | [↗](parameter-atlas.md#api-scale_gene_fill_manual) | [↗](chord.md) |
| [`scale_link_alpha_continuous`](../../man/scale_link.Rd) | 点连接标度 | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_alpha_continuous) | [↗](chord.md) |
| [`scale_link_alpha_manual`](../../man/scale_link.Rd) | 点连接标度 | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_alpha_manual) | [↗](chord.md) |
| [`scale_link_color_gradientn`](../../man/scale_link.Rd) | 点连接标度；同义拼写 scale_link_colour_gradientn | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_color_gradientn) | [↗](chord.md) |
| [`scale_link_color_manual`](../../man/scale_link.Rd) | 点连接标度；同义拼写 scale_link_colour_manual | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_color_manual) | [↗](chord.md) |
| [`scale_link_colour_gradientn`](../../man/scale_link.Rd) | 点连接标度 | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_colour_gradientn) | [↗](chord.md) |
| [`scale_link_colour_manual`](../../man/scale_link.Rd) | 点连接标度 | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_colour_manual) | [↗](chord.md) |
| [`scale_link_linetype_manual`](../../man/scale_link.Rd) | 点连接标度 | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_linetype_manual) | [↗](chord.md) |
| [`scale_link_linewidth_continuous`](../../man/scale_link.Rd) | 点连接标度 | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_linewidth_continuous) | [↗](chord.md) |
| [`scale_link_linewidth_manual`](../../man/scale_link.Rd) | 点连接标度 | [用法](scale-families.md#point-links) | 弦图/共线性；匹配点连接角色 | [↗](parameter-atlas.md#api-scale_link_linewidth_manual) | [↗](chord.md) |
| [`scale_region_fill_manual`](../../man/scale_gene_fill_manual.Rd) | 基因、注释及区域填色标度 | [用法](scale-families.md#rings) | 兼容骨架坐标（见支持表）；匹配 region_fill | [↗](parameter-atlas.md#api-scale_region_fill_manual) | [↗](chord.md) |
| [`scale_ribbon_alpha_continuous`](../../man/scale_ribbon_alpha_continuous.Rd) | 区间连接透明度、轮廓及线型标度 | [用法](scale-families.md#categories) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_alpha_continuous) | [↗](chord.md) |
| [`scale_ribbon_alpha_manual`](../../man/scale_ribbon_alpha_continuous.Rd) | 区间连接透明度、轮廓及线型标度 | [用法](scale-families.md#categories) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_alpha_manual) | [↗](chord.md) |
| [`scale_ribbon_color_manual`](../../man/scale_ribbon_alpha_continuous.Rd) | 区间连接透明度、轮廓及线型标度；同义拼写 scale_ribbon_colour_manual | [用法](scale-families.md#categories) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_color_manual) | [↗](chord.md) |
| [`scale_ribbon_colour_manual`](../../man/scale_ribbon_alpha_continuous.Rd) | 区间连接透明度、轮廓及线型标度 | [用法](scale-families.md#categories) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_colour_manual) | [↗](chord.md) |
| [`scale_ribbon_fill_gradientn`](../../man/scale_ribbon_fill_stepsn.Rd) | 区间连接填色标度 | [用法](scale-families.md#numeric) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_fill_gradientn) | [↗](chord.md) |
| [`scale_ribbon_fill_identity`](../../man/scale_ribbon_fill_stepsn.Rd) | 区间连接填色标度 | [用法](scale-families.md#identity) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_fill_identity) | [↗](chord.md) |
| [`scale_ribbon_fill_manual`](../../man/scale_ribbon_fill_stepsn.Rd) | 区间连接填色标度 | [用法](scale-families.md#categories) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_fill_manual) | [↗](chord.md) |
| [`scale_ribbon_fill_stepsn`](../../man/scale_ribbon_fill_stepsn.Rd) | 区间连接填色标度 | [用法](scale-families.md#numeric) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_fill_stepsn) | [↗](chord.md) |
| [`scale_ribbon_linetype_manual`](../../man/scale_ribbon_alpha_continuous.Rd) | 区间连接透明度、轮廓及线型标度 | [用法](scale-families.md#categories) | 弦图/共线性；匹配连接带角色 | [↗](parameter-atlas.md#api-scale_ribbon_linetype_manual) | [↗](chord.md) |
| [`scale_seq_color_manual`](../../man/scale_seq_colour_manual.Rd) | 序列颜色标度；同义拼写 scale_seq_colour_manual | [用法](scale-families.md#categories) | 兼容骨架坐标（见支持表）；匹配序列颜色角色 | [↗](parameter-atlas.md#api-scale_seq_color_manual) | [↗](chord.md) |
| [`scale_seq_colour_manual`](../../man/scale_seq_colour_manual.Rd) | 序列颜色标度 | [用法](scale-families.md#categories) | 兼容骨架坐标（见支持表）；匹配序列颜色角色 | [↗](parameter-atlas.md#api-scale_seq_colour_manual) | [↗](chord.md) |
| [`scale_seq_position_continuous`](../../man/scale_seq_position_continuous.Rd) | 基因组序列位置标度 | [用法](views.md#focus) | 兼容骨架坐标（见支持表）的位置轴；与直线标尺不同 | [↗](parameter-atlas.md#api-scale_seq_position_continuous) | [↗](chord.md) |
| [`scale_seq_ring_manual`](../../man/scale_seq_ring_manual.Rd) | 显式序列环标度 | [用法](scale-families.md#rings) | 弦图/圆图；显式正半径；共线性拒绝 | [↗](parameter-atlas.md#api-scale_seq_ring_manual) | [↗](chord.md) |
| [`scale_track_colour_gradientn`](../../man/scale_track_fill_gradientn.Rd) | 设置轨道填色或颜色映射 | [用法](scale-families.md#tracks) | 兼容骨架坐标（见支持表）；ID 必须匹配 position_track | [↗](parameter-atlas.md#api-scale_track_colour_gradientn) | [↗](tracks.md) |
| [`scale_track_colour_manual`](../../man/scale_track_fill_gradientn.Rd) | 设置轨道填色或颜色映射 | [用法](scale-families.md#tracks) | 兼容骨架坐标（见支持表）；ID 必须匹配 position_track | [↗](parameter-atlas.md#api-scale_track_colour_manual) | [↗](tracks.md) |
| [`scale_track_fill_gradientn`](../../man/scale_track_fill_gradientn.Rd) | 设置轨道填色或颜色映射 | [用法](scale-families.md#tracks) | 兼容骨架坐标（见支持表）；ID 必须匹配 position_track | [↗](parameter-atlas.md#api-scale_track_fill_gradientn) | [↗](tracks.md) |
| [`scale_track_fill_manual`](../../man/scale_track_fill_gradientn.Rd) | 设置轨道填色或颜色映射 | [用法](scale-families.md#tracks) | 兼容骨架坐标（见支持表）；ID 必须匹配 position_track | [↗](parameter-atlas.md#api-scale_track_fill_manual) | [↗](tracks.md) |
| [`scale_track_value_continuous`](../../man/scale_track_value_continuous.Rd) | 设置命名轨道的数值高度标度 | [用法](scale-families.md#tracks) | 兼容骨架坐标（见支持表）；ID 必须匹配 position_track | [↗](parameter-atlas.md#api-scale_track_value_continuous) | [↗](tracks.md) |
| [`segment_ggchord_data`](../../man/segment_ggchord_data.Rd) | 为分段线性图准备独立窗口 | [用法](views.md#segments) | 显式输入、窗口与来源的工作流 | [↗](parameter-atlas.md#api-segment_ggchord_data) | [↗](views.md) |
| [`stat_ribbon_bundle`](../../man/stat_ribbon_bundle.Rd) | 以统计图层聚束区间连接 | [用法](ribbons.md#statistics-bundle) | 弦图/共线性；显式聚束 | [↗](parameter-atlas.md#api-stat_ribbon_bundle) | [↗](chord.md) |
| [`stat_ribbon_density`](../../man/stat_ribbon_density.Rd) | 以统计图层计算局部连接密度 | [用法](ribbons.md#statistics-density) | 弦图/共线性；一行观测仍是一条带 | [↗](parameter-atlas.md#api-stat_ribbon_density) | [↗](chord.md) |
| [`theme`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-theme) | [↗](chord.md) |
| [`theme_ggchord`](../../man/theme_ggchord.Rd) | ggchord 主题 | [用法](scale-families.md#themes) | 受支持图的相应主题/图例；按角色及设备解释 | [↗](parameter-atlas.md#api-theme_ggchord) | [↗](chord.md) |
| [`theme_ggchord_dark`](../../man/theme_ggchord.Rd) | ggchord 主题 | [用法](scale-families.md#themes) | 受支持图的相应主题/图例；按角色及设备解释 | [↗](parameter-atlas.md#api-theme_ggchord_dark) | [↗](chord.md) |
| [`theme_ggchord_minimal`](../../man/theme_ggchord.Rd) | ggchord 主题 | [用法](scale-families.md#themes) | 受支持图的相应主题/图例；按角色及设备解释 | [↗](parameter-atlas.md#api-theme_ggchord_minimal) | [↗](chord.md) |
| [`theme_ggchord_plasmid`](../../man/theme_ggchord.Rd) | ggchord 主题 | [用法](scale-families.md#themes) | 受支持图的相应主题/图例；按角色及设备解释 | [↗](parameter-atlas.md#api-theme_ggchord_plasmid) | [↗](chord.md) |
| [`theme_ggchord_publication`](../../man/theme_ggchord.Rd) | ggchord 主题 | [用法](scale-families.md#themes) | 受支持图的相应主题/图例；按角色及设备解释 | [↗](parameter-atlas.md#api-theme_ggchord_publication) | [↗](chord.md) |
| [`translate_ggchord_cds`](../../man/translate_ggchord_cds.Rd) | 通过可选Biostrings翻译真实CDS并保留codon来源 | [用法](aligned-sequences.md#cds-translation) | 显式源/对齐列映射；独立Cartesian细节；已观测完整codon | [↗](parameter-atlas.md#api-translate_ggchord_cds) | [↗](formats.md) |
| [`unit`](../../man/grid-helpers.Rd) | 物理单位与方向箭头 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-unit) | [↗](chord.md) |
| [`validate_ggchord_data`](../../man/validate_ggchord_data.Rd) | 绘图前校验 ggchord 输入 | [用法](data.md#validation) | 显式校验/修复/成员工作流；保留源表及报告 | [↗](parameter-atlas.md#api-validate_ggchord_data) | [↗](real-data.md) |
| [`validate_ggchord_hierarchy`](../../man/validate_ggchord_hierarchy.Rd) | 校验显式 genome 与 contig 关系 | [用法](data.md#hierarchy) | 显式校验/修复/成员工作流；保留源表及报告 | [↗](parameter-atlas.md#api-validate_ggchord_hierarchy) | [↗](real-data.md) |
| [`view_ggchord`](../../man/view_ggchord.Rd) | 按预期导出尺寸预览图形 | [用法](output.md#size) | 所有受支持图；目标尺寸临时预览 | [↗](parameter-atlas.md#api-view_ggchord) | [↗](chord.md) |
| [`waiver`](../../man/ggplot2-helpers.Rd) | 适用于 ggchord 的部分 ggplot2 辅助函数 | [上游](parameter-atlas.md#upstream) | 上游对象/接口；遵循 ggplot2/grid 的原始契约 | [↗](parameter-atlas.md#api-waiver) | [↗](chord.md) |

<a id="reproduce"></a>

## 怎样复现 {#reproduce}

按主讲页从准备代码开始运行；字段含义在[基础表](data.md)，实际格式转换在[导入](formats.md)。所有主讲页的共享片段在独立环境执行，重要对照图来自真实 ggchord 渲染。已标明的简化几何示意不作为生物学证据。原始输入、完整源码和复现说明随对应页面保留；安装方法见[安装页](installation.md)。

<details class="maintenance-record"><summary>维护覆盖记录</summary>

L2 表示最小示例与输入/输出说明，L3 另含关键参数对照；均不保证所有参数组合。

| API | 等级 |
| --- | --- |
| `adapt_ggchord` | L2 |
| `aes` | L1 |
| `after_scale` | L1 |
| `after_stat` | L1 |
| `annotate` | L1 |
| `arrow` | L1 |
| `as_ggchord_alignment` | L2 |
| `as_ggchord_feature` | L2 |
| `as_ggchord_link` | L2 |
| `as_ggchord_seq` | L2 |
| `build_ggchord_cds` | L2 |
| `bundle_ggchord_ribbons` | L3 |
| `clean_ggchord_data` | L2 |
| `coord_chord` | L3 |
| `coord_circular` | L3 |
| `coord_collinear` | L3 |
| `coord_dotplot` | L2 |
| `coord_radial_synteny` | L2 |
| `coord_sankey` | L2 |
| `deduplicate_ggchord_ribbons` | L2 |
| `diagnose_ggchord` | L2 |
| `element_blank` | L1 |
| `element_line` | L1 |
| `element_rect` | L1 |
| `element_text` | L1 |
| `expansion` | L1 |
| `export_ggchord_fusion` | L2 |
| `export_ggchord_layout` | L2 |
| `filter_ggchord_ribbons` | L2 |
| `filter_restriction_sites` | L2 |
| `find_common_features` | L2 |
| `find_ggchord_loci` | L2 |
| `find_primer_bindings` | L2 |
| `find_restriction_sites` | L2 |
| `focus_ggchord_data` | L2 |
| `focus_ggchord_track` | L2 |
| `fusion_panel` | L2 |
| `fusion_ribbon` | L2 |
| `geom_alignment_segment` | L2 |
| `geom_feature` | L3 |
| `geom_feature_label` | L3 |
| `geom_feature_label_repel` | L3 |
| `geom_feature_plasmid` | L3 |
| `geom_gene` | L3 |
| `geom_gene_label` | L3 |
| `geom_gene_label_repel` | L3 |
| `geom_link_line` | L3 |
| `geom_link_ribbon` | L3 |
| `geom_primer` | L3 |
| `geom_primer_label_repel` | L3 |
| `geom_residue` | L2 |
| `geom_restriction_site` | L3 |
| `geom_seq` | L3 |
| `geom_seq_break` | L2 |
| `geom_seq_center_label` | L3 |
| `geom_seq_label` | L3 |
| `geom_seq_region` | L3 |
| `geom_seq_scale` | L3 |
| `geom_track` | L2 |
| `geom_track_bar` | L3 |
| `geom_track_col` | L2 |
| `geom_track_line` | L3 |
| `geom_track_tile` | L3 |
| `get_chord_layout` | L2 |
| `ggchord` | L2 |
| `ggchord_fusion` | L2 |
| `ggsave` | L1 |
| `ggtitle` | L1 |
| `guide_ggchord_colorbar` | L2 |
| `guide_ggchord_colourbar` | L2 |
| `guide_ggchord_legend` | L2 |
| `guides` | L1 |
| `labs` | L1 |
| `last_plot` | L1 |
| `link_ggchord_clusters` | L2 |
| `magnify_ggchord` | L2 |
| `map_alignment_position` | L2 |
| `map_ggchord_position` | L2 |
| `margin` | L1 |
| `merge_ggchord_ribbons` | L2 |
| `optimize_ggchord_layout` | L2 |
| `position_feature_stack` | L3 |
| `position_identity` | L1 |
| `position_plasmid` | L2 |
| `position_strand` | L2 |
| `position_track` | L3 |
| `PositionIdentity` | L1 |
| `prepare_ggchord_data` | L2 |
| `primer_catalog` | L2 |
| `project_ggchord_features` | L2 |
| `project_ggchord_links` | L2 |
| `read_aligned_fasta` | L2 |
| `read_bed` | L2 |
| `read_bedgraph` | L2 |
| `read_bedpe` | L2 |
| `read_bigwig` | L2 |
| `read_blast` | L2 |
| `read_chain` | L2 |
| `read_delta` | L2 |
| `read_fasta_lengths` | L2 |
| `read_feature_relations` | L2 |
| `read_genbank` | L2 |
| `read_gff3` | L2 |
| `read_maf` | L2 |
| `read_paf` | L2 |
| `read_signature_matches` | L2 |
| `read_vcf` | L2 |
| `register_ggchord_alignment` | L2 |
| `rel` | L1 |
| `scale_feature_fill_manual` | L2 |
| `scale_feature_fill_plasmid` | L2 |
| `scale_feature_shape_manual` | L2 |
| `scale_feature_shape_plasmid` | L2 |
| `scale_gene_fill_manual` | L2 |
| `scale_link_alpha_continuous` | L2 |
| `scale_link_alpha_manual` | L2 |
| `scale_link_color_gradientn` | L2 |
| `scale_link_color_manual` | L2 |
| `scale_link_colour_gradientn` | L2 |
| `scale_link_colour_manual` | L2 |
| `scale_link_linetype_manual` | L2 |
| `scale_link_linewidth_continuous` | L2 |
| `scale_link_linewidth_manual` | L2 |
| `scale_region_fill_manual` | L2 |
| `scale_ribbon_alpha_continuous` | L2 |
| `scale_ribbon_alpha_manual` | L2 |
| `scale_ribbon_color_manual` | L2 |
| `scale_ribbon_colour_manual` | L2 |
| `scale_ribbon_fill_gradientn` | L2 |
| `scale_ribbon_fill_identity` | L2 |
| `scale_ribbon_fill_manual` | L2 |
| `scale_ribbon_fill_stepsn` | L2 |
| `scale_ribbon_linetype_manual` | L2 |
| `scale_seq_color_manual` | L2 |
| `scale_seq_colour_manual` | L2 |
| `scale_seq_position_continuous` | L2 |
| `scale_seq_ring_manual` | L2 |
| `scale_track_colour_gradientn` | L2 |
| `scale_track_colour_manual` | L2 |
| `scale_track_fill_gradientn` | L2 |
| `scale_track_fill_manual` | L2 |
| `scale_track_value_continuous` | L2 |
| `segment_ggchord_data` | L2 |
| `stat_ribbon_bundle` | L3 |
| `stat_ribbon_density` | L3 |
| `theme` | L1 |
| `theme_ggchord` | L2 |
| `theme_ggchord_dark` | L2 |
| `theme_ggchord_minimal` | L2 |
| `theme_ggchord_plasmid` | L2 |
| `theme_ggchord_publication` | L2 |
| `translate_ggchord_cds` | L2 |
| `unit` | L1 |
| `validate_ggchord_data` | L2 |
| `validate_ggchord_hierarchy` | L2 |
| `view_ggchord` | L2 |
| `waiver` | L1 |

</details>
