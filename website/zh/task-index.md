按想做的事情选择课程。A 为弦图，B 为共线性图，C 为单序列圆图；带条件的条目先查看[兼容边界](compatibility.html)。基础内容完成后，再按需选择进阶内容。可用浏览器查找任务关键词，站内搜索也会索引本页。

## 功能与参数索引 {#tasks}

| 任务 | 函数 | 关键参数 | 模式 | 课程与层级 |
| --- | --- | --- | --- | --- |
| 改变骨架样式 | [`geom_seq()`](reference/geom_seq.html) | `seq_style, seq_backbone_gap, seq_backbone_width` | A / B / C | [基础](plasmid-tuning.html#style) |
| 指定半径与曲率 | [`geom_seq()`](reference/geom_seq.html) | `seq_radius, seq_curvature` | A | [进阶](chord-geometry.html) |
| 改变序列顺序和方向 | [`geom_seq()`](reference/geom_seq.html) | `seq_order, seq_orientation` | A / B | [进阶](chord-geometry.html) |
| 改变连接路径与避障 | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `link_ctrl_point, link_gap, link_avoid` | A / B | [进阶](chord-routing.html) |
| 将色带改为有向点连接 | [`geom_link_line()`](reference/geom_link_line.html) | `link_type, qpos/spos, arrow` | A / B | [进阶](chord-branches.html#point-links) |
| 调整基因箭头形状 | [`geom_gene()`](reference/geom_gene.html) | `arrow_head_style, short_feature` | A / B / C (bp) | [进阶](feature-geometry.html#arrow-heads) |
| 强调一个局部区间 | [`geom_seq_region()`](reference/geom_seq_region.html) | `region_side, region_width, region_offset` | A / B / C | [进阶](chord-routing.html) |
| 根据相似度配色 | [`scale_ribbon_fill_stepsn()`](reference/scale_ribbon_fill_stepsn.html) | `ribbon_fill` | A / B | [基础](adjustments.html) |
| 分类配色与图例 | [`scale_seq_colour_manual()`](reference/scale_seq_colour_manual.html) | `values, guide` | A / B / C | [进阶](styling.html) |
| 筛选、去重、合并比对 | [`filter_ggchord_ribbons()`](reference/filter_ggchord_ribbons.html) | `min_pident, min_length` | A / B | [进阶](ribbons.html) |
| 汇总连接密度 | [`stat_ribbon_density()`](reference/stat_ribbon_density.html) | `bins, weight` | A / B | [进阶](ribbon-statistics.html#density) |
| 显式聚束连接 | [`stat_ribbon_bundle()`](reference/stat_ribbon_bundle.html) | `min_bundle, group_by` | A / B | [进阶](ribbon-statistics.html#bundle) |
| 左中右对齐序列 | [`coord_collinear()`](reference/coord_collinear.html) | `align` | B | [基础](collinear-layout.html#alignment) |
| 以共同位点锚定 | [`coord_collinear()`](reference/coord_collinear.html) | `anchor` | B | [进阶](collinear-layout.html#anchor) |
| 自动选择显示方向 | [`coord_collinear()`](reference/coord_collinear.html) | `orientation, min_support` | B | [进阶](collinear-layout.html#orientation) |
| 连接非相邻行 | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `adjacent_only, pairs, nonadjacent_linewidth` | B | [进阶](collinear-routing.html) |
| 显示局部窗口与断点 | [`segment_ggchord_data()`](reference/segment_ggchord_data.html) | `loci; coord_collinear: segments, segment_gap` | B | [进阶](views.html#segments) |
| 绘制固定长度标尺 | [`geom_seq_scale()`](reference/geom_seq_scale.html) | `length, location, offset` | B | [进阶](collinear-tuning.html#scale) |
| 改变圆图朝向 | [`coord_circular()`](reference/coord_circular.html) | `rotation, direction, gap, fit` | C | [基础](plasmid-tuning.html) |
| 堆叠特征注释 | [`position_feature_stack()`](reference/position_feature_stack.html) | `spacing, side, base_position` | A / B / C | [进阶](annotations.html) |
| 自动排布特征文字 | [`geom_feature_label_repel()`](reference/geom_feature_label_repel.html) | `label_layout, label_fit, label_min_size` | A / B / C (conditional) | [进阶](labels.html) |
| 搜索引物结合位置 | [`find_primer_bindings()`](reference/find_primer_bindings.html) | `set, circular` | DNA | [进阶](primers.html#input) |
| 绘制引物及文字 | [`geom_primer()`](reference/geom_primer.html) | `position` | C (bp) | [进阶](primers.html#primer-arcs) |
| 搜索和筛选酶切位点 | [`find_restriction_sites()`](reference/find_restriction_sites.html) | `methylation, enzymes` | DNA | [进阶](recognition.html#restriction) |
| 调整酶切文字引线 | [`geom_restriction_site()`](reference/geom_restriction_site.html) | `leader, label_style, label_side` | C (bp) | [进阶](restriction-sites.html#restriction-labels) |
| 绘制 GC 柱图与折线 | [`geom_track_bar()`](reference/geom_track_tile.html) | `track_value, baseline` | A / B / C | [进阶](tracks.html) |
| 绘制 tile 色带 | [`geom_track_tile()`](reference/geom_track_tile.html) | `track_fill, position` | A / B / C | [进阶](tracks.html#tile) |
| 设置轨道宽度和次序 | [`position_track()`](reference/position_track.html) | `track, side, order, width, gap` | A / B / C | [进阶](tracks.html#sharing) |
| 导入 FASTA、GFF3 和比对 | [`read_fasta_lengths()`](reference/read_fasta_lengths.html) | `sequence_unit` | general | [基础](formats.html) |
| 绘制原生蛋白比较 | [`ggchord()`](reference/ggchord.html) | `seq_data: sequence_unit` | A / B (aa) | [进阶](protein-alignments.html) |
| 预览并导出论文图 | [`view_ggchord()`](reference/view_ggchord.html) | `width, height, units` | general | [基础](output.html) |
| 诊断布局与导出坐标 | [`export_ggchord_layout()`](reference/export_ggchord_layout.html) | `plot` | general | [进阶](faq.html) |
| 让相同端点的连接共享主干 | [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | `link_branch, link_branch_fraction` | A / B | [进阶](chord-branches.html#branches) |
| 调整中心名称和长度 | [`geom_seq_center_label()`](reference/geom_seq_center_label.html) | `name, length, separator, name_style, length_style` | C | [进阶](circular-layers.html#center) |
| 固定标签方向与内外侧 | [`geom_feature_label()`](reference/geom_feature_label.html) | `label_orientation, label_side, label_rotation` | A / B / C | [进阶](circular-labels.html#fixed-labels) |
| 自动标签缩放与省略 | [`geom_feature_label_repel()`](reference/geom_feature_label_repel.html) | `label_fit, label_max_lines, external` | A / B / C | [进阶](circular-labels.html#fitting) |
| 选择连续与分段角色标度 | [`scale_ribbon_fill_gradientn()`](reference/scale_ribbon_fill_stepsn.html) | `colours, values, limits, guide` | A / B | [进阶](scale-families.html) |
| 诊断窗口和查询目录 | [`diagnose_ggchord()`](reference/diagnose_ggchord.html) | `plot, width, height; catalog: aliases` | general | [进阶](workflow-tools.html) |
| 查看全部接口的主讲位置 | [`ggchord()`](reference/ggchord.html) | `all exported APIs` | general | [进阶](api-curriculum.html) |
| 查找参数类型与完整默认值 | [`ggchord()`](reference/ggchord.html) | `all public formals` | general | [进阶](parameter-atlas.html) |

## 三条独立学习路线 {#paths}

弦图：骨架 → 比对与基因 → 布局 → 连接路径 → 标签 → 四基因组案例。共线性图：骨架 → 输入与比对 → 行距 → 对齐与方向 → 跨行与局部窗口 → 四基因组案例。圆图：骨架 → 特征与标签 → 样式 → 注释布局 → 引物与酶切 → 轨道 → 完整 pUC19 → 可选参考图库。

固定样式、映射、标度、图例和主题统一在[配色与图例](styling.html)主讲；[输出](output.html)集中说明预览、尺寸和保存。函数参考提供完整签名，任务表提供教学入口，并不表示所有参数已经达到同样教学深度。
