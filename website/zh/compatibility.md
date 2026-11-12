换坐标系时先确认单位、图层与参数。下表结合当前 Rd、源码契约与实际渲染；支持表示存在对应实现，不保证每个参数组合。三种布局可以使用原生 bp/aa 单位；基因、引物和酶切仍有核酸限制。

## 图层与三类坐标 {#layers}

| 图层与参考 | 弦图 | 共线性图 | 单序列圆图 | 条件与边界 |
| --- | --- | --- | --- | --- |
| [`geom_seq()`](reference/geom_seq.html) | 支持 | 条件支持 | 支持 | 共线性拒绝显式 radius、curvature、ring |
| [`geom_seq_label()`](reference/geom_seq_label.html) | 支持 | 条件支持 | 支持 | 共线性不接受显式径向标签 radius |
| [`geom_gene()`](reference/geom_gene.html) | 条件支持 | 条件支持 | 条件支持 | 要求核酸单位和源正负链，aa 用通用特征 |
| [`geom_feature()`](reference/geom_feature.html) | 支持 | 支持 | 支持 | 要求合法 ID、区间和一致的原生单位 |
| [`geom_feature_plasmid()`](reference/geom_feature_plasmid.html) | 支持 | 支持 | 支持 | 通用特征预设；跨模式使用须重查圆图 Position 默认值 |
| [`geom_gene_label()`](reference/geom_gene_label.html) | 支持 | 条件支持 | 支持 | 与基因共享 Position；共线性 radial/tangent 对应竖排/横排 |
| [`geom_feature_label()`](reference/geom_feature_label.html) | 支持 | 条件支持 | 支持 | 与特征共享 Position；共线性 radial/tangent 对应竖排/横排 |
| [`geom_gene_label_repel()`](reference/geom_gene_label_repel.html) | 支持 | 条件支持 | 支持 | 共线性拒绝 external=FALSE，按目标尺寸检查隐藏结果 |
| [`geom_feature_label_repel()`](reference/geom_feature_label_repel.html) | 支持 | 条件支持 | 支持 | 共线性拒绝 external=FALSE，检查适配与隐藏 |
| [`geom_link_ribbon()`](reference/geom_link_ribbon.html) | 支持 | 条件支持 | 不支持 | 共线性按显示行判断邻接；分段同行带省略并报告 |
| [`geom_link_line()`](reference/geom_link_line.html) | 支持 | 条件支持 | 不支持 | 要求点坐标；不支持显式分段同行点连接 |
| [`geom_seq_region()`](reference/geom_seq_region.html) | 支持 | 支持 | 支持 | 要求合法源区间、ID 和局部侧向 |
| [`geom_seq_scale()`](reference/geom_seq_scale.html) | 不支持 | 支持 | 不支持 | 直线标尺仅适用共线性模式 |
| [`geom_seq_break()`](reference/geom_seq_break.html) | 不支持 | 条件支持 | 不支持 | 要求共线性窗口边界，不修改源序列 |
| [`geom_seq_center_label()`](reference/geom_seq_center_label.html) | 条件支持 | 不支持 | 条件支持 | 仅一条显示序列，长度使用其原生 bp/aa 单位 |
| [`geom_primer()`](reference/geom_primer.html) | 条件支持 | 不支持 | 条件支持 | 仅一条显示 DNA 序列及合法结合坐标 |
| [`geom_primer_label_repel()`](reference/geom_primer_label_repel.html) | 条件支持 | 不支持 | 条件支持 | 仅一条显示 DNA 序列，与引物共享位置 |
| [`geom_restriction_site()`](reference/geom_restriction_site.html) | 条件支持 | 不支持 | 条件支持 | 仅一条显示 DNA 序列及实际计算的切割坐标 |
| [`geom_track()`](reference/geom_track.html) | 支持 | 支持 | 支持 | 源空间 Stat/Position；见[首批矩阵](foundation.html#tracks) |
| [`geom_track_bar()`](reference/geom_track.html) | 支持 | 支持 | 支持 | 统计观测数，给定数值请使用 col |
| [`geom_track_col()`](reference/geom_track.html) | 支持 | 支持 | 支持 | 显式 position_track、数值高度与数值域 |
| [`geom_track_line()`](reference/geom_track.html) | 支持 | 支持 | 支持 | 显式 position_track；闭合须明确请求，不跨缺失样本连接 |
| [`geom_track_tile()`](reference/geom_track.html) | 支持 | 支持 | 支持 | 显式 position_track；不要求数值高度 |

23 个代表性图层分别在三种主要模式实际渲染，共 69 组明确标注的合成几何输入：59 组成功，10 组按契约拒绝，包含自动引物标签。中心标签在单序列弦图及圆图支持 bp/aa；引物与酶切要求 DNA。帮助已澄清这些边界。此处验证模式边界，不表示所有参数组合、生物识别或设备均已验收。参见[完整运行记录](../assets/source/tools/teaching-audit/layer-compatibility.csv)。

## 读懂条件支持 {#conditions}

`R/coord-collinear.R` 的 `ggchord_validate_coord_contract()` 显式检查圆图参数、自动方向与分段冲突。通用数据契约和核酸专属限制由运行时校验负责。完整源码随网站放在 `assets/source/`；函数参考从当前 roxygen/Rd 生成。

Position 只控制局部放置，不把不支持的图层转换成另一种坐标图层。共线性图的 outside/inside 是上方/下方，圆图是局部法线的外/内；复制配置前应重新读图。显示翻转与基因生物学链的区别见[四种方向](collinear-layout.html#orientation)。

## 版本与复现 {#versions}

支持版本见[安装页](installation.html)，更新现有脚本时请参阅[迁移说明](../assets/MIGRATION-Hans.md)。单圆图只接受单条显示序列；圆形蛋白只是显示布局，不声明蛋白分子为圆形。混合 bp/aa 数据必须先投影或分别作图。


## 新增坐标系 {#new-views}

新增契约见[Dotplot/Radial支持矩阵](view-expansion.html#support)。Dotplot使用独立identity HSP层；Radial将已有Feature与Track框架接入环形邻接。圆图专有层、共享直线比例尺和branch会明确报错。新增图层检查与此前三视图证据分开记录。


| 图层与参考 | Dotplot | Radial |
| --- | --- | --- |
| [`geom_alignment_segment()`](reference/geom_alignment_segment.html) | 支持identity HSP线段 | 不支持 |


## 静态跨面板组合 {#fusion}

[Fusion](fusion.html)接受Circular/Chord/Collinear中央图与最多四个Collinear外围面板，各面板保留原有图层契约。跨面板完整前缘采用一基闭区间边界和连续曲线映射；要求同单位来源区间完整可见，面板与实例路由无歧义。Fusion的标准标度只作用于跨面板图层。外围Sankey/Radial/Dotplot、分面及隐式第三方坐标适配均拒绝。保存时重新测量物理面板，不代表具备全局避障或通用字号缩放。


## 局部对齐字符单元 {#aligned-cells}

| 图层/参考 | Chord / Collinear / Circular / Dotplot / Radial | 普通 Cartesian ggplot |
| --- | --- | --- |
| [`geom_residue()`](reference/geom_residue.html) | 不支持生物学骨架投影 | 支持局部字符单元及 identity/reverse 位置标度 |

[字符对齐工作流](formats.html#aligned-characters) 使用普通 ggplot x/y 坐标，源位置和对齐列保留为不同字段。按最终物理字形度量选择可读文字或明确的色格/报错策略，不改变源数据、不缩小字号。不代表已支持 Fusion 面板注册、自动 CDS 翻译或新比对计算。

## 可选目标尺寸策略 {#adaptive}

[`adapt_ggchord()`](styling.html#adaptive) 支持单个 Chord、Circular、Collinear、
Sankey、Radial、Dotplot 图及显式启用策略的 Fusion 子面板。它调整默认主题样式，
在用户显式约束内排列原生图例；默认不启用。分面、任意自定义 grob、Fusion 外层
区域分配和全局标签求解器不在此范围内。
