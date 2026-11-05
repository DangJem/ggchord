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
| [`geom_track_bar()`](reference/geom_track_tile.html) | 支持 | 支持 | 支持 | 显式 position_track、数值高度与数值域 |
| [`geom_track_line()`](reference/geom_track_tile.html) | 支持 | 支持 | 支持 | 显式 position_track；闭合须明确请求，不跨缺失样本连接 |
| [`geom_track_tile()`](reference/geom_track_tile.html) | 支持 | 支持 | 支持 | 显式 position_track；不要求数值高度 |

21 个 geom 分别在三种模式实际渲染，共 63 组明确标注的合成几何输入：53 组成功，10 组按契约拒绝，包含自动引物标签。中心标签在单序列弦图及圆图支持 bp/aa；引物与酶切要求 DNA。帮助已澄清这些边界。此处验证模式边界，不表示所有参数组合、生物识别或设备均已验收。参见[完整运行记录](../assets/source/tools/teaching-audit/layer-compatibility.csv)。

## 读懂条件支持 {#conditions}

`R/coord-collinear.R` 的 `ggchord_validate_coord_contract()` 显式检查圆图参数、自动方向与分段冲突。通用数据契约和核酸专属限制由运行时校验负责。完整源码随网站放在 `assets/source/`；函数参考从当前 roxygen/Rd 生成。

Position 只控制局部放置，不把不支持的图层转换成另一种坐标图层。共线性图的 outside/inside 是上方/下方，圆图是局部法线的外/内；复制配置前应重新读图。显示翻转与基因生物学链的区别见[四种方向](collinear-layout.html#orientation)。

## 版本与复现 {#versions}

本表面向开发版 0.16.7.9000；[安装页](installation.html)说明本地源码安装与正式版入口。单圆图只接受单条显示序列；圆形蛋白只是显示布局，不声明蛋白分子为圆形。混合 bp/aa 数据必须先投影或分别作图。
