<a id="api"></a>

## 接口迁移 {#api}

自 v0.50.0 起，序列布局、注释分道、轨道和图例使用重构后的接口。旧脚本可以逐层调整：序列的排列和方向由 `layout_sequences()` 指定，特征分道由 `position_feature_stack()` 指定，定量图层可用 `track_band()` 共用轨道。

[完整迁移表与替代写法](../../MIGRATION-Hans.md) · [按用途查函数](../../man/)

<a id="pages"></a>

## 旧文档地址 {#pages}

原有绘图页、案例和图像锚点仍可通过收藏的地址访问。从[绘图指南](modes.md)可按图形查找，从[数据与操作](formats.md)可按工作步骤查找。

## 已移除的图层与位置函数 {#removed-layers}

- 将 `geom_gene()` 和 `geom_feature_plasmid()` 分别改为 `geom_feature(preset = "gene")` 和 `geom_feature(preset = "plasmid")`，继续传入对应的图层数据。[特征预设](../../man/geom_feature.Rd#arg-preset)
- 将 `geom_feature_label_repel()` 改为 `geom_feature_label(placement = "repel")`。用 `label_layout` 设置放置策略，用 `overflow` 设置标签过长时的处理方式。需要指定实体层时使用 `anchor_layer`。[特征标签](../../man/geom_feature_label.Rd)
- 将 `position_plasmid(offset)` 改为 `position_offset(offset)`。需要分道时用 `position_feature_stack()` 指定起始位置和扩展方向。[偏移](../../man/position_offset.Rd) · [特征分道](../../man/position_feature_stack.Rd)
- 用 `track_band("depth", width = w)` 声明物理带，在各定量图层中传入 `track = "depth"`。旧的 `position_track()` 曾负责轨带的物理放置，新的 `position` 负责图层内数据的堆叠、并排或微移。[轨带](../../man/track_band.Rd) · [定量图层](../../man/geom_track.Rd)
