[English](MIGRATION.md) · [现代汉语](MIGRATION-Hans.md)

# 迁移到 ggchord 候选 API

## 范围与时间 {#scope}
下列写法是当前候选替代，不是 v1.0.0 长期兼容承诺。此前 0.16.x 版本均本地完成、
尚未发布。表中的旧调用用于解释迁移，不应作为当前例子执行。
历史证据只提供删除上界时，不猜测精确首删 tag。

## 替代表 {#table}
| 旧写法（勿执行） | 当前写法 | 删除 / 守卫 | 含义 |
| --- | --- | --- | --- |
| `geom_ribbon()` | `geom_link_ribbon()` | ≤ 0.16.1 | 连接图层改名，不影响 ggplot2::geom_ribbon()。 |
| `ggchord(..., rotation = 20)` | `ggchord(...) + coord_chord(rotation = 20)` | ≤ 0.16.1 | 全局坐标使用 + 组合。 |
| `geom_feature(feature_shape = "block", feature_width = .05)` | `geom_feature(shape = "block", width = .05)` | 0.16.1 | 映射角色 feature_shape 保留。 |
| `geom_gene(gene_offset = c("+" = .1, "-" = .2))` | `geom_gene(position = position_strand(c("+" = .1, "-" = -.2)))` | 0.16.1 | 显式 Position 值带符号。 |
| `geom_feature(feature_offset = .1)` | `geom_feature(position = position_strand(.1))` | 0.16.1 | 使用 Position 对象。 |
| `geom_gene(gene_width = .05)` | `geom_gene(width = .05)` | 0.16.1 | 固定几何统一使用 width。 |
| `geom_link_ribbon(ribbon_gap = .01)` | `geom_link_ribbon(link_gap = .01)` | 0.16.1 | 使用 link_gap，控制点同样使用 link_ctrl_point。 |
| `geom_gene_label_repel(gene_label_layout = "aligned")` | `geom_gene_label_repel(label_layout = "auto")` | 0.16.1 | 同时迁移参数名和取值。 |
| `geom_gene_label_repel(label_rotation = 90)` | `geom_gene_label(label_rotation = 90)` | 0.16.4 guard | 原先忽略的手动控制现报错，也可删除它使用自动布局。 |
| `geom_feature_label_repel(label_orientation = "horizontal")` | `geom_feature_label(label_orientation = "horizontal")` | 0.16.4 guard | 未知链 feature 使用固定 feature 标签层。 |
| `geom_restriction_site(leader = "trunk")` | `geom_restriction_site(leader = "radial")` | 0.16.1 | radial 默认值保持。 |
| `theme_ggchord(base_size = 10)` | `theme_ggchord() + theme(text = element_text(size = 10))` | ≤ 0.16.1 | 文字外观属于 theme 元素。 |

## 运行当前替代 {#replacement}
先运行[准备代码](README-Hans.md)及[圆图](README-Hans.md)中的 `features` 表：


```r
p_migrated <- ggchord(seqs[1, ]) + geom_seq() +
  geom_feature(data = features, shape = "block", width = .05) +
  geom_feature_label(data = features, label_rotation = 0) + coord_circular()
p_offsets <- ggchord(seqs, gene_data = genes) + geom_seq() +
  geom_gene(position = position_strand(c("+" = .1, "-" = -.2)))
```


旧负链 offset 幅度在显式 Position 中需要负号。映射仍使用
`aes(feature_shape = type)` 配合 shape scale，固定参数则为 `shape`。
自动标签使用 `label_layout`，手动几何控制使用相应固定标签层。
`label_layout = "aligned"` 改为 `"auto"`。已删除名称明确报错，不保留静默执行路径。

## 保留的支持 {#retained}
五个 colour/color 拼写别名、全部上游重导出、gene/feature/plasmid 快捷入口及 getter/export
两类几何接口均保留。各坐标默认差异和已记录结果字段保持。NULL/NA/缺参语义应查阅
API 契约 (`help("ggchord_api_contract", package = "ggchord")`)和当前参考；NULL 不统一表示关闭。
