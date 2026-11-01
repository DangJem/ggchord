## 保留注释语义 {#stack}
先运行[共享准备代码](index.html#quick-start)和[圆图例子](modes.html#circular)。
这里保留 pUC19 中实际重叠的注释和原本未知的链方向。
几何和填色分别映射同一类别；显式固定 `shape` 会覆盖 shape 映射。

{{code:featurestack}}
{{figure:featurestack}}

`position_feature_stack()` 解决注释轨道重叠。标签层使用同一 Position，让文字锚点
跟随轨道。标签子集或重排时优先提供唯一 `feature_id`。严格的 gene 快捷层要求
已知 `+/-`，通用注释保留未知/双向语义。

## 按含义选择形状 {#shapes}
| 几何 | 常见用途 | 重要区别 |
| --- | --- | --- |
| `arrow`, `compact_arrow` | 有方向的编码/注释区间 | 已知方向来自证据，不由显示猜测 |
| `block`, `capped_line` | 区域或不填充区间 | 选择形状不推断生物学功能 |
| `promoter_arrow`, `primer_arrow`, `primer_arc` | 已声明的启动子/引物式注释 | 计算匹配与绘图分开 |
| `marker`, `chevron`, `lollipop` | 自定义区间/类点强调 | 保留真实区间元数据，显示形状属于映射 |

`geom_gene()`、`geom_feature_plasmid()` 是有用预设，默认值并不相同。
`position_strand()` 提供带符号位移，`position_plasmid()` 提供圆图放置预设，
后者不构成第四类坐标模式。

## 高亮区间而不改写观测 {#decoration}
{{code:decorations}}
{{figure:decorations}}

区域叠层有自己的数据和 `region_fill` 标度。`geom_seq_label()` 命名序列，
`geom_seq_center_label()` 在单圆图中心提供名称/长度，`geom_seq_break()` 标记裁剪或
分段末端。`geom_seq_scale()` 是仅用于共线性图的独立 bp 参照条，与控制序列基因组轴
的 `scale_seq_position_continuous()` 不同。普通 ggplot2 注释应使用
[位置映射](views.html#mapping)，不能把 bp 直接当作笛卡尔 x/y。
