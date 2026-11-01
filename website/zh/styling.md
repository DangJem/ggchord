## 从数据角色构建样式 {#roles}
先运行[共享准备代码](index.html#quick-start)。例子明确声明类别：ribbon 填色映射
query ID，序列颜色映射序列 ID，基因填色映射链。手动 scale 的名称与实际取值匹配。

{{code:styling}}
{{figure:styling}}

| 角色 | 映射数据 | 常用 scale 家族 |
| --- | --- | --- |
| `seq_colour`, `seq_ring` | 序列/类别 | 序列颜色或 ring 手动标度 |
| `ribbon_fill`, `ribbon_alpha`, `ribbon_colour`, `ribbon_linetype` | 比对元数据或统计 | 手动、渐变/分段/identity 和透明度标度 |
| `feature_fill`, `feature_shape`, `gene_fill`, `region_fill` | 注释类别或链 | 相应角色的手动/预设标度 |
| `link_colour`, `link_alpha`, `link_linewidth`, `link_linetype` | 点连接元数据 | 连续或手动 link 标度 |
| `track_value`, `track_fill`, `track_colour` | 命名轨道的数值/类别 | 命名数值域与颜色/填色标度 |

## 固定样式和用户 scale 优先 {#priority}
常量样式写在 `aes()` 外，变化样式写在映射中并配合合适的标度。
显式 geom 固定值优先，省略该固定参数才保持映射。用户 scale 优先于自动默认。
同类多层独立保留数据、映射和参数。按 `accver` 命名的几何向量，不表示 scale
也接受逐序列列表。

## 组织可读图例 {#guides}
`guide_ggchord_legend()` 控制离散角色图例，`guide_ggchord_colourbar()` 控制连续颜色
图例。colour/color 拼写别名保留，本文优先用 `colour`。图层 `show.legend`
决定显示，`guides()` 选择 guide，`theme(legend.position = ...)` 决定位置。
过宽图例可分行或增加导出宽度，不需要删减生物学数据来迁就。

## 为场景选择主题 {#themes}
例子在 `theme_variants` 中保存 minimal、dark 和 plasmid 备选，按场景添加到相应图中。
上图使用 `theme_ggchord_publication()`，plasmid 预设用于圆图流程。
主题控制非数据外观，不选择坐标系，也不替代数据驱动的填色 scale。
文字外观用 `theme(text = element_text(...))`；已删除的 `base_size` 式快捷参数见
[迁移指南](migration.html)。
