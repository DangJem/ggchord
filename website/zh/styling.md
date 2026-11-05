<div class="task-intro"><p>让颜色表达明确含义：先决定类别或连续数值，再选择映射、标度和图例。网站强调色不决定数据配色。</p></div>

<span id="prepare"></span>

## 从数据角色构建样式 {#roles}

使用四条弦图骨架和一条实际比对。把连接带填色映射到 query ID，再用一个命名颜色值给该类别设色。序列颜色继续由默认序列标度处理。

{{code:styling}}
{{figure:styling}}

| 角色 | 映射数据 | 常用标度家族 |
| --- | --- | --- |
| `seq_colour`, `seq_ring` | 序列/类别 | 序列颜色或 ring 手动标度 |
| `ribbon_fill`, `ribbon_alpha`, `ribbon_colour`, `ribbon_linetype` | 比对元数据或统计 | 手动、渐变/分段/identity 和透明度标度 |
| `feature_fill`, `feature_shape`, `gene_fill`, `region_fill` | 注释类别或链 | 相应角色的手动/预设标度 |
| `link_colour`, `link_alpha`, `link_linewidth`, `link_linetype` | 点连接元数据 | 连续或手动 link标度|
| `track_value`, `track_fill`, `track_colour` | 命名轨道的数值/类别 | 命名数值域与颜色/填色标度 |

## 固定样式和用户标度优先 {#priority}
常量样式写在 `aes()` 外，变化样式写在映射中并配合合适的标度。
显式几何层固定值优先，省略该固定参数才保持映射。用户标度优先于自动默认。
同类多层独立保留数据、映射和参数。按 `accver` 命名的几何向量，不表示标度
也接受逐序列列表。

## 组织可读图例 {#guides}
`guide_ggchord_legend()` 控制离散角色图例，`guide_ggchord_colourbar()` 控制连续颜色
图例。colour/color 拼写别名保留，本文优先用 `colour`。图层 `show.legend`
决定显示，`guides()` 选择图例，`theme(legend.position = ...)` 决定位置。
过宽图例可分行或增加导出宽度，不需要删减生物学数据来迁就。

## 为图形选择主题 {#themes}
本页使用默认主题。需要调整文字、图例或背景时，在图后添加 `theme()`；例如 `theme(text = element_text(size = 12))` 设置基础字号。`theme_ggchord_publication()` 提供出版样式，`theme_ggchord_plasmid()` 用于单序列圆图。这些主题控制非数据外观，不改变坐标或数据配色。

## 完成后检查 {#check-result}
同一 query 的连接颜色一致；连续相似度应使用连续标度。固定颜色放在 aes() 外，图例必须能说明映射。

<div class="next-step">继续：<a href="adjustments.html">调整图形 →</a></div>

各角色 manual/continuous/steps/identity、线型、透明度和主题的完整可运行家族例子见[标度家族](scale-families.html)。
