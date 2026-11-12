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

## 按目标尺寸排版 {#adaptive}

添加 `adapt_ggchord()` 后才启用尺寸策略。此例使用随包记录的序列长度与真实 HSP；
参考画布单位为英寸，自动字号下限单位为 mm。诊断描述指定的 8 × 6 英寸 PNG；
标准 `ggsave()` 使用相同尺寸即可重现此排布。

```r
library(ggchord)
library(dplyr)
seqs <- seq_data_example |>
  select(accver, length)
plot <- ggchord(seqs, ribbon_data_example) +
  geom_seq() +
  geom_link_ribbon() +
  adapt_ggchord(
    mode = "adaptive",
    reference_size = c(8, 6),
    text_scale_range = c(0.85, 1.35),
    min_text_size = 2.5
  )
report <- diagnose_ggchord(plot, width = 8, height = 6, device = "png")
report$guide_layout$guides
report$adaptive_layout$styles
```

`fixed` 固定默认主题和图例部件的参考物理尺寸，`scale` 在倍率范围内调整主题驱动
的默认样式，`adaptive` 再联合选择原生图例的行列、侧边和方向并重测。它复用已有
标签布局；`feasible` 仅表示图例空间结果，不保证所有标签零碰撞。
不添加策略时保留既有绘图行为。

显式 geom 样式、主题字号、图例位置/方向/行列/键尺寸、用户标度及 `show.legend`
优先。连续和分段色条保持完整，类别及源记录不丢弃。约束超出画布预算时，诊断返回
`insufficient_space` 并在渲染时警告；此时增加导出尺寸或调整显式图例设置。
有界搜索不保证全局最优。仅导出布局时会记录策略，最终图例测量需要目标设备。

在组合 Fusion 前向需要的子图添加策略；各面板使用自己的分配尺寸，外层分配与
跨面板图例另行控制。分面、任意自定义 grob、所有几何常量的整体缩放，以及字体
嵌入认证不在本策略范围内。

## 完成后检查 {#check-result}
同一 query 的连接颜色一致；连续相似度应使用连续标度。固定颜色放在 aes() 外，图例必须能说明映射。

<div class="next-step">继续：<a href="adjustments.html">调整图形 →</a></div>

各角色 manual/continuous/steps/identity、线型、透明度和主题的完整可运行家族例子见[标度家族](scale-families.html)。
