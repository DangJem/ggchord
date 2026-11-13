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

{{code:courseblockw}}

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

<div class="next-step">继续：<a href="styling.html#roles">调整图形 →</a></div>

各角色 manual/continuous/steps/identity、线型、透明度和主题的完整可运行家族例子见[标度家族](scale-families.html)。

## 使用标准 theme 细调外观 {#native-theme}

`ggchord()` 已自动提供 `theme_ggchord()`。只有场景需要时才切换四个预设：Minimal简化装饰，Publication采用紧凑字体，Dark调整深色背景上的非数据装饰，Plasmid强调圆形注释。数据palette、科学连接及序列方向保持独立；背景需要其他数据对比度时，应明确使用数据scale或Geom固定颜色。

| 要调整的内容 | 所属组件 |
| --- | --- |
| 背景、字体、非数据标签/引线外观 | `theme()` |
| 序列颜色、Feature填色、Track数值范围/刻度/单位 | `aes()`及相应Scale |
| 骨架形态、Ribbon箭头/渐变模式、数据对象固定样式 | Geom |
| Track物理band、Feature堆叠与偏移 | Position |
| 序列顺序/方向及投影 | Coord / 显式数据流程 |
| 图例行列、顺序、方向、单个位置 | Guide；全局摆放用Theme |

注册的 `ggchord.axis.*` 物理设置接受单个绝对单位，与 `theme_ggchord()` 现有短参数使用同一配置。后续细调优先写 `theme(ggchord.axis.gap = unit(1, "mm"))`。一次主题构造只使用一种拼写，同义入口同时提供会报错；后续 `theme()` 按熟悉的ggplot2逐字段合并。`ggchord.seq.label`、`ggchord.gene.label`、`ggchord.feature.label`及其 `.segment` 线元素控制真实注释图层，显式Geom样式/映射仍优先。普通Cartesian的 `axis.text.x` 不代替自动序列轴；Dotplot保留独立原生x/y轴。

### Track共享外观与局部覆盖 {#track-theme}

以下真实GC比例来自未改动的随包pUC19序列，需可选dplyr。Track数据Scale负责比例值域和单位，其局部原生Guide负责轴装饰，普通Theme负责共享外观，没有虚构测序深度或read support。

{{code:courseblockx}}

`ggchord.track.axis.line/text/title` 控制共享数值轴外观；`.text.offset` 与 `.title.offset` 使用非负绝对单位。预留图幅时采用解析后的字号/字体/角度实测。GC Guide的显式文字样式优先于后来添加的全局9pt默认；局部 `element_blank()` 只隐藏对应装饰。`guide_axis()` 支持局部 `axis.line/text/title`、`angle` 与 `check.overlap`，其余选项明确报错，不静默忽略。

同Track ID的图层共享band、数值域和一个数值参考轴；不同ID已有独立数值/fill/colour Scale和颜色Guide。按ID添加 `scale_track_*()`，无需为每个Track注册主题元素。颜色Guide的局部外观可以通过对应Track色标中的原生 `guide_legend(theme = ...)` 或 `guide_colourbar(theme = ...)` 设置；全局图例字体/背景/位置仍由Theme负责。主题切换不改变序列顺序、HSP筛选、统计或单位。标签/Guide自动排版是独立策略，会按最终导出尺寸重新测量。

<div class="task-intro"><p>把常见调整放到正确的 ggplot2 组件：数据映射、几何位置、坐标和主题各负其责。</p></div>

<span id="mapping-prepare"></span>

## 职责分配 {#grammar}
[绘图基本概念](modes.html#prepare)解释图层、映射、标度和坐标的分工。下面只演示如何把这种分工用到比对相似度。

## 映射连接相似度 {#mapping}


{{code:adjustments}}

{{figure:mappingexample}}

## 放置注释 {#position}
使用 `position_strand(c("+" = .1, "-" = -.2))` 指定带符号的逐链位移。
`position_feature_stack()` 处理重叠注释，`position_track()` 分配命名定量轨道。
逐序列几何参数可用单值、等长无名向量或按 `accver` 命名的向量；推荐命名向量。
这套语法不扩展到标度、图例或 theme。

每层独立拥有 `data`、`mapping` 和参数，同类注释/连接层不会互相替换。
`show.legend` 和 `inherit.aes` 遵循 ggplot2。输入角色映射在几何生成前处理，
后期视觉映射遵循 ggplot2 求值。详见 [接口说明](reference/ggchord_api_contract.html)。

## 完成后检查 {#check-result-supplement}
旋转整体图形使用 coord；逐序列几何使用明确命名值；颜色随数据变化时使用映射和标度。
