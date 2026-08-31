# ggchord 设计路线图

本文档记录 ggchord 的实际实现状态、已确认的问题、目标 API 以及版本计划。
状态以当前代码为准，不再把“已有设计”标记成“已经完成”。

当前正式版本为 **v0.9.0**，代码库版本为 **v0.10.0 开发版**。
以正确性和 ggplot2 语法一致性为中心的基础重构已于 2026-08-31
随 v0.9.0 发布；v0.10.0 在用户确认前不创建正式标签或 Release。

总体优先级：

1. 正确性与可重复性；
2. ggplot2 风格的公开接口；
3. 绘图表现力；
4. 导出、性能和生态稳定性。

---

## 一、版本计划

| 版本 | 主题 | 状态 |
| --- | --- | --- |
| v0.7.0 | 数据验证、清理和基础测试 | 已发布；视觉回归并未真正建立 |
| v0.8.0 | 数据导入、ribbon 预处理、序列分组 | 已发布；遗留正确性问题已在 v0.9.0 修复 |
| v0.9.0 | grammar 基础、scale/theme/guide/coord、图层独立性 | 已正式发布（2026-08-31） |
| v0.10.0 | 发表级默认视觉、固定标签、feature 形状、显式多环和密集 ribbon | 开发中 |
| v0.11.0 | 布局导出和大数据性能 | 规划中 |
| v1.0.0 | API 冻结、完整文档和长期兼容承诺 | 规划中 |

### v0.9.0 正式版状态

v0.9.0 已完成：

- 数据清理、验证和导入中的已知正确性修复；
- plot-owned 布局、稳定 `layer_id` 和逐图层 geometry registry；
- 全部公开 geom 的 `data` / `mapping` 接入及同类图层隔离；
- role-specific aesthetic、scale、theme、guide 和 `coord_chord()`；
- `aligned`、`radial`、`arc` 三种确定性基因标签布局；
- 首批静态默认样式、静态设备验收及 Plotly 实验代码清理；
- 精简 testthat、roxygen/Rd 一致性检查和 `R CMD check`。

以下内容不是 v0.9.0 发布阻断项：

- 旧 scale 类 geom 参数在 v0.9.0 中继续作为带警告的兼容入口，到 v1.0.0 删除；
- 自定义 Stat/Geom 的进一步拆分属于内部演进，不改变 v0.9.0 的公开契约；
- 全量文档、网页和图片重写留到 v1.0.0；
- feature shape、多环、密集 ribbon 和正式布局导出按后续版本推进。

### v0.9.0 正式版验收记录（2026-08-31）

- testthat 保持为 15 个公开行为场景、90 项断言：0 fail、0 warning、0 skip；
- `R CMD check` 包含中英文 vignette 重建：0 error、0 warning；唯一 NOTE 为隔离
  环境无法联网校验系统时间，与包代码无关；
- 用户四序列示例在 `aligned`、`radial`、`arc` 下分别以 6×4 和 12×8 英寸构建，
  均为 0 标签重叠、0 指示线交叉、0 二次裁切损失；
- 相同设备和输入重复构建的标签与线段数据完全一致；
- 80 标签 `aligned` 基准的三次构建中位数约 0.36 秒，低于 2.5 秒验收线；
- PNG、PDF、SVG 均成功导出；默认、minimal、dark、publication 和灰度效果已
  人工检查；临时 SVG 验收工具不加入包依赖；
- 中英文 README 没有因本轮一般功能收尾发生变化。

---

## 二、设计原则

### 2.1 参数职责

遵循 ggplot2 的基本分工：

- `aes()` 决定数据列映射到什么视觉属性；
- `scale_*()` 决定映射后的取值、范围、变换、breaks、labels 和 guide；
- `geom_*()` 决定几何形状、固定样式和局部布局；
- `coord_chord()` 决定全局旋转、比例、坐标范围、裁切和自动适配；
- `theme_*()` 决定非数据元素的外观；
- `guide_*()` 决定图例的排列和绘制方式。

固定颜色或透明度仍可写在 geom 中。例如 `geom_ribbon(fill = "red")` 是固定
样式，不属于 scale。只有由数据映射而来的颜色、透明度、线型和图例配置必须
进入 `scale_*()`。

### 2.2 图层独立性

每个 ggchord 图层必须拥有独立的 `layer_id`、`data`、`mapping` 和参数。
同一张图应支持多个 region、feature、gene、highlight 或 label 图层，后添加的
图层不得覆盖前一图层的数据和样式。

### 2.3 默认效果

不添加自定义 scale、theme 或 coord 时，默认绘图效果应尽量保持不变。接口迁移
可以改变参数归属，但不应无理由改变默认配色、布局和标签位置。

### 2.4 兼容策略

- v0.9.0 保留旧的 scale 类参数并给出一次性弃用警告；
- 旧参数在内部转换成对应 scale；
- 同时使用旧参数和新 scale 时明确报冲突；
- v1.0.0 删除已经完成迁移的旧参数；
- 当前被静默忽略的 `data` 和 `mapping` 直接修复，不需要弃用过程；
- 已经完成精简的 `geom_gene_label_repel()` 不恢复旧的随机排斥参数。

---

## 三、v0.9.0 实现审计

本节同时保留问题来源与处理目标，便于解释 API 为什么发生变化。除明确写为
“后续项”的内容外，本节所列 v0.9.0 问题均已处理，不再视为待办。

### 3.1 核心构建流程

#### v0.9.0 前问题（已处理）

1. `compute_chord_geometry()` 每种图层只保存一份参数，多个同类图层由最后一个
   覆盖。两个 `geom_seq_region()` 最终会绘制完全相同的最后一组 region。
2. 多数 geom 虽然声明了 `mapping` 和 `data`，但实际使用固定占位数据和固定
   mapping，用户输入被忽略。
3. `ggplot_build.ggchord()` 会重建所有图层、插入 scale 并覆盖 coordinate，导致
   标准 ggplot2 扩展行为受到限制。
4. `get_chord_layout()` 返回包环境里“最近一次构建”的布局，多图交错构建时含义
   不明确。
5. `ggchord()` 重复执行基础验证和 `validate_ggchord_data()`，验证规则存在两个来源。

#### v0.9.0 处理结果

- 为每个图层分配稳定 `layer_id`；
- sequence 基础布局只计算一次；
- 其他图层按各自的 data、mapping 和 params 独立生成几何；
- 计算结果按 `layer_id` 注入对应图层；
- 几何数据保留原始数据列，允许标准 `aes()` 映射；
- 新增 `get_chord_layout(plot, build = TRUE)`，无参数调用先弃用再删除；
- 自定义 `ggplot_build()` 已缩小为布局注入和标准 build 桥接；进一步迁移到自定义
  Stat/Geom 是后续内部优化，不阻断 v0.9.0。

### 3.2 scale

v0.9.0 前，scale 在 `make_ggchord_scales()` 中临时创建，颜色、limits、breaks、
图例标题、图例位置和 key 尺寸分散在 geom 参数中。当前实现已公开角色专用
aesthetic 和 scale；默认 scale 只在用户没有提供相同角色 scale 时注入。

已公开的角色专用 aesthetic：

```r
seq_colour
group_colour
ribbon_fill
ribbon_alpha
ribbon_colour
ribbon_linetype
gene_fill
feature_fill
region_fill
```

首批公开 scale：

```r
scale_seq_colour_manual()
scale_seq_color_manual()       # 美式拼写别名
scale_group_colour_manual()
scale_group_color_manual()     # 美式拼写别名

scale_ribbon_fill_stepsn()
scale_ribbon_fill_gradientn()
scale_ribbon_fill_manual()
scale_ribbon_fill_identity()
scale_ribbon_alpha_continuous()
scale_ribbon_alpha_manual()
scale_ribbon_colour_manual()
scale_ribbon_color_manual()    # 美式拼写别名
scale_ribbon_linetype_manual()

scale_gene_fill_manual()
scale_feature_fill_manual()
scale_region_fill_manual()
```

所有 scale 已接受与 ggplot2 对应 scale 一致的常用参数：`name`、`breaks`、
`labels`、`limits`、`values`/`colours`、`na.value`、`oob`、`transform` 和 `guide`。
默认 scale 仅在用户没有添加对应角色 scale 时插入。

基因组坐标轴增加：

```r
scale_seq_position_continuous(
  breaks = waiver(),
  minor_breaks = waiver(),
  labels = waiver()
)
```

它负责各序列的位置 breaks、minor breaks 和标签格式。tick 长度、轴间距和文字
方向仍属于 geom/theme。

#### 参数迁移表

| 旧参数 | 新归属 |
| --- | --- |
| `seq_colors` | `scale_seq_colour_manual(values = ...)` |
| `seq_group_colors` | `scale_group_colour_manual(values = ...)` |
| `ribbon_colors` | ribbon fill scale |
| `ribbon_color_limits/breaks/name` | ribbon fill scale 的 `limits/breaks/name` |
| `ribbon_alpha_range` | ribbon alpha scale 的 `range` |
| `ribbon_outline_colors` | ribbon colour scale |
| `ribbon_linetypes` | ribbon linetype scale |
| `ribbon_direction_*` 视觉取值 | 对应的 alpha/colour/linetype scale |
| `gene_colors` | gene fill scale |
| `gene_order` | gene fill scale 的 `limits` 或 `breaks` |
| `feature_colors/order` | feature fill scale |
| region 分类颜色 | region fill scale |
| 轴 major/minor 数量和标签格式 | sequence position scale |

以下参数改为 `aes()`：

- `ribbon_color_by`；
- `ribbon_alpha_by`；
- `ribbon_outline_by`；
- `ribbon_linetype_by`；
- ribbon direction；
- gene strand/annotation；
- feature type/category；
- region category。

### 3.3 theme

v0.9.0 前 `ggchord()` 硬编码标题、边距、背景、网格和图例样式。当前默认外观已
提取为：

```r
theme_ggchord()
theme_ggchord_minimal()
theme_ggchord_dark()
theme_ggchord_publication()
```

`theme_ggchord()` 必须复现当前默认外观。只保留少量通用主题，不提供带特定期刊名称的主题。

可注册的 ggchord 专用主题元素：

```r
ggchord.axis.line
ggchord.axis.ticks
ggchord.axis.text
ggchord.seq.label
ggchord.group.label
ggchord.gene.label
ggchord.gene.label.segment
```

数据图层的填充色不进入 theme。显式 geom 参数覆盖 theme 默认值。

`ggchord()` 中以下参数逐步迁移：

- `title` → `labs(title = ...)`；
- `panel_margin` → `theme(plot.margin = ...)`；
- `show_legend` → `theme(legend.position = "none")`。

### 3.4 默认配色与发表级外观

v0.9.0 建立了视觉参数的 scale/theme/guide 基础，v0.10.0 继续校准默认
视觉，使用户不增加参数也能得到适合论文初稿的图形。这项工作不追求
装饰性，而强调可读性、层级、打印和色觉友好。

#### 配色原则

- sequence、group、gene 和 feature 的离散色板优先使用色觉友好且灰度可区分的
  配色；
- ribbon 连续色板使用亮度单调的方案，避免彩虹色板和难以解释的颜色跳变；
- strand、direction 等二分类采用对比明确、黑白打印仍可借助形状或线型识别的
  组合；
- region 和 highlight 使用可与 ribbon、gene 区分的强调色；
- 所有默认色板都允许由对应 `scale_*()` 完整替换；
- 不使用特定期刊、商业软件或机构的专有配色名称。

#### 默认参数与图例

- 重新校准 ribbon alpha、outline、sequence linewidth、gene 边框、标签字号、
  leader line 和图形留白；
- 默认图例按 sequence、ribbon、gene/feature、annotation 的阅读顺序排列；
- 优化 legend key glyph，使 sequence 方向、gene strand、ribbon fill/outline 在
  图例中与图中语义一致；
- 连续色条根据位置自动选择横向或纵向，并提供合适的长度、刻度和标题间距；
- 减少重复图例和无意义图例，默认标题使用简洁且可发表的文字；
- 默认白色背景适合 PDF/SVG 输出，同时保证用户主题可完全覆盖。

验收使用典型双序列、多序列、不同 radius/curvature、正反 orientation 和密集
标签示例，人工比较屏幕、PDF、SVG 和灰度打印效果。默认参数的变化必须记录在
NEWS 中。

### 3.5 guide

优先复用 ggplot2 的 `guide_legend()`、`guide_colourbar()` 和 `guides()`。新增的
辅助函数只做薄封装：

```r
guide_ggchord_legend()
guide_ggchord_colourbar()
```

封装只提供合适的 key glyph、默认方向和尺寸，不新建复杂 Guide ggproto。
新代码应使用 `guides()` 管理图例。`legend_position`、`legend_key_width`、
`legend_key_height` 在 v0.9.0 中仅保留为带警告的兼容入口，并计划在 v1.0.0
删除：

```r
guides(
  ribbon_fill = guide_ggchord_colourbar(position = "left"),
  gene_fill = guide_ggchord_legend(position = "right")
)
```

### 3.6 coord

旧 `coord_chord(layout = NULL)` 的 `layout` 未使用，构建时还会被新的
`coord_fixed()` 覆盖。当前已实现并采用以下接口：

```r
coord_chord(
  rotation = 45,
  ratio = 1,
  xlim = NULL,
  ylim = NULL,
  expand = FALSE,
  clip = "off",
  fit = c("labels", "geometry", "manual")
)
```

- 新代码由 coord 管理 `rotation`；`ggchord(rotation = ...)` 在 v0.9.0 中仅作兼容入口；
- 用户给出的 x/y limits 必须优先；
- `fit = "labels"` 使用当前设备感知的文字边界；
- `fit = "geometry"` 只适配几何数据；
- `fit = "manual"` 要求显式 limits；
- 自动范围对 x/y 分别紧贴内容，`coord_fixed()` 负责保持物理单位等比；
- 默认 `expand = FALSE`，因为自动范围已包含安全边距；
- 用户添加的 coord 不得被静默覆盖。

### 3.7 geom 参数和行为

下列 v0.9.0 接口整理均已完成；其中旧参数仍按 2.4 的兼容策略工作。
v0.10.0 不新增独立的手动标签 geom，而是增强现有 `geom_gene_label()`。

#### `geom_seq()`

- `data`、`mapping` 已生效；
- `seq_colors`、`seq_group_colors` 已迁入 scale，旧参数保留迁移警告；
- 序列显示文字由 `geom_seq_label(labels = ...)` 控制，scale labels 独立；
- group label 已有独立 `geom_seq_group_label()`，隐式旧行为暂时兼容；
- `linewidth` 和 arrow 已通过不同静态输出设备检查；
- `seq_order`、`seq_orientation`、`seq_gap`、`seq_radius`、
  `seq_curvature`、`seq_group` 和 `seq_group_gap` 保留为布局参数。

#### `geom_ribbon()`

- `ribbon_*_by` 已由角色 aes 接替，旧参数保留迁移警告；
- palette、limits、breaks、name、range 和分类取值已有对应 scale；
- 固定透明度可使用标准 `alpha`，旧 `ribbon_alpha` 仍兼容；
- outline 的新映射使用 `ribbon_colour` / `ribbon_linetype` 及对应 scale；
- direction 作为计算列保留，可由 aes/scale 显示；
- `ribbon_gap`、`ribbon_ctrl_point` 继续作为几何参数。

#### `geom_gene()` 与 `geom_feature()`

- gene/feature 分类列已可由 aes 指定，颜色和顺序由 scale 指定；
- `geom_feature()` 已统一为 `mapping = NULL, data = NULL` 的 ggplot2 风格签名；
- category 被 label 覆盖的问题已修复；
- 多个 gene/feature 图层已经独立；
- `geom_gene()` 中旧 label 参数会立即给出迁移错误；
- feature 的 block、chevron、lollipop 留到 v0.10.0，不标记为完成。

#### gene label

- `geom_gene_label()` 是唯一的固定位置和手工微调图层；
- `geom_gene_label_repel()` 已提供 `aligned`、`radial`、`arc` 三种确定性模式；
- 所有 label geom 的 `data`、`mapping` 已生效；
- leader 的 colour、linewidth、alpha 默认由 theme 元素控制；
- 固定文本 size 不再创建全局 `scale_size_identity()`；
- v0.10.0 新增水平/径向/切向文字、内外侧以及隐藏/微调/允许重叠策略；
- 逐序列、逐链的 rotation/radial/circumferential offset 继续提供手工自由；
- 不规划 `geom_gene_label_manual()`，避免两个固定标签图层职责重叠。

#### `geom_axis()` 和 `geom_seq_label()`

- axis breaks/labels 已移入 `scale_seq_position_continuous()`；
- 轴线、tick 和文字已分别接收经过筛选的样式参数；
- 无意义的 axis `show_legend` 已删除并给出明确错误；
- sequence label 的 data/mapping 已生效；
- 文字大小不再污染其他文字图层的 size scale。

#### region 和 highlight

- `region_color` 已正确绘制；
- `region_side = "auto"` 已根据真实局部法线选择方向；
- category 可通过 `region_fill` aesthetic 和 scale 产生图例；
- region 数据示例和 mapping 语义已修正；
- highlight selection 参数已有类型、长度、范围和有限性检查；
- highlight legend 使用映射而不是常量填充；
- 多个 region/highlight 图层已经独立工作。

### 3.8 数据验证、清理和导入

#### v0.9.0 前问题（已处理）

- `clean_ggchord_data(unknown_id = "keep")` 对未知 gene ID 报错；
- ribbon 的默认 interval 排序会丢失反向比对方向；
- unknown ribbon ID 与 unknown gene ID 的严重等级不一致；
- 大型 duplicate group 被跳过时没有在报告中说明；
- `deduplicate_ggchord_ribbons(keep = "first")` 按坐标排序而非输入顺序；
- 去重代表替换后报告引用可能过期；
- filter 一行只能记录最后一个移除原因；
- merge 后附加统计列可能沿用第一行的过时值；
- `read_blast(format = "outfmt7")` 未真正校验格式；
- outfmt 7 被错误假设为固定 17 列；
- GFF3 未显式处理 `##FASTA`；
- 多文件导入未保留来源文件。

#### v0.9.0 处理结果

- 统一结构验证入口，`ggchord()` 只调用一套规则；
- 清理 ribbon 时保留原始 direction 或 `sstrand`；
- filter report 支持一行多个原因；
- deduplicate 明确定义有向/无向 pair 和 tie 规则；
- merge 增加附加字段汇总策略；
- outfmt 7 解析 `# Fields:`，custom 格式继续允许显式列名；
- GFF3 在 `##FASTA` 停止，并完善 attributes 解析；
- 多文件读取可选添加 `.source_file`。

### 3.9 代码与文档清理

- 已删除完全注释且自 v0.3.0 不再使用的 `R/helpers_plot.R`、Plotly 实验代码、
  随机排斥辅助函数和确认无调用的旧内部边距函数；
- 新增内部函数优先使用 `@noRd`；既有 internal Rd 不在 v0.9.0 做无意义的批量
  重写，统一留给 v1.0.0 文档重构；
- v0.9.x 只更新与正确性、破坏性变更和公开 API 直接相关的 NEWS、roxygen 和 Rd；
- 按用户确认保持中英文 README 原有主体内容，不加入一般功能更新；
- 全量 source 文档、vignette、网页手册和图片统一留到 v1.0.0 发布前重构，避免
  同一内容在开发阶段被反复编辑。

---

## 四、v0.9.0 — grammar 与正确性

**状态：已于 2026-08-31 正式发布。**

### A. 已确认 bug 修复

- 修复 clean、deduplicate、feature category、region outline/legend、BLAST format；
- 修复默认静态绘图中的明显样式问题；
- 为每个 bug 保留一个最小回归测试。

状态：已完成。

### B. 图层独立性

- layer ID 和逐层 geometry registry；
- data/mapping 生效；
- 多个同类图层互不覆盖；
- `get_chord_layout(plot)` 替代全局最近布局。

状态：已完成；无参调用仅作为 v0.9.0 兼容入口保留并警告。

### C. scale/theme/guide/coord

- 实现第三节列出的首批接口；
- 建立首批色觉友好默认色板、图例和主题元素；v0.10.0 继续视觉校准；
- 旧参数进入弃用期；
- 当前阶段只维护必要的 NEWS 和生成文档，README 保持稳定；完整文档重构留到 v1.0.0。

状态：已完成。README 按确认保持稳定，必要变更已记录在 NEWS 和生成的 Rd 中。

### D. 标签布局

- 保留已完成的 `aligned`、`radial`、`arc`；
- 继续保证不同 radius、curvature、gap、orientation、rotation 和设备尺寸；
- 标签实现迁入逐层架构，不恢复随机 force 参数。

状态：已完成；三种模式通过不同序列方向、半径、曲率、间距、旋转和静态设备
组合验收。

---

## 五、v0.10.0 — 表达能力与默认视觉

**状态：开发中；用户确认后再发布正式版。**

### A. 固定标签与发表级默认视觉

- 增强 `geom_gene_label()`，使其同时承担简洁默认标签和精确手工微调；
- Identity 色条使用紧凑的物理尺寸，不再随设备高度无限拉伸；
- sequence/gene 图例符号表达实际方向，并减少过粗的线条和箭头；
- strand 图例与 `geom_gene()` 使用一致的等宽箭身和渐缩箭头；
- `aligned` 标签轨道可向附近空闲扇区再平衡，并用受限的外侧角部侧列疏散偏斜聚集；
- 默认白底、标题、轴线、标签、ribbon 透明度和图例间距统一校准；
- 文档和默认验收图使用 4:3 画布；实际设备尺寸仍由 RStudio 或 `ggsave()` 控制；
- 用 4×3、6×4、8×6、12×8 英寸以及 PNG/PDF/SVG 验收。

视觉设计只借鉴通用原则：[Circos](https://genome.cshlp.org/content/19/9/1639)
的环形信息层级与克制 ribbon、[clinker](https://academic.oup.com/bioinformatics/article/37/16/2473/6129045)
的发表级基因箭头、[DNA Features Viewer](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/)
的注释冲突处理，以及 [SnapGene](https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map)
和 [Geneious](https://manual.geneious.com/en/latest/Sequences.html) 的局部特征标签与拥挤隐藏思路。
不复制第三方资产、专有配色或具体视觉实现。

### B. 通用 feature geometry

实现真正独立的 feature shape：

```r
geom_feature(aes(feature_shape = type))
scale_feature_shape_manual(values = c(
  CDS = "arrow", tRNA = "block", repeat = "chevron"
))
```

首批只考虑 `arrow`、`block`、`chevron`、`lollipop`。`geom_gene()` 保持为 gene
arrow 的便捷封装，但验收标准是几何和视觉等价，不承诺内部字节完全一致。

### C. 显式多环

不实现根据序列数自动猜测的 `seq_ring = "auto"`。如果实际需求充分，提供用户
显式的 ring 映射；每个 ring 的 radius、gap 和 ribbon 锚点必须可解释、可导出。

### D. 密集 ribbon

不在 `geom_ribbon()` 中加入 `ribbon_reduce`。使用职责更明确的接口：

```r
stat_ribbon_bundle()
stat_ribbon_density()
bundle_ggchord_ribbons()
```

抽样属于数据处理，应使用显式预处理函数，而不是由 geom 静默丢弃数据。

---

## 六、v0.11.0 — 导出与性能

### A. 布局导出

只保留一个正式名称：

```r
export_ggchord_layout(
  plot,
  include = c("seq", "ribbon", "gene", "feature", "labels", "axis"),
  original_data = FALSE
)
```

返回带 `layer_id`、`source_row` 和原始映射列的布局副本。坐标契约必须记录 rotation、
ratio、单位和是否已经应用 coord transform。不再规划重复 alias。

### B. 性能

- sequence 基础布局缓存一次；
- 每层只计算自己的几何；
- 大型 ribbon 的 bundle/density 有独立 benchmark；
- 不使用影响绘图结果的全局可变缓存。

---

## 七、v1.0.0 — 稳定 API

- 删除已完成迁移的旧参数；
- 冻结核心 aes、scale、theme、guide、coord 和 layout export 契约；
- 基于最终 v1.0.0 API 重写全部中英文说明文档、README、vignette、roxygen、Rd、
  网页手册、FAQ 和迁移指南；
- 删除旧的示例、旧图片及过时表述，从最终 API 重新编写最小且连贯的示例；
- 重新生成所有带图片输出的示例，清理不再被引用的 `man/figures` 和网页图片；
- 文档结构遵循“安装 → 数据格式 → 最小出图 → 常用调整 → 进阶功能 → FAQ”，
  让新用户能够快速上手；
- README 继续保持稳定简洁，把详细功能说明留给 vignette 和 reference；
- CRAN check、最小视觉快照和跨平台检查；
- 明确支持的 R 和 ggplot2 版本；
- 无严重已知正确性问题后才发布 1.0.0。

交互式绘图不属于 v0.9.0–v1.0.0 的发布范围。相关实现、依赖、说明和网页入口在
正式版前从发行内容中清理；v1.0.0 发布后再基于稳定的 scale、guide、coord 和
layout export 契约，单独评估更成熟的交互方案，不在本路线图中预先承诺接口。

---

## 八、testthat 精简策略

测试只保留公开 API 的最小行为，不再通过大量精确坐标断言锁定内部实现。

### 保留

- `ggchord()` 创建对象并能完成一次基础 build；
- 核心 geom 组合能 build；
- 三种 label layout 能 build；
- validate、clean、三个 import、三个 ribbon utility 的单一成功路径；
- region、highlight、feature、sequence group 的单一成功路径；
- 每个已确认严重 bug 修复后保留一个最小回归用例；
- 一个最基本的错误输入测试。

### 删除

- 直接测试未导出的内部函数；
- 精确浮点坐标、文字框尺寸和迭代次数；
- 多个设备尺寸的重复组合；
- guide ggproto 内部字段；
- 同一参数的多种等价输入格式穷举；
- 重复构建、序列化、文件输出等低价值集成测试；
- 易受机器性能影响的固定秒数测试；
- 与实现细节耦合的内部 scale 标记测试。

测试目标是快速发现“公开函数完全不可用”，而不是证明所有排列组合。复杂布局的
质量主要通过少量视觉示例、人工验收和专门 benchmark 检查，不塞入默认 testthat。

---

## 九、验收规则

- 精简后的 `testthat::test_local()` 必须通过；
- `R CMD check` 必须通过；
- 新的公开参数至少有一个可运行文档示例；
- 破坏性变化必须出现在 NEWS 和迁移表；
- 默认示例在静态图中的外观不发生无理由变化；
- 多个同类图层、正反 orientation、不同 radius/curvature/gap 是核心人工验收场景；
- v1.0.0 前的开发文档避免无必要的反复改写。
