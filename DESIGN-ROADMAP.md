# ggchord 设计路线图

本文档记录 ggchord 的实际实现状态、已确认的问题、目标 API 以及版本计划。
状态以当前代码为准，不再把“已有设计”标记成“已经完成”。

当前版本为 **v0.9.0 开发版**。v0.9.0 正式发布前的首要目标不是继续增加
图层数量，而是完成一次以正确性和 ggplot2 语法一致性为中心的基础重构。

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
| v0.8.0 | 数据导入、ribbon 预处理、序列分组 | 已发布；仍有若干正确性问题待修复 |
| v0.9.0 | grammar 基础、scale/theme/guide/coord、图层独立性 | 开发中 |
| v0.10.0 | feature 形状、手动标签、显式多环和密集 ribbon | 规划中 |
| v0.11.0 | 布局导出和大数据性能 | 规划中 |
| v1.0.0 | API 冻结、完整文档和长期兼容承诺 | 规划中 |

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

## 三、当前实现审计

### 3.1 核心构建流程

#### 已确认问题

1. `compute_chord_geometry()` 每种图层只保存一份参数，多个同类图层由最后一个
   覆盖。两个 `geom_seq_region()` 最终会绘制完全相同的最后一组 region。
2. 多数 geom 虽然声明了 `mapping` 和 `data`，但实际使用固定占位数据和固定
   mapping，用户输入被忽略。
3. `ggplot_build.ggchord()` 会重建所有图层、插入 scale 并覆盖 coordinate，导致
   标准 ggplot2 扩展行为受到限制。
4. `get_chord_layout()` 返回包环境里“最近一次构建”的布局，多图交错构建时含义
   不明确。
5. `ggchord()` 重复执行基础验证和 `validate_ggchord_data()`，验证规则存在两个来源。

#### 改进方向

- 为每个图层分配稳定 `layer_id`；
- sequence 基础布局只计算一次；
- 其他图层按各自的 data、mapping 和 params 独立生成几何；
- 计算结果按 `layer_id` 注入对应图层；
- 几何数据保留原始数据列，允许标准 `aes()` 映射；
- 新增 `get_chord_layout(plot, build = TRUE)`，无参数调用先弃用再删除；
- 中期逐步迁移到自定义 Stat/Geom，缩小自定义 `ggplot_build()` 的职责。

### 3.2 scale

当前 scale 在 `make_ggchord_scales()` 中临时创建，颜色、limits、breaks、图例标题、
图例位置和 key 尺寸又分散在各 geom 参数中。普通 `scale_fill_*()` 只能可靠影响
gene，不能独立控制 ribbon、region 和 feature。

目标是公开角色专用 aesthetic：

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

scale_ribbon_fill_stepsn()
scale_ribbon_fill_gradientn()
scale_ribbon_fill_manual()
scale_ribbon_fill_identity()
scale_ribbon_alpha_continuous()
scale_ribbon_alpha_manual()
scale_ribbon_colour_manual()
scale_ribbon_linetype_manual()

scale_gene_fill_manual()
scale_feature_fill_manual()
scale_region_fill_manual()
```

所有 scale 应接受与 ggplot2 对应 scale 一致的常用参数：`name`、`breaks`、
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

`ggchord()` 当前硬编码标题、边距、背景、网格和图例样式。目标是将默认外观提取为：

```r
theme_ggchord()
theme_ggchord_minimal()
theme_ggchord_dark()
theme_ggchord_publication()
```

`theme_ggchord()` 必须复现当前默认外观。只保留少量通用主题，不提供带特定期刊
名称的主题。

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

grammar 重构时同步调整默认视觉，使用户不增加参数也能得到适合论文初稿的图形。
这项工作不追求装饰性，而强调可读性、层级、打印和色觉友好。

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
`legend_position`、`legend_key_width`、`legend_key_height` 从 geom 删除，改为：

```r
guides(
  ribbon_fill = guide_ggchord_colourbar(position = "left"),
  gene_fill = guide_ggchord_legend(position = "right")
)
```

### 3.6 coord

当前 `coord_chord(layout = NULL)` 的 `layout` 未使用，构建时又被新的
`coord_fixed()` 无条件覆盖。目标接口：

```r
coord_chord(
  rotation = 45,
  ratio = 1,
  xlim = NULL,
  ylim = NULL,
  expand = TRUE,
  clip = "off",
  fit = c("labels", "geometry", "manual")
)
```

- `rotation` 从 `ggchord()` 迁入 coord；
- 用户给出的 x/y limits 必须优先；
- `fit = "labels"` 使用当前设备感知的文字边界；
- `fit = "geometry"` 只适配几何数据；
- `fit = "manual"` 要求显式 limits；
- 用户添加的 coord 不得被静默覆盖。

### 3.7 geom 参数和行为

#### `geom_seq()`

- `data`、`mapping` 当前无效；
- `seq_colors`、`seq_group_colors` 移入 scale；
- `seq_labels` 拆分为序列显示文字和图例 labels，避免一参两用；
- group label 建议成为独立的 `geom_seq_group_label()`；
- `linewidth` 和 arrow 应在不同静态输出设备中保持一致；
- 保留 `seq_order`、`seq_orientation`、`seq_gap`、`seq_radius`、
  `seq_curvature`、`seq_group` 和 `seq_group_gap` 等布局参数。

#### `geom_ribbon()`

- `ribbon_*_by` 改为 aes；
- palette、limits、breaks、name、range 和分类取值移入 scale；
- `ribbon_alpha` 与 `alpha` 合并为标准 `alpha`；
- outline 参数改用 `colour`、`linewidth`、`linetype`；
- direction 始终作为计算列提供，由 aes/scale 决定是否展示；
- `ribbon_gap`、`ribbon_ctrl_point` 保留为几何参数。

#### `geom_gene()` 与 `geom_feature()`

- gene/feature 分类列由 aes 指定，颜色和顺序由 scale 指定；
- `geom_feature()` 统一成 `mapping = NULL, data = NULL` 的 ggplot2 风格签名；
- 修复 category 被 label 覆盖的问题；
- 多个 gene/feature 图层必须独立；
- `geom_gene()` 中已迁移的旧 label 参数应从 `...` 真正移除；
- feature 的 block、chevron、lollipop 尚未实现，不再标记为完成。

#### gene label

- 保留 `geom_gene_label()` 作为固定位置和手工微调图层；
- 保留 `geom_gene_label_repel()` 的 `aligned`、`radial`、`arc` 三种确定性模式；
- 修复所有 label geom 的 `data`、`mapping`；
- leader 的 colour、linewidth、alpha 默认由 theme 元素控制；
- 文本 size 作为固定值时不再映射到全局 `scale_size_identity()`；
- 真正的逐标签手动坐标契约放入 v0.10.0。

#### `geom_axis()` 和 `geom_seq_label()`

- axis breaks/labels 移入 `scale_seq_position_continuous()`；
- 轴线、tick 和文字分别接收样式，不能共享同一个未筛选的 `...`；
- 删除无意义的 axis `show_legend`；
- sequence label 的 data/mapping 必须生效；
- 文字大小不再污染其他文字图层的 size scale。

#### region 和 highlight

- 修复 `region_color` 未绘制；
- `region_side = "auto"` 必须根据真实局部法线选择方向；
- category 通过 `region_fill` aesthetic 和 scale 产生图例；
- 修正文档中把 region 数据误传给 mapping 的示例；
- highlight 的 selection 参数增加类型、长度、范围和有限性检查；
- highlight legend 使用映射而不是常量填充；
- 多个 region/highlight 图层独立工作。

### 3.8 数据验证、清理和导入

#### 已确认问题

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

#### 改进方向

- 统一结构验证入口，`ggchord()` 只调用一套规则；
- 清理 ribbon 时保留原始 direction 或 `sstrand`；
- filter report 支持一行多个原因；
- deduplicate 明确定义有向/无向 pair 和 tie 规则；
- merge 增加附加字段汇总策略；
- outfmt 7 解析 `# Fields:`，custom 格式继续允许显式列名；
- GFF3 在 `##FASTA` 停止，并完善 attributes 解析；
- 多文件读取可选添加 `.source_file`。

### 3.9 代码与文档清理

- 删除完全注释且自 v0.3.0 不再使用的 `R/helpers_plot.R`；
- 删除已被确定性标签算法取代的随机排斥辅助函数；
- 内部函数优先使用 `@noRd`，减少无意义的 internal Rd 页面；
- 删除未使用变量和重复赋值；
- README 中“Full ggplot2 integration”的表述在重构完成前改为更准确的说明；
- v0.9.x 只更新与正确性、破坏性变更和公开 API 直接相关的说明；
- 一般功能更新只写入 NEWS，不加入 README；
- README 保持简短稳定，只保留定位、安装、最小示例和文档入口；
- 全量 source 文档、roxygen、Rd、vignette、网页手册和图片统一留到 v1.0.0
  发布前重构，避免同一内容在开发阶段被反复编辑。

---

## 四、v0.9.0 — grammar 与正确性

### A. 已确认 bug 修复

- 修复 clean、deduplicate、feature category、region outline/legend、BLAST format；
- 修复默认静态绘图中的明显样式问题；
- 为每个 bug 保留一个最小回归测试。

### B. 图层独立性

- layer ID 和逐层 geometry registry；
- data/mapping 生效；
- 多个同类图层互不覆盖；
- `get_chord_layout(plot)` 替代全局最近布局。

### C. scale/theme/guide/coord

- 实现第三节列出的首批接口；
- 以发表级可读性为目标更新默认色板、图例和默认参数；
- 旧参数进入弃用期；
- 当前阶段只维护精简后的中英文 README 和必要 NEWS；完整文档重构留到 v1.0.0。

### D. 标签布局

- 保留已完成的 `aligned`、`radial`、`arc`；
- 继续保证不同 radius、curvature、gap、orientation、rotation 和设备尺寸；
- 标签实现迁入逐层架构，不恢复随机 force 参数。

---

## 五、v0.10.0 — 表达能力

### A. 通用 feature geometry

实现真正独立的 feature shape：

```r
geom_feature(aes(feature_shape = type))
scale_feature_shape_manual(values = c(
  CDS = "arrow", tRNA = "block", repeat = "chevron"
))
```

首批只考虑 `arrow`、`block`、`chevron`、`lollipop`。`geom_gene()` 保持为 gene
arrow 的便捷封装，但验收标准是几何和视觉等价，不承诺内部字节完全一致。

### B. 手动标签布局

为逐标签覆盖提供独立数据契约，例如：

```r
geom_gene_label_manual(
  data = label_positions,
  aes(gene_id = gene_id, x = x, y = y, label = label)
)
```

自动模式的结果可以由 layout export 导出、手工调整后再输入。具体列名和优先级在
实现前单独设计。

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
