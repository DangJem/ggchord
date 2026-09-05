# ggchord 设计路线图

本文档记录 ggchord 的实际实现状态、已确认的问题、目标 API 以及版本计划。
状态以当前代码为准，不再把“已有设计”标记成“已经完成”。

当前正式版本为 **v0.11.0**（发布提交 `1d85383`）。以正确性和 ggplot2 语法一致性为中心的基础重构
已于 2026-08-31 随 v0.9.0 发布；表达能力、密集 ribbon、静态预览和默认视觉
收尾已于 2026-09-01 随 v0.10.0 发布。

### 2026-09-05 收尾验收

- `fc83464b` 包含 radial 默认布局、auto 混合求解器、有符号曲率、预览适配和图例边距；
  使用前应核对加载目录，Git worktree 的提交不会自动更新主目录或已加载的 R 会话。
- 常规复验入口为 `Rscript tools/validate-release.R <mode> <临时输出目录>`；
  `mode` 为 `geometry`、`labels`、`output`、`benchmark` 或 `examples`。
- `labels` 保留 default/unequal 两种场景、radial/auto 两种模式；auto 还断言左右列
  共用横坐标，而其余标签确实由 radial 求解。输出记录设备、会话、计数、耗时和图像。
- 用户决定本轮停止极端参数压力测试；大幅正负弯曲、极短弧、超密集标签不作为
  v0.11.0 收尾门槛。先前极端渲染中观察到的边界表现仍需后续定位，不算已解决。
- 常规标签四组均保留 32 个标签，0 文字碰撞、0 实际引线交叉；默认 auto
  的 17 个标签使用左右列、15 个使用 radial。实际耗时仅作记录，不设机器相关门槛。
- PNG/SVG 自动预览、PNG 四种单位的显式 11×7 英寸等价尺寸、PDF 导出和
  图例五种位置已复验。过宽图例给出明确警告，显式分行后的图像已检查。
- 本机 R 4.6.0、ggplot2 4.0.3、svglite 2.2.2；历史 SVG 依赖缺失已不适用。
- 1000/5000 ribbon build 约 0.63/1.50 秒；32/80 标签常规基准已运行。
- 完整包检查包含英语/汉语 vignette 重建，初次 0 error / warning，唯一 NOTE
  为远程系统时间不可校验；最终检查关闭该远程校验，其余检查保持启用，
  结果为 `Status: OK`（0 error / warning / note，含全部测试与 vignette 重建）。
- 主目录原先停在 `fb0a398`，未提交内容包含旧 aligned 实现。同步前完整保留
  Git stash 与备份分支；不将旧求解器重新覆盖到已验证的 radial/auto 实现上。
- 修正随包英语/汉语教程中的已删除参数；保留旧说明图片并明确其历史属性。
  完整网站、README 与说明图片重写、六边形 logo 仍留到正式版阶段。

### 2026-09-05 标签、画幅与正式版待办

- 当前默认 `gene_label_layout = "radial"`：每条序列独立偏移轮廓，水平文字，
  法线直线或“短法线段＋长连接段”；拥挤时沿轮廓展开。
- `auto` 替代 `aligned`：左右垂直列，其余位置使用 radial；不再维护 aligned 接口。
- 画幅以标签、图例和几何实际边界为依据，预览自动收紧高度；显式导出尺寸优先。
- 修复曲率的负值截断及 0/1 分支不连续；后续扩展大幅正负弯曲下 ribbon、feature、
  刻度、文字和区域高亮的几何压力矩阵，区分正常的形状相交与数值错误。
- 正式版网站导航参考 [gggenomes](https://thackl.github.io/gggenomes/index.html)：
  Get started、Reference、Articles（悬停展开多页）、Changelog；发布前另行讨论信息架构。
- 正式版统一说明图片来源和生成脚本，建立可复用资产清单；网站、README、手册尽量
  共用相同图片，只为不同介质所需的尺寸或格式生成变体，避免重复生成和提交。
- 正式版设计 R 包六边形 logo；届时讨论图形、配色、字标及 SVG/PNG 导出，本轮不制作。

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
| v0.10.0 | 发表级视觉、feature 几何、显式 ribbon 聚合、布局导出与静态预览 | 已正式发布（2026-09-01） |
| v0.11.0 | 高级轨道、显式多环、ggplot2 grammar 内核与 API 收敛 | 已正式发布（2026-09-05） |
| v0.12.0 | link 几何命名、单点连接与局部平滑避障 | 开发中；已完成规范 ribbon 命名迁移 |
| v0.13.0 | 单序列圆图、限制性酶切位点与注释表达 | 规划中 |
| v0.14.0 | 直线型多序列共线性坐标与同步布局 | 规划中 |
| v0.15.0 | 数据准备适配器、规模化诊断与生态互操作 | 规划中 |
| v0.16.0 | 定量多轨道、无障碍性与输出质量审计 | 规划中 |
| v0.99.0 | API 冻结候选与发布前审计 | 规划中 |
| v1.0.0 | API 冻结、完整文档和长期兼容承诺 | 规划中 |

### 九项已确认方案的路线图审计

| 项目 | 落实位置 | 状态 |
| --- | --- | --- |
| 一、建议的版本路线 | 本节 v0.11.0–v1.0.0 分期 | 已落实，并扩展至 v0.16.0 和 v0.99.0 |
| 二、v0.11.0 收尾 | 第六节 | 已正式发布（2026-09-05） |
| 三、单基因组圆图 | v0.13.0 A | 已规划 |
| 四、限制性酶切位点 | v0.13.0 B | 已规划，包含邻近位点共享主干 |
| 五、直线型多序列共线性 | v0.14.0 | 已规划为专用 coord，不静默忽略冲突参数 |
| 六、基因标签布局与自适应文字 | v0.11.0 F | 已实现紧凑序列轨道、换行/省略和内侧对齐修复 |
| 七、示例数据调整 | v0.11.0 F、v0.13.0 C | 已建立确定性生成规则，后续补单基因组 fixture |
| 八、前期数据准备 | v0.15.0 | 已规划 R 适配器与报告式流程，外部比对不在绘图 build 中运行 |
| 九、与类似工具相比的后续功能空间 | v0.12.0–v0.16.0 | 已规划 link、多轨道、focus/sync、可访问性和输出质量 |

该审计表是后续更新的最小覆盖线；版本可以延后，但不得在没有说明的情况下删除其中任一项。

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

- 旧 scale 类 geom 参数在 v0.9.0–v0.10.0 中作为兼容入口，并已在 v0.11.0 删除；
- 自定义 Stat/Geom 的进一步拆分属于内部演进，不改变 v0.9.0 的公开契约；
- 全量文档、网页和图片重写留到 v1.0.0；
- feature shape、多环、密集 ribbon 和正式布局导出按后续版本推进。

### v0.9.0 正式版验收记录（2026-08-31）

- testthat 保持为 15 个公开行为场景、90 项断言：0 fail、0 warning、0 skip；
- `R CMD check` 包含英语和汉语 vignette 重建：0 error、0 warning；唯一 NOTE 为隔离
  环境无法联网校验系统时间，与包代码无关；
- 用户四序列示例在 `aligned`、`radial`、`arc` 下分别以 6×4 和 12×8 英寸构建，
  均为 0 标签重叠、0 指示线交叉、0 二次裁切损失；
- 相同设备和输入重复构建的标签与线段数据完全一致；
- 80 标签 `aligned` 基准的三次构建中位数约 0.36 秒，低于 2.5 秒验收线；
- PNG、PDF、SVG 均成功导出；默认、minimal、dark、publication 和灰度效果已
  人工检查；临时 SVG 验收工具不加入包依赖；
- 英语和汉语 README 没有因本轮一般功能收尾发生变化。

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

- v0.9.0–v0.10.0 已完成旧接口的过渡和迁移说明；
- v0.11.0 完成第一次集中收敛，直接删除旧的 scale 类 geom 参数，
  不再保留软弃用或静默兼容；
- 删除的参数必须立即给出可执行的迁移方向；
- `data`、`mapping`、`position`、`show.legend` 和 `inherit.aes` 采用标准
  ggplot2 语义；
- 已完成精简的 `geom_gene_label_repel()` 不恢复随机 force/seed 参数；
- v1.0.0 发布后再开始长期兼容承诺。

### 2.5 设计借鉴、来源说明与许可

- 说明文档和网页教程应明确写出某项功能借鉴了哪个软件、包或论文的什么
  思路，包括函数语义、参数组织、布局策略和视觉风格；
- 公开函数和参数仍使用通用且自洽的命名，不把第三方品牌名写进 API；
- 可以直接吸收成熟软件的信息层级、交互逻辑和视觉规则，不必回避说明借鉴关系；
- 代码、图标、图片、配色文件或其他资产只在其许可允许时直接复用，并在
  `LICENSE` / `NOTICE` 或对应页面记录来源、版本和许可；无明确授权的专有资产
  只借鉴其可抽象的设计原则；
- 当前已完成功能的归因矩阵：整体 grammar 参考 ggplot2；环形层级参考 Circos；
  基因箭头与注释可读性参考 clinker 和 DNA Features Viewer；径向/贴弧标签参考
  SnapGene 和 Geneious；focus、方向同步、feature stacking 参考 gggenomes；
  静态尺寸预览的工作流参考 ggview。

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
ribbon_fill
ribbon_alpha
ribbon_colour
ribbon_linetype
gene_fill
feature_fill
feature_shape
region_fill
```

首批公开 scale：

```r
scale_seq_colour_manual()
scale_seq_color_manual()       # 美式拼写别名

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
scale_feature_shape_manual()
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

- sequence、gene 和 feature 的离散色板优先使用色觉友好且灰度可区分的
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
- `seq_colors` 已迁入 scale，旧参数保留迁移警告；
- 序列显示文字由 `geom_seq_label(labels = ...)` 控制，scale labels 独立；
- `linewidth` 和 arrow 已通过不同静态输出设备检查；
- `seq_order`、`seq_orientation`、`seq_gap`、`seq_radius`、
  `seq_curvature` 保留为布局参数。
- v0.10.0 移除序列分组、组间距、组标签和组配色整套接口；这些概念
  与序列顺序、半径和显式多环的职责重叠，且没有足够清晰的通用语义。

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
- feature 的 block、chevron、lollipop 已在 v0.10.0 完成。

#### gene label

- `geom_gene_label()` 是唯一的固定位置和手工微调图层；
- `geom_gene_label_repel()` 现提供 `radial`、`auto`、`arc` 三种确定性模式；
- 所有 label geom 的 `data`、`mapping` 已生效；
- leader 的 colour、linewidth、alpha 默认由 theme 元素控制；
- 固定文本 size 不再创建全局 `scale_size_identity()`；
- v0.10.0 新增水平/径向/切向文字、内外侧以及隐藏/微调/允许重叠策略；
- 逐序列、逐链的 rotation/radial/circumferential offset 继续提供手工自由；
- 不规划 `geom_gene_label_manual()`，避免两个固定标签图层职责重叠。

#### 自动轴和 `geom_seq_label()`

- `geom_axis()` 已删除，序列轴作为默认内部坐标装饰自动绘制；
- axis breaks/labels 已移入 `scale_seq_position_continuous()`；
- 轴线、主刻度、次刻度和文字由 `theme_ggchord()` 分别控制；
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
- 按用户确认保持英语和汉语 README 原有主体内容，不加入一般功能更新；
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

状态：已完成；无参调用在 v0.9.0–v0.10.0 兼容并警告，v0.11.0 已删除。

### C. scale/theme/guide/coord

- 实现第三节列出的首批接口；
- 建立首批色觉友好默认色板、图例和主题元素；v0.10.0 继续视觉校准；
- 旧参数在 v0.9.0–v0.10.0 进入弃用期，并于 v0.11.0 删除；
- 当前阶段只维护必要的 NEWS 和生成文档，README 保持稳定；完整文档重构留到 v1.0.0。

状态：已完成。README 按确认保持稳定，必要变更已记录在 NEWS 和生成的 Rd 中。

### D. 标签布局

- 保留 `radial`、`auto`、`arc`，`aligned` 已由 `auto` 替代；
- 继续保证不同 radius、curvature、gap、orientation、rotation 和设备尺寸；
- 标签实现迁入逐层架构，不恢复随机 force 参数。

状态：已完成；三种模式通过不同序列方向、半径、曲率、间距、旋转和静态设备
组合验收。

---

## 五、v0.10.0 — 表达能力与默认视觉

**状态：已于 2026-09-01 正式发布。**

### A. 固定标签与发表级默认视觉

- 增强 `geom_gene_label()`，使其同时承担简洁默认标签和精确手工微调；
- Identity 色条以 8×6 英寸画布上的 50 mm 长、3.6 mm 厚为参考，
  随实际输出设备等比缩放；缩放设有上下限，避免极端画布中太大或太小；
- sequence/gene 图例符号表达实际方向，并减少过粗的线条和箭头；
- 默认白底、标题、轴线、标签、ribbon 透明度和图例间距统一校准；
- 所有自动生成的图例 key、文字、标题、边距和色条同步缩放；用户显式给出的
  guide 尺寸不被覆盖；
- 文字碰撞框和自动坐标范围按真实设备短边估算，不再将小画布误当成 6 英寸画布；
  默认图例背景透明，不遮挡小画布上合法伸入边距的标签；
- v0.10.0 当时的默认 Viewer 和主要验收图使用 1:1 方形画布；v0.11.0 已改为
  11 英寸宽、内容适配高度。实际设备尺寸仍由 RStudio 或
  `ggsave()` 控制，非方形期刊版式继续作为兼容性场景；
- 用 4×4、6×6、8×8 英寸为主，并以 6×4、8×6、12×8 英寸和 PNG/PDF/SVG
  补充验收。

视觉设计明确借鉴：[Circos](https://genome.cshlp.org/content/19/9/1639)
的环形信息层级与克制 ribbon、[clinker](https://academic.oup.com/bioinformatics/article/37/16/2473/6129045)
的发表级基因箭头、[DNA Features Viewer](https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/)
的注释冲突处理，以及 [SnapGene](https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map)
和 [Geneious](https://manual.geneious.com/en/latest/Sequences.html) 的局部特征标签与拥挤隐藏思路。
实现时可直接吸收成熟的视觉规则；若复用具体代码或资产，按 2.5 节记录许可和来源。

### B. 通用 feature geometry

实现真正独立的 feature shape：

```r
geom_feature(aes(feature_shape = type))
scale_feature_shape_manual(values = c(
  CDS = "arrow", tRNA = "block", repeat_region = "chevron"
))
```

首批只考虑 `arrow`、`block`、`chevron`、`lollipop`。`geom_gene()` 保持为 gene
arrow 的便捷封装，但验收标准是几何和视觉等价，不承诺内部字节完全一致。
四种形状已在 v0.10.0 实现；`lollipop` 的圆头在局部切线—法线
正交坐标中直接生成，不再因序列半径或曲率显示为扁椭圆。

v0.10.0 同时移除 `geom_seq_group_label()`、`scale_group_colour_manual()`
及 `geom_seq()` 中的所有 `seq_group*` 参数。旧参数立即报错，不保留
隐式分组布局。

### C. 密集 ribbon

不在 `geom_ribbon()` 中加入 `ribbon_reduce`，也不让 geom 静默抽样或改变图义。
v0.10.0 提供两个显式数据/布局工具：

```r
bundle_ggchord_ribbons()
optimize_ggchord_layout()
```

前者按序列对、方向、归一化中点网格和用户分组聚合，并保留源行映射；后者在
不改变输入的前提下确定性优化序列顺序和方向。`geom_ribbon()` 继续严格保持
“一行数据对应一条 ribbon”。无法消除的非平面连接用 bundling、透明度和明确
诊断缓解，不承诺所有 ribbon 都能无交叉。

### D. 局部障碍感知 ribbon

- `ribbon_gap = NULL` 是默认自动模式；
- 每个 ribbon 端点按自己的基因组区间判断，不在整条序列上统一留白；
- 只有 `geom_gene()` / `geom_feature()` 的实体多边形算障碍，文字和连线不算；
- 无障碍端点贴近 `geom_seq()`，有障碍端点保留实际所需安全间距；
- 显式设置数值 `ribbon_gap` 时关闭自动判断，严格使用用户数值。

### E. 布局导出

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
该单一正式接口已并入 v0.10.0。

### F. 性能

- sequence 基础布局缓存一次；
- 每层只计算自己的几何，未添加的 ribbon/gene 实体不生成多边形；
- 大型 ribbon 使用 `tools/benchmark-layout.R` 独立验收，不将机器耗时写入 testthat；
- 不使用影响绘图结果的全局可变缓存；
- 不整包导入 ggplot2、grid 或 grDevices；内部调用使用命名空间限定，仅为 S3
  注册保留不可避免的精确导入，减少 IDE 自动补全干扰。
- 选择性重导出构建 chord 图直接需要的 ggplot2 辅助函数，包括 `aes()`、
  `labs()`、`ggtitle()`、`guides()`、`theme()`、主题元素构造器、
  `annotate()`、`ggsave()` 和映射辅助函数；不重导出 Cartesian geom、
  facet、coord、通用 scale 或完整主题预设。

### G. 输出尺寸预览（v0.10.0 历史设计）

新增 `view_ggchord()`：使用标准 `ggsave()` 参数把图临时渲染为 PNG/SVG，再在
IDE Viewer 或浏览器显示；默认 8×6 英寸，`viewer = "none"` 只返回临时文件路径。
它只用于检查真实导出尺寸，不引入 canvas、专用保存函数、Shiny、Plotly 或
HTML widget。正式保存仍使用 `ggsave()`，临时目录仅保留最近 20 次预览。

### v0.10.0 正式版验收记录（2026-09-01）

- 精简后的 testthat 全部通过，`R CMD check` 为 `Status: OK`；
- 1000 ribbon 基准：原始 build 约 0.535 秒、聚合约 0.019 秒、布局优化约
  0.076 秒；5000 ribbon 基准：原始 build 约 1.410 秒、聚合约 0.035 秒、
  近似布局优化约 0.029 秒；
- 1000/5000 ribbon 均保持输入不变、权重和源行可追踪，超过 2000 条时报告
  `approximate = TRUE`；
- 按发布约束未重建 README、vignette、网站或任何说明图片；完整文档和图片
  统一留到 v1.0.0 前重构。

---

## 六、v0.11.0 — 高级轨道与 grammar 收敛（开发中）

v0.11.0 合并原 v0.12.0 计划，作为 v1.0.0 前的集中 API 收敛版本。当前实现包括：

### A. 高级布局

- `stat_ribbon_bundle()` / `stat_ribbon_density()` 已公开
  `after_stat(bundle_n)`、`after_stat(bundle_weight)` 和 `after_stat(density)`；
- `focus_ggchord_data()` 已支持多 locus 同步裁切、`trim` / `drop` 边界策略和
  来源追踪；
- `position_feature_stack()` 已实现确定性最少径向轨道；
- `seq_ring` 和 `scale_seq_ring_manual()` 已实现显式多环，不猜测环数或半径；
- 继续保留 `geom_ribbon_highlight()`，用于安全、可追踪的叠加筛选。

### B. 标准图层与构建

- 每个公开 geom/stat 返回一个标准 `LayerInstance`；
- 公开 `geom_axis()` 已删除，内部 `GeomChordAxis` 由构建流程自动注入；
- 自动轴默认显示，`theme_ggchord(axis = element_blank())` 会同时从绘制范围和
  标签障碍物中移除整套轴；
- `geom_gene_label_repel()` 使用组合 Geom/gTree 绘制标签和指示线；
- `+.ggchord` 不再展开自定义图层列表；
- `prepare_ggchord_plot()` 成为布局、几何、默认标度和 Coord 的唯一准备入口；
- `ggplot_build.ggchord()` 仅在准备后进入标准 ggplot2 build；
- `coord_chord()` 现在返回真正的 `CoordChord` ggproto；
- 删除全局“最近一次布局”，`get_chord_layout(plot)` 必须显式接收 plot。

### C. 参数和标度收敛

所有公开 geom/stat 统一采用 `mapping`、`data`、`position`、`show.legend`、
`inherit.aes` 和 `...`。`data` 可以是数据框或接收该图层默认数据的函数。

| 已删除入口 | v0.11.0 用法 |
| --- | --- |
| `ggchord(title/rotation/panel_margin/show_legend)` | `labs()` / `coord_chord()` / `theme()` |
| `show_legend` | `show.legend` |
| geom 的 `legend_position` / `legend_key_*` | `guides()` + `guide_ggchord_*()` |
| geom 的颜色、顺序、limits、breaks、name | role-specific `scale_*()` |
| `ribbon_*_by` 和 direction 样式参数 | `aes(ribbon_* = ...)` + 对应 scale |
| `ribbon_alpha`、ribbon outline 参数 | `alpha`、`colour`、`linewidth`、`linetype` |
| feature 的 type/category/label 字符串列名 | `aes(feature_type/feature_fill/feature_label = ...)` |
| `geom_seq(seq_labels/seq_colors)` | `geom_seq_label(labels = ...)` / seq scale |
| `geom_axis()` 及其全部参数 | 自动轴；外观和物理间距用 `theme_ggchord(axis.*)`，breaks/labels 用位置 scale |
| `geom_seq_region(regions/region_* style)` | `data`、`fill`、`colour`、`alpha` |

旧参数不再软弃用或静默忽略，而是立即给出迁移错误。用户提供的 scale 永远优先；
只提供映射而没有 scale 时，根据实际数据推断连续或离散默认标度；同一 aesthetic
出现不兼容类型时明确报错。

### D. aesthetic、主题与命名空间

- role aesthetic 使用 quosure/tidy evaluation，在构建中只解析一次并保留原始列；
- 控制坐标的 role aesthetic 禁止 `after_stat()` / `after_scale()`；视觉 aesthetic
  支持标准 staged mapping；
- `theme_ggchord()`、`theme_ggchord_minimal()`、
  `theme_ggchord_publication()`、`theme_ggchord_dark()` 使用完全相同的点号参数；
  `text` 是父级字体入口，轴、序列标签、基因标签和图例字号通过 `rel()` 继承；
- 自动轴由 `axis.*` 控制；轴间距、主次刻度长度和文字偏移使用物理单位，位置
  scale 继续独占 breaks、minor breaks、labels、limits 和 transform；
- `legend.seq.*`、`legend.ribbon.*`、`legend.gene.*`、
  `legend.feature.*`、`legend.region.*` 独立控制五种角色图例，并逐项继承通用
  `legend.*`；role 父级 `element_blank()` 只隐藏该角色；
- 英式 `colour/colours` 为规范拼写，美式 `color/colors`、role aesthetic、scale
  函数和 `guide_ggchord_colorbar()` 为完整等价别名；同次调用不得混用；
- 删除 `base_size`、`base_family`、旧下划线主题参数以及重复的
  `theme_ggchord_elements()`；
- minimal 弱化辅助线，publication 使用紧凑期刊排版，dark 保留完整深色适配；
- `theme()`、`ggtitle()`、`labs()` 等仍是 ggplot2 原函数的精确精选重导出，
  不创建同名包装器，也不导出无关 geom/facet/coord/完整主题。

### E. 接口价值审计结论

- 保留核心 geom、role scale、guide、coord、数据验证/清理/导入、ribbon 整理、
  Viewer、focus、stack、ring 和布局导出；
- `get_chord_layout()` 用于完整内部调试，`export_ggchord_layout()` 用于稳定导出；
- color/colour 双拼写别名遵循 ggplot2 生态，保留；
- 删除全局布局环境、旧 scale 兼容 helper、图层列表展开和重复 build 路径；
- 内部纯几何、文字测量、碰撞检测和验证 helper 保留。

### F. fc83464b 的历史验收记录（2026-09-05）

- 默认布局已改为 `radial`，沿各自序列偏移轮廓排布水平文字；`auto` 将左右
  垂直列与其他方向的 radial 组合；`aligned` 已移除；
- `geom_gene_label_repel()` 增加
  `gene_label_segment_overlap = "fade" | "clip" | "show"` 和
  `gene_label_segment_overlap_alpha`；默认将穿过其他标签的指示线区段
  以低不透明度连续绘制，同时保留硬裁切和完整显示选项；
- `geom_gene_label_repel()` 增加 `gene_label_fit` 和
  `gene_label_max_lines`，对初始碰撞或物理宽度过大的文字执行换行、省略或
  组合策略；显式
  `gene_label_wrap` 仍具有最高优先级；
- `geom_gene_label(gene_label_side = "auto")` 的水平文字按相对自身序列曲线的
  实际位移决定对齐，不再错误地使用相对全图原点的象限；
- 单链 gene 数据的 Strand guide 只生成实际存在的箭头 key，避免 guide 行数不匹配；
- `data/gene_data_example.rda` 改为由未修改的 `examples/gene_track.tsv` 确定性
  生成的紧凑绘图 fixture：每条序列八个分散且长短有层次的可见特征、正负展示链各四个，并保留
  `source_strand` 说明原始注释方向；生成规则位于
  `data-raw/generate_example_data.R`；
- `data/ribbon_data_example.rda` 同样由未修改的 BLAST 原始文件确定性
  生成：去除小于 300 bp 的短区段，密集序列对最多保留 3 条空间分散的
  代表连接，稀疏序列对原样保留；数量从 31 条降为 13 条；
- testthat 保持轻量公开行为测试；
- 主图、三种自动标签、组合 Geom、stat、feature shape、region、highlight、ring、
  data function、标度推断、自动轴、角色图例、双拼写别名、布局导出、标准
  `ggsave()` 和 Viewer 均已通过烟雾测试；
- 上轮 `devtools::check(document = FALSE, manual = FALSE, vignettes = FALSE)`
  达到 `Status: OK`（0 error / warning / note）；222 项测试通过；
- `radial` / `auto` 各四组半径、曲率、方向与旋转矩阵共八组临时渲染，均保留
  32 个标签，文字冲突与实际连线交叉均为 0；极短密集弧仍可能需要整组外移；
- 径向求解增加整组扇形展开和上下额外间距平衡；相同 11×7 英寸示例中，
  顶部平均引线缩短约 18%，上下平均长度差缩小约 62%；
- 自动预览直接导出已测量的 gtable，避免重新布局触发长宽比振荡；1650 像素宽
  默认预览左右外部安全边距约 4 / 12 像素，图例不裁切；
- 默认 Viewer 为 11 英寸宽、按内容自适应高度；默认示例约为 11×7.48 英寸；
  单圆密集标注及宽幅 PNG/PDF 已生成到系统临时目录，不写入说明图片目录；
- 前一轮 `aligned`、`radial`、`arc` 和小画布的临时 PNG 及 PDF 已完成视觉
  验收，不写入仓库；SVG 代码路径保留可选 `svglite` / Cairo 检查，
  当时因未安装 `svglite` 且 X11 动态库缺失，未做本地 SVG 视觉渲染；
  最新收尾环境已具备 `svglite` 并完成 SVG 导出；
- README、vignette、网站及说明图片在开发阶段保持不变。

### G. 灵活参数格式的保留边界

早期“灵活参数格式”仍然受支持，但只应用于确实具有逐序列或逐链语义的几何
参数：

- 序列几何继续接受单值、与序列数相同的无名向量、按 `seq_id` 命名的向量，
  以及现有等价列表形式；典型参数包括 `seq_radius`、`seq_gap`、
  `seq_curvature` 和 `seq_orientation`；
- 固定 gene/feature 几何的 offset、width 和固定标签位置参数继续允许按序列和
  `+`/`-` 链分别指定；
- 不把这种列表语法扩展到 aesthetic、scale、guide 或 theme。数据驱动差异使用
  `aes()` 与 `scale_*()`，非数据外观使用 `theme_ggchord()`，避免同一个参数同时
  承担数据映射和外观覆盖；
- v0.11.0 只补齐一致的校验、缺省值和错误提示，不再增加新的输入容器类型；
  v0.99.0 再决定是否保留“按序号命名的列表”等低频形式。按 `seq_id` 命名的
  向量作为长期推荐写法。

旧教程中已经移除的 `axis_label_size` 等名称不代表当前公开接口；自动轴字号和
物理间距分别由 `theme_ggchord(axis.text = ...)` 与 `axis.*` 单位参数管理。

---

## 七、v0.12.0 — link 几何与局部平滑避障

### 开发进度（2026-09-06）

- 从正式版 `1d85383` 的独立 worktree 开始，开发版本号为 `0.11.0.9000`。
- 已实现 `geom_link_ribbon()` 规范入口；旧 `geom_ribbon()` 明确警告，
  保留显式参数、颜色别名和原有 ribbon 图层、标度、主题及布局导出契约。
- stat、测试、验收脚本及直接相关包内文档已迁移；网站与说明图片不改动。
- 后续实施顺序：单点 link 图层及角色标度/图例主题 → 图层障碍 provider →
  局部平滑剖面及 uniform/none 回退 → 针对性几何验收和完整包检查。
- 单点 link 及下述 B–D 仍是待实现设计；本次命名迁移不改变现有避障行为。
- 命名迁移验收：250 项测试通过；包含双语 vignette 重建的包检查为
  `Status: OK`（0 error / warning / note；未生成 PDF manual，关闭远程系统时间校验）。

### A. 连接图层命名

ggplot2 已导出 `geom_ribbon()`，ggchord 继续使用同名函数会造成搜索路径遮蔽、
文档跳转歧义和 IDE 自动补全干扰。推荐建立以 `link` 为总概念的两个明确图层：

```r
geom_link_ribbon()  # 区间到区间的带状连接
geom_link_line()    # 单点到单点的线状连接
```

- `geom_link_ribbon()` 作为现有 ggchord `geom_ribbon()` 的规范名称；
- `geom_link_line()` 接受 `qaccver/saccver/qpos/spos`，不把零宽区间伪装成 ribbon；
- 线状连接支持直线/曲线、方向箭头、`link_colour`、`link_alpha`、
  `link_linewidth` 和 `link_linetype`；
- ribbon 继续使用 `ribbon_fill` 等带状专属 aesthetic，避免为了改函数名
  而无必要地重命名全部标度和用户映射；
- 新增 `legend.link.*` 只管理线状连接，`legend.ribbon.*` 仍管理带状连接；
- v0.12.0 先新增规范名称并迁移内部实现，旧 `geom_ribbon()` 只保留一个开发周期的
  明确弃用提示；v0.13.0 移除其导出，在 v1.0.0 前彻底消除同名冲突。

不推荐用单一 `geom_link(type = ...)` 同时承担线和带：两者的必需数据列、
图例 key、fill 语义和避障边界都不同，分开函数更利于自动补全和错误提示。

### B. 当前 ribbon 避障范围

当前 `ribbon_gap = NULL` 的自动避障具体规则为：

- 只收集 `geom_gene()` 和 `geom_feature()` 实体多边形所在的基因组区间；
- 只把位于 ribbon 一侧、且与 query/subject 端点区间重叠的实体作为障碍；
- 文字、指示线、轴、序列标签、region 和 highlight 不算障碍；
- query 和 subject 两端独立计算；任一端区间命中障碍时，该端整条前沿
  使用最大所需间距，未命中时使用贴近序列的小间距。

这个实现安全但较生硬，而且新增几何必须修改中央构建函数才能参与避障。

### C. 可扩展障碍协议

将障碍物改为图层自主注册的几何契约，内部至少统一输出：

```text
seq_id, start, end, side, normal_min, normal_max, priority, source_layer
```

- gene、feature、未来的多轨道实体和用户扩展 geom 可以实现同一 obstacle provider；
- 文字和连线默认不注册，透明 region 默认也不注册；
- 协议保存局部法向范围和真实 side，不再只用 `gene_offset + gene_width / 2`
  近似所有 shape；
- v0.12.0 先作为内部契约，等第三方 geom 需求稳定后再决定是否公开。

### D. 局部平滑避障

默认不再将整条 ribbon 前沿统一后移。新算法沿 query/subject 区间采样距离剖面：

1. 无障碍位置使用贴近序列的 `close_gap`；
2. 命中障碍的局部使用该 shape 真实法向范围加安全间距；
3. 障碍边界两侧用余弦或保形三次曲线过渡，并对必需间距做上包络，
   防止平滑后反而穿过障碍；
4. 根据区间长度、障碍边界和曲率自适应采样，检查前沿自交和方向翻转；
5. query 和 subject 分别生成剖面，再与中部 Bezier 连接。

预计提供 `link_avoid = "smooth" | "uniform" | "none"`：`"smooth"` 作为新默认，
`"uniform"` 保留当前保守形式，`"none"` 不考虑实体障碍。显式数值 gap 始终优先，
不被自动剖面覆盖。参数最终名称在原型验收后冻结。

---

## 八、v0.13.0 — 单序列圆图与生物学注释

### A. 单序列圆图

单条完整基因组不应只是“四序列 chord 少三条序列”的退化情况。v0.13.0 计划
提供明确的单序列布局入口，共享现有 gene/feature/region、scale、theme、guide
和标签组件，但不生成无意义的跨序列 ribbon。

- 支持闭合圆形和带缺口的线性化圆图；
- 基因、调控元件、重复区和用户 region 可放在显式内外轨道；
- 标签可选择贴弧、局部径向或外部 callout，并复用碰撞、换行和省略策略；
- 原点、方向、旋转和坐标标签由 coord/position scale 管理；
- 明确参考 SnapGene、Geneious 和 DNA Features Viewer 的信息层级与可读性，
  API 保持通用命名，资产复用按 2.5 节的许可规则执行。

正式命名需要在实现前做最小原型比较；不把单序列行为偷偷塞入
`coord_chord()` 的条件分支。

### B. 限制性酶切位点

计划拆成数据计算和几何表达两层：

```r
sites <- find_restriction_sites(sequence, enzymes = c("EcoRI", "BamHI"))

ggchord(seq_data, gene_data = gene_data) +
  geom_seq() +
  geom_restriction_site(data = sites)
```

- `find_restriction_sites()` 返回酶名、识别序列、切割位置、黏性/平末端及来源；
- 酶数据库版本必须可追踪，用户也可传入自定义 motif；
- `geom_restriction_site()` 只负责刻线、标签和指示线；
- 邻近标签采用确定性的共享主干（trunk）后分叉，不把不同切点合并成一个数据点；
- 可按酶、切割次数和窗口过滤，重复名称仍保留全部位点；
- 序列搜索优先使用轻量实现；大型序列可选用 Biostrings，但不设为强制依赖。

### C. 示例体系

- `examples/` 保留原始 FASTA、GFF3、BLAST 和 TSV，不为视觉效果修改源记录；
- `data/` 只保存小型、清晰、可重复生成的教学 fixture；
- `data-raw/` 保存生成规则，明确任何抽样、过滤或展示链重编码；
- 补充单基因组、无 ribbon、单链 gene、密集标签和酶切位点的最小 fixture；
- 开发阶段不生成文档图片，统一在 v1.0.0 文档重构时出图。

---

## 九、v0.14.0 — 直线型多序列共线性

新增专用 coord 是合理方向，但不能“静默忽略”弦图参数。计划先验证
`coord_collinear()` 原型：

```r
ggchord(seq_data, ribbon_data, gene_data) +
  geom_seq() +
  geom_link_ribbon() +
  geom_gene() +
  coord_collinear()
```

- 同一套 sequence/gene/feature/ribbon 数据和 role aesthetics 可在弦图与直线图间切换；
- sequence 变为水平轨道，ribbon 变为相邻或跨轨道连接，gene/feature 使用法线偏移；
- `seq_order`、`seq_orientation`、position scale 和 `focus_ggchord_data()` 保留
  对应语义；
- `seq_radius`、`seq_curvature`、`seq_ring` 等弦图专属参数若被显式设置，立即
  给出可执行的冲突信息，不能静默忽略；
- 可无损转换的参数才转换，例如 `seq_gap` 对应轨道间距；
- 先支持全局坐标和相邻比较，再评估类似 gggenomes `sync()` 的局部锚定、
  feature stack 和局部翻转；
- 直线模式共享 scale/theme/guide，不另建一套重复 geom API。

该版本重点是坐标契约和语义一致性，不承诺在第一版实现自由 track 容器或任意
facet 组合。

---

## 十、v0.15.0 — 数据准备与生态互操作

ggchord 应降低“已有结果转成图”的门槛，但不把 BLAST/DIAMOND 等外部程序本身
打包进绘图包：

- 保留并增强 `read_fasta_lengths()`、`read_gff3()`、`read_blast()`；
- 增加 GenBank/GBFF、BED、PAF 和常见 synteny 表的可选适配器；
- 提供统一列名映射与 `as_ggchord_*()` 转换器，接受 Biostrings、GenomicRanges
  或常见 R 包产物时不强制依赖这些包；
- 增加 `prepare_ggchord_data()` 报告式工作流，串联导入、验证、过滤、去重、
  merge、bundle 和布局建议，但每一步仍可独立调用；
- 对 BLASTN/BLASTP、DIAMOND、MMseqs2 等只提供输入格式说明、可执行文件检查和
  结果导入，不在核心包中隐式安装、下载或运行外部二进制；
- 对小型纯 R 比对可评估可选适配器，但必须明确其适用规模，不能让绘图函数在
  build 阶段启动序列比对；
- 增加复杂 ribbon 的诊断摘要，例如连接密度、非平面交叉下界、推荐 bundling
  bins 和适合的 order/orientation，而不自动改变图义。

生态比较重点参考 gggenomes 的 tidy 多表输入、focus/sync 和 feature stacking，
参考 clinker、genoPlotR、gggenes 的共线性表达，参考 circlize 的多轨道扩展性，
并继续保持 ggchord 的优势：同一 grammar 下的弦图、显式预处理、可导出布局和
角色专用 scale/theme/guide。

---

## 十一、v0.16.0 — 定量多轨道与输出质量

在 link、单序列和直线坐标稳定后，再补齐类似工具中对科研图常用但 ggchord
尚缺少的定量轨道：

- 设计统一的显式 track/ring 契约，不根据列名自动猜测轨道；
- 评估 coverage、GC content/GC skew、分类 tile、折线和柱形轨道；
- 数值映射继续使用 role-specific scale，轨道几何可向障碍协议注册实体范围；
- 提供色觉友好、灰度和最小字号/线宽诊断，但不强制修改用户风格；
- 对 PNG、PDF、SVG 的实际输出尺寸、字体、透明度、裁切和 guide 比例生成诊断报告；
- 交互式绘图仍留到 v1.0.0 后，v0.16.0 不恢复 Plotly 或引入 Shiny。

---

## 十二、v0.99.0 与 v1.0.0 — 稳定 API

v0.99.0 只做发布候选审计，不再扩张功能：

- 清点全部导出函数、参数、aesthetic、S3 类和返回对象；
- 删除无调用或职责重复的接口，冻结保留接口；
- 对不同 R、ggplot2 版本和主流平台执行兼容检查；
- 建立少量稳定视觉基准和性能基准；
- 修复 release blocker，不再加入新布局模式。

v1.0.0 完成发布与文档：

- 删除已完成迁移的旧参数；
- 冻结核心 aes、scale、theme、guide、coord 和 layout export 契约；
- 基于最终 v1.0.0 API 重写全部英语和现代汉语说明文档、README、vignette、
  roxygen、Rd、网页手册、FAQ 和迁移指南；面向用户的语种名称统一写作“英语”
  和“汉语”，不使用“中文”作为语种标签；
- 网站每一个页面都建立英语/汉语一对一对应路由和显眼的语种切换，页面层级、
  代码示例、参数和锚点尽量对齐，CI 检查缺失翻译或孤立页；
- 英语和汉语页共用同一份图片资产和代码生成源，不为翻译版重复保存图片；
- 网站视觉重构同时优化字体层级、行宽、对比度、代码折叠、页内目录、
  搜索、移动端和键盘可访问性，以阅读效率而不是装饰复杂度为验收标准；
- 删除旧的示例、旧图片及过时表述，从最终 API 重新编写最小且连贯的示例；
- 重新生成所有带图片输出的示例，清理不再被引用的 `man/figures` 和网页图片；
- 文档结构遵循“安装 → 数据格式 → 最小出图 → 常用调整 → 进阶功能 → FAQ”，
  让新用户能够快速上手；
- README 继续保持稳定简洁，把详细功能说明留给 vignette 和 reference；
- 每项主要布局、函数簇和默认风格在英语/汉语文档中都列出设计参考对象、
  借鉴内容和可用的来源链接，资产许可遵循 2.5 节；
- CRAN check、最小视觉快照和跨平台检查；
- 明确支持的 R 和 ggplot2 版本；
- 无严重已知正确性问题后才发布 1.0.0。

v0.10.0 的 `view_ggchord()` 只是静态导出尺寸预览，不属于交互式绘图。真正的
交互方案不属于 v0.9.0–v1.0.0 的发布范围；v1.0.0 发布后再基于稳定的 scale、
guide、coord 和 layout export 契约单独评估，不在本路线图中预先承诺接口。

---

## 十三、testthat 精简策略

测试只保留公开 API 的最小行为，不再通过大量精确坐标断言锁定内部实现。

### 保留

- `ggchord()` 创建对象并能完成一次基础 build；
- 核心 geom 组合能 build；
- 三种 label layout 能 build；
- validate、clean、三个 import、三个 ribbon utility 的单一成功路径；
- region、highlight、feature 的单一成功路径；
- 每个已确认严重 bug 修复后保留一个最小回归用例；
- 一个最基本的错误输入测试。

### 删除

- 单独锁定无公开行为意义的内部函数；已确认几何缺陷可使用内部检测器
  验证公开布局结果，避免用精确坐标替代碰撞、方向和数据保留断言；
- 精确浮点坐标、文字框尺寸和迭代次数；
- 默认 testthat 中多个设备尺寸的重复组合（完整矩阵放在独立验收脚本）；
- guide ggproto 内部字段；
- 同一参数的多种等价输入格式穷举；
- 重复构建、序列化、文件输出等低价值集成测试；
- 易受机器性能影响的固定秒数测试；
- 与实现细节耦合的内部 scale 标记测试。

测试目标是快速发现“公开函数完全不可用”，而不是证明所有排列组合。复杂布局的
质量主要通过少量视觉示例、人工验收和专门 benchmark 检查，不塞入默认 testthat。

---

## 十四、验收规则

- 精简后的 `testthat::test_local()` 必须通过；
- `R CMD check` 必须通过；
- 新的公开参数至少有一个可运行文档示例；
- 破坏性变化必须出现在 NEWS 和迁移表；
- 默认示例在静态图中的外观不发生无理由变化；
- 多个同类图层、正反 orientation、不同 radius/curvature/gap 是核心人工验收场景；
- v1.0.0 前的开发文档避免无必要的反复改写。
