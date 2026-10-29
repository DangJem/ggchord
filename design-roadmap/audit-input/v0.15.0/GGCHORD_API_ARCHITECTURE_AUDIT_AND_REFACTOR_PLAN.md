# ggchord API 全面审计与分阶段重构执行规范

> **文档定位**：交付 Codex 的架构审计、决策与实施指导文件；同时是开发者评审用设计规范。
>
> **审计基线**：`DangJem/ggchord`，`main`，`bdbb3326847fcbb71928eed41647c209bd545871`，`DESCRIPTION Version: 0.15.0`；核对日期：2026-10-30。
>
> **总体结论**：三种主要绘图类别（弦图、单序列圆形图、线性共线性图）已建立，不在本阶段增加第四种大类。**从 v0.16.0 即开始分阶段重构**，可以在重构中补齐确有价值的新功能/API。**v1.0.0 是架构、公开 API、教程、网站、示例、验收均已稳定后的正式发布版本，不是开始重构的版本。**
>
> **证据等级**：`C`=已从当前源码明确确认；`H`=高价值设计假设、需要实验/调用证据；`D`=待项目所有者决策。严禁把 H/D 自动当作 bug 或授权修改。
>
> **禁止误用**：本文件不是授权 Codex 立即删除接口、运行完整发布验收、推送/发版或直接修改 `main`。先核查 repo `AGENTS.md`、当前路线图与具体版本文档，记录用户要求与既有文档冲突，征求取舍；只在用户明确要求实施并确认范围后分批修改。

## 0. Codex 必须先理解的目标和边界

用户并不只是希望清理旧函数名，而是希望针对**所有函数、参数、aesthetic、图层、Position、Stat、Scale、Guide、Coord、Theme、导入与数据处理 API、S3 类与返回对象，以及所有内部布局函数和状态**，重新审视“职责、可发现性、功能冗余、参数默认、错误处理、公开尺度、语义准确、组合性和扩展能力”。项目大部分代码由 AI 辅助完成，因此不可假定现存设计是有意为之的最优解；也不可以仅凭字面名称相近就合并生物学含义不同的能力。

**优先顺序**：生物学与坐标正确性 > 行为确定性/可追溯性 > 清晰一致的语义和用户体验 > ggplot2 生态兼容 > 组合性/扩展性 > 性能 > 接口数目最少。参数数量本身不是缺陷，公开行为无充分理由的复杂度才是问题。

**审计与实施分离**：每一类 API 在动手前必须产生“现状、可验证问题、目标设计、替代方案、迁移、验收”六项记录。禁止一次性大改全部模块；禁止仅因未见引用就删除 API；禁止只在 README 上改名而没有修改源码/帮助/测试/示例。

## 1. 源码事实与范围

- 当前 `main` 的 `NAMESPACE`：**130** 个 `export()`；其中 **22** 个是来自 `ggplot2`/`grid` 的重导出，剩余 **108** 个是 ggchord 自有导出名称（包含别名，并不等于 108 个独立函数）；另有 **9** 个 S3 注册方法。
- 当前源码树中有 **83** 个 `R/*.R`，合计约 **1,322,681 字节**；`tests/testthat/` 有 10 个 `test-*.R` 文件（仅为当前基线，后续以实时清点为准）。
- 当前旧路线图 `design-roadmap/DESIGN-ROADMAP.md` 把原 v0.16.0 范围并入 v0.15.0，随后直接安排 v0.99.0。**这已不符合用户最新决策**：需要重新规划从 v0.16.0 开始的多版本重构；保留原 v0.16.0 合并历史，不可静默覆盖历史记录。
- `DESCRIPTION` 的标题仍为 “Multi-Sequence Alignment Chord Diagram Visualization Tool”，虽然描述已明确涵盖 chord、collinear、single circular maps；公开定位与包名文案需要统一，包名 `ggchord` 本身可以保留。
- 已有规范：`AGENTS.md` 规定主分支不能直接开发，版本分支从 `main` 创建；普通开发只做相关轻量检查，完整 `R CMD check`/设备矩阵/发布仅在明确发布授权后执行。最新用户决策与旧路线图冲突应先提出明确修订方案。

**参考文件**：`DESCRIPTION`、`NAMESPACE`、`AGENTS.md`、`R/` 全部源码、`tests/testthat`、`NEWS.md`、`design-roadmap/DESIGN-ROADMAP.md`、`design-roadmap/versions/v0.99.0.md`、`design-roadmap/versions/v1.0.0.md`。历史网站/教程/README 不作为当前 API 的权威证据。

## 2. 全项目设计原则：先统一语义，再统一名称

### 2.1 用户心智模型

ggchord 是一个 **ggplot2 风格的基因组比较与注释可视化框架**，其三种坐标模式共享统一的数据与图层语义；这不是三个彼此无关的图形工具。API 应使用户能够按自然问题找到入口：

1. 我有什么数据？序列 / 区间注释 / 两端连接 / 定量观测。
2. 我要画什么？骨架 / feature 或 gene / 连接 / 标签 / 定量轨道 / 比例尺与轴。
3. 应怎么排列？`coord_*` 选择全局布局，`position_*` 选择局部摆放，`scale_*` 管理数据映射，`theme_*` 管理非数据样式。
4. 要选择哪些数据？显式过滤、聚焦、分段或转换，不通过绘图函数悄悄改变生物学关系。
5. 要交付什么？可复现的静态预览、保存、布局导出与诊断。

### 2.2 必须明确的生物学契约

- **原始坐标**：必须记录 1-based inclusive 与外部格式 0-based half-open 的转换规则；核查单碱基区间、零宽区间、越界、环形跨原点、反向端点及分段窗口；不制造不存在的序列长度或片段。
- **`strand`、`orientation`、`direction`**：链方向（生物学事实）、显示方向（布局选择）、两段连接方向（相对关系）必须严格区分；任何 `auto` 都不能被解释为已证实的生物学方向。
- **`accver` 与 `qaccver`/`saccver`**：ID 严格原样保留；用户提供的 query/subject 顺序不应被图层或适配器擅自反转。
- **gene vs feature**：gene 是具有明确方向性的生物学区间的一种角色；feature 还可包含未知链、复合区间及其他注释。决定是否共用公开 geom 必须依据语义和人机交互，而非仅几何是否都是箭头。
- **link vs alignment**：`GInteractions`、断点连线、同源关系不必具有 BLAST 身份率或 alignment length；保持 link 的广义语义，任何 ribbon 统计都必须说明额外字段需求。
- **track**：数据值的数值域与实际绘图带宽是不同对象；物理尺寸以输出设备为准，不能把 `mm` 传给数据 scale 或用 scale 模拟轨道位置。
- **信息不丢失**：转换、清洗、折叠、聚合、分段、聚焦和导出要保留源行、原 ID、原坐标和变更报告；默认不做不可逆的数据变换。

### 2.3 API 命名与分层规则

- `ggchord()`：创建绘图与承载基础数据，不承接各图层专属设置；构造器 `...` 不应成为隐藏的参数通道。
- `geom_*()`：一种可识别的视觉物件；优先 `mapping`、`data`、`position`、`show.legend`、`inherit.aes` 和真正的固定外观参数，避免塞入导入/聚合/设备碰撞参数。
- `stat_*()`：数据计算、统计汇总；与原始值的统计字段需有可描述、可验证的含义。
- `position_*()`：序列局部空间内位置与堆叠；保留层级差异，不把所有 Position 化为一个“大模式枚举”。
- `coord_*()`：整个画布的坐标模式、方向与区域结构；避免因某个 geom 的例外需求向坐标函数追加特例参数。
- `scale_*()`：映射值域、调色盘、变换与图例；不能控制 track 实际物理位置。
- `guide_*()`：图例的显示与组织；自定义 guide 应严格遵守 ggplot2 扩展契约。
- `theme_*()`：非数据元素与预设；专用主题元素只有在标准 ggplot2 主题无法表达时添加。
- `read_*()` / `as_ggchord_*()` / `validate_*()` / `clean_*()` / `project_*()`：导入、结构转换、检查、修复和坐标投影的职责分别清晰；不隐式串连。
- **公开/内部边界**：公共 API 能被用户稳定发现并有文档；内核的 registry、缓存、绘图中间态、标签候选、私有图元尽量不导出；返回结构只承诺明确定义的稳定字段。
- 参数语法：明确单值、向量、按 `accver` 命名的向量、按链命名的向量、列表与 R `unit()` 的适用层级，不要每个 geom 定义不同“猜测规则”。
- `NULL`、`NA`、缺参、`"auto"`、`"none"`：**必须有统一政策**，分别代表自动、缺测、未提供、自动推断、禁用时，文档必须列出区别。
- 参数顺序：常用参数在前，进阶参数在后；数量过多时考虑对象/独立图层，但禁止为了“签名变短”简单将所有参数藏入 `...`。

## 3. 已定位的问题与审计假设（第一批 issue ledger）

> **C=已在 v0.15.0 代码中定位；H=待完整调用/测试验证；D=待用户设计取舍。** 具体行号可能随提交变化，Codex 应以固定基线重新核对。

| 编号 | 级别 | 证据 | 对象 | 发现和处理要求 |
|---|---|---|---|---|
| F01 | P1/C | `R/clean.R` | `clean_ggchord_data()` | `ribbon_out <- if (is.data.frame(ribbon_data)) ribbon_data else NULL`，`gene_out` 同理；非 NULL 但不合类型的输入可能被静默当作未提供。改成显式检查并给出可行动错误；验证 tibble/GRanges/NULL 的预期行为。 |
| F02 | P1/C | `R/clean.R` | `reversed_interval="sort"` | 报告 `original_value`/`new_value` 在交换端点之前填入，可能无法区分变更前后。交换后再写新值，并验证方向记录。 |
| F03 | P1/H | `R/data-prepare.R`,`R/clean.R` | 清洗默认值 | 独立清洗默认排序反向端点，prepare 中默认为 `keep`；差异有正确性原因，但 UI 不一致，应给出公开契约和是否调整默认值的决策。 |
| F04 | P1/C | `R/layout.R` | `compute_chord_layout()` | 巨型签名并通过 `eval(ggchord_layout_*_step, envir=context)` 使用共享词法环境。输入/输出和阶段依赖隐式，重构为显式上下文与阶段函数前须画出状态依赖图；不能贸然改数值/绘图。 |
| F05 | P1/H | `R/ggchord-build.R` | `compute_chord_geometry_single()` | 按图层 `ggchord_params$type` 汇集状态，旧式局部变量与新图层可能交叉耦合；审查同类多层、图层顺序、clone 与缓存所有权。 |
| F06 | P2/D | `R/theme.R` | `theme_ggchord()` | 庞大显式参数集合；审查哪些属于 ggplot2 基础元素、真正自定义元素、角色图例，并制定清晰分组/文档；不简单拆分成多个无法组合的 theme。 |
| F07 | P2/D | `R/zzz.R`,`R/layout-export.R` | `get_chord_layout()` / `export_ggchord_layout()` | 原始内部布局列表与受控表格导出的公共边界混合。建议优先把受控 export 作为稳定接口，评估旧 getter 是否转为内部调试入口、迁移助手或明确不稳定接口。 |
| F08 | P2/C | `R/geom-feature.R` | `geom_feature_plasmid()` | 已是 `geom_feature()` 的预设薄封装，区别在 position 和默认外观。可保留便利入口，但保证只有一个语义实现与明确文档，不要独立复制绘图逻辑。 |
| F09 | P2/D | `R/geom-gene.R`,`R/geom-feature.R` | gene / feature | 检查哪些几何、标签、审美角色确需分离；可以共享内部实现，公开入口是否保留须通过代表性使用案例判断。 |
| F10 | P2/C | `R/geom-gene.R`,`R/geom-feature.R` | `gene_offset`、`feature_offset` | 已有弃用警告，Position 是新机制。建立迁移表与删除阶段，并检查旧示例是否还引用。 |
| F11 | P2/C | `R/geom-feature.R` | `shape` / `feature_shape` | 同一形状固定设置有双参数入口，且同时提供时报错；决定唯一首选，别名弃用后删除。区分固定形状参数与映射 `feature_shape` aesthetic。 |
| F12 | P2/H | `R/geom-feature-label.R`,`R/geom-gene.R` | label 命名 | `gene_label_*` 与 `label_*` 参数有语义相近却命名体系不同的部分，仔细比较默认值/返回属性，制定命名统一与迁移策略。 |
| F13 | P2/H | `R/geom-ribbon.R`,`R/geom-link-line.R` | link / ribbon 几何参数 | `ribbon_gap`/`link_gap`、`ribbon_ctrl_point`/`link_ctrl_point` 与共同 `link_avoid/branch` 需定义统一意义；不因名字相近合并几何入口。 |
| F14 | P2/D | `R/scales.R`,`R/link-scales.R`,`R/guides.R` | color / colour | 存在别名；推荐保留 ggplot2 规范 `colour` 为标准，`color` 的保留/删除须测用户负担及迁移。 |
| F15 | P2/D | `R/reexports.R`,`NAMESPACE` | 22 个重导出 | 逐一判断便利价值、IDE 自动补全污染和外部依赖风险。`PositionIdentity`、`last_plot()` 等优先审查，不代表可立即删。 |
| F16 | P2/H | `R/data-import.R`,`R/data-formats.R`,`R/data-genbank.R` | `file`/`files` | 统一多文件的输入规范，兼顾路径向量、glob、来源追溯，避免多对同义参数。 |
| F17 | P2/H | `R/data-genbank.R` | `read_genbank()` | 当前生成/调用 Python/Biopython 解析脚本，审查用户环境依赖、外部进程风险、离线场景、字符编码、严格位置/反向互补契约及失败信息；不是说必须删除该方法。 |
| F18 | P2/H | `R/geom-track.R`,`R/layout-track.R` | track | 检查多 track 同 ID 声明冲突、`bounds`、`closed`、`baseline` 的坐标语义、不同设备的带宽、多个同 ID 图层及 scale 独立性。 |
| F19 | P2/H | `R/coord-chord.R`,`R/coord-circular.R`,`R/coord-collinear.R` | coord 共有参数 | `ratio`/`fit`/`expand`/`clip`/`xlim`/`ylim` 的共同定义与不同行为差异需要一张参数矩阵。 |
| F20 | P2/H | `R/data-convert.R`,`R/data-relations.R` | Bioconductor 适配 | `GRanges`、`Seqinfo`、`GInteractions` 的元数据、NA strand、seqlevels、坐标端点保真；避免强依赖和隐式生物学推断。 |
| F21 | P2/H | `R/layout-circular-registry.R`,`R/layout-external-annotation.R`,`R/helpers-label-routing.R` | annotation registry | 比较圆图/线性图标签冲突算法共享程度，评估 registry 的内部稳定性；不公开候选定位中间态。 |
| F22 | P2/H | `R/layout-export.R`,`R/diagnostics.R`,`R/viewer.R` | 导出/预览/诊断 | 分离绘图 geometry export、视觉尺寸预览、诊断信息；统一单位、稳定字段、来源追踪与 cache 的生命期。 |
| F23 | P3/D | `R/ggchord.R`,`DESCRIPTION` | 包的定位 | `ggchord` 软件包名称可保留；DESCRIPTION Title、网站首页、README 需由“仅弦图”改为涵盖三类模式，避免过度包装成万能组学可视化。 |
| F24 | P3/H | `R/zzz.R` 等 | 迁移错误分支 | 大量拒绝过时参数的逻辑需要统一清单，审查是否在稳定版中继续保留无谓兼容检查或提供错误提示。 |
| F25 | P2/H | 全部 `R/` | 内部未用函数 | 必须由源码调用图 + 动态反射/注册检索 + ggproto 回调 + `quote()` 间接引用共同确认；绝不凭 `grep` 零引用删除。 |
| F26 | P2/H | `R/ggchord-frame.R`,`R/render-context.R` | 缓存/设备状态 | `+` 组合后、重复绘制、复制图、改变设备尺寸、多图并存时必须独立且可复现；任何缓存都属于图对象，不泄漏成全局状态。 |
| F27 | P2/H | `R/stat-ribbon.R` | 统计图层 | 检查 `after_stat()` staged mapping 真正遵循 ggplot2 契约，不要把已聚合数据误当原始行坐标。 |
| F28 | P2/H | 多个源文件 | `...` 透传 | 逐个检查 `...` 是向 ggplot2 参数合法透传，还是掩盖未知参数、拼写错误与弃用行为；定义严格错误策略。 |
| F29 | P2/H | `R/geom-seq-scale.R`,`R/scales.R` | 名称可发现性 | `geom_seq_scale()` 是独立比例尺图元，`scale_seq_position_continuous()` 是坐标映射，用户可能混淆；需要命名或文档中的显著区别。 |
| F30 | P2/H | `R/common-features.R`,`R/restriction-sites.R`,`R/primers.R` | 内置序列识别算法 | 明确搜索算法和静态绘图的界限、许可/数据库来源、结果置信度、性能、图形输入格式；避免用户误当实验验证。 |

### 3.1 如何认定“无用 API”

必须同时提供：`NAMESPACE` 是否公开、包内直接及间接引用、帮助和示例、测试、下游潜在使用、动态派发/注册/ggproto/`quote()` 中的使用、替代者是否覆盖全部行为、删除是否影响既有文件与序列化对象。分类为：**确定无用 / 仅内部 / 仅别名 / 可被新接口完全替换 / 暂不能删除**。对于未调查过的函数不得使用“无用”定性。

### 3.2 如何认定冗余和过度抽象

比较 **生物学输入域 + 输出类型 + 错误行为 + 默认值 + 所属图层职责 + 边界情况**，而不是只比较函数名称或实现代码重复。可以保留“面向用户的薄封装”而共享底层引擎；反之，同一个函数控制不同领域且参数互斥繁多时应考虑拆分。

## 4. API 全量盘点：按实际导出入口分组

> 以下为审计的完整**导出名称分组目录**，并非每个名称都被判定应保留。`api_export_baseline.csv` 提供逐行清单；内部函数/所有参数用配套 `audit_ggchord_api.R` 在 Codex checkout 上生成，不要据此省略人工审查。

### A. 核心与坐标

`ggchord`；`coord_chord`、`coord_circular`、`coord_collinear`。

- 重审构造器 `seq_data`、`ribbon_data`、`gene_data` 的核心数据分工；未来 `features` 与 `links` 是否仅靠图层传入，是否需要新的显式“多表数据容器”，必须用最小真实案例决定。
- Coord 只描述全局坐标；`coord_circular(direction,gap,rotation)` 与 `coord_collinear(align,anchor,orientation,min_support,segments,segment_gap)` 所有参数需逐个解释、检查默认值与互斥关系。
- `coord_chord` 的弦布局与 `coord_circular` 的单序列圆图不是同一生物学图形，只因“都是圆形”不应合并。

### B. 图层与图形对象（21）

`geom_seq`、`geom_seq_label`、`geom_seq_region`、`geom_seq_break`、`geom_seq_center_label`、`geom_seq_scale`；`geom_gene`、`geom_gene_label`、`geom_gene_label_repel`；`geom_feature`、`geom_feature_plasmid`、`geom_feature_label`、`geom_feature_label_repel`；`geom_link_line`、`geom_link_ribbon`；`geom_primer`、`geom_primer_label_repel`、`geom_restriction_site`；`geom_track_tile`、`geom_track_line`、`geom_track_bar`。

- 对每个 geom 列出支持的三个 coord、必需数据列、缺省美学、必需映射、固定参数、Position、图例、空数据、重复/覆盖、输出组件，以及是否 `ggplot2::layer` 原生语义。
- 细分：序列骨架/显示辅助、基因与 feature、link、实验注释、数值轨道。不要把比例尺 `geom_seq_scale()` 错归 scale。
- 检查所有 `geom_*_label[_repel]` 是否可以通过统一参数命名减少用户认知负担，保留必要的视觉语义差别。

### C. 布局/位置（自有 4 + 重导出 1）

`position_strand`、`position_plasmid`、`position_feature_stack`、`position_track`；重导出 `position_identity`。

- 逐项核查 `offset`、`spacing`、`side`、`base_position`、`track`、`order`、`width`、`gap` 在三坐标中的含义、单位及命名格式。
- Position 类对象支持字符串快捷写法是否长期必要；审查 `position = "identity"` 与 `ggplot2::position_identity()` 的一致性。

### D. 数据映射/比例尺（33）

`scale_feature_fill_manual`、`scale_feature_fill_plasmid`、`scale_feature_shape_manual`、`scale_feature_shape_plasmid`、`scale_gene_fill_manual`；

`scale_link_alpha_continuous`、`scale_link_alpha_manual`、`scale_link_color_gradientn`、`scale_link_color_manual`、`scale_link_colour_gradientn`、`scale_link_colour_manual`、`scale_link_linetype_manual`、`scale_link_linewidth_continuous`、`scale_link_linewidth_manual`；

`scale_region_fill_manual`；

`scale_ribbon_alpha_continuous`、`scale_ribbon_alpha_manual`、`scale_ribbon_color_manual`、`scale_ribbon_colour_manual`、`scale_ribbon_fill_gradientn`、`scale_ribbon_fill_identity`、`scale_ribbon_fill_manual`、`scale_ribbon_fill_stepsn`、`scale_ribbon_linetype_manual`；

`scale_seq_color_manual`、`scale_seq_colour_manual`、`scale_seq_position_continuous`、`scale_seq_ring_manual`；

`scale_track_colour_gradientn`、`scale_track_colour_manual`、`scale_track_fill_gradientn`、`scale_track_fill_manual`、`scale_track_value_continuous`。

- 审查所有 `*_color_*` 等于 `*_colour_*` 的别名与其使用率/兼容价值；`link` vs `ribbon` 的 aesthetic 命名规则须形成角色矩阵。
- 解释每个 scale 可操作的 aesthetic、是否连续/离散、默认 guide、是否适用于三坐标；核查 `scale_seq_ring_manual()` 的数据映射与局部空间界限。

### E. 图例/主题/统计

Guide：`guide_ggchord_legend`、`guide_ggchord_colourbar`、`guide_ggchord_colorbar`。

Theme：`theme_ggchord`、`theme_ggchord_dark`、`theme_ggchord_minimal`、`theme_ggchord_plasmid`、`theme_ggchord_publication`。

Stat：`stat_ribbon_bundle`、`stat_ribbon_density`。

- 每个主题的独立参数应与基础 `theme_ggchord()` 的参数签名/行为一致，不应通过运行时动态替换函数体产生令人困惑的帮助接口；保留现有预设的视觉稳定性。
- Stat 不应接收属于 geom 的物理尺寸/设备状态。用户映射计算值时必须可靠，且保留真实来源。

### F. 数据输入/转换（9）

`read_blast`、`read_fasta_lengths`、`read_gff3`、`read_bed`、`read_paf`、`read_genbank`；`as_ggchord_seq`、`as_ggchord_feature`、`as_ggchord_link`。

- 所有读取器逐个审查单文件/多文件、来源字段、压缩、异常行、坐标语义、NA/空值、格式歧义与可选依赖。
- 三个适配器要覆盖 data.frame / 常见 Bioconductor 结构，区分 link 和 alignment，保留 `GRanges` / `GInteractions` 原始元数据。

### G. 数据计算、工作流与输出（26）

`build_ggchord_cds`、`bundle_ggchord_ribbons`、`clean_ggchord_data`、`deduplicate_ggchord_ribbons`、`diagnose_ggchord`、`export_ggchord_layout`、`filter_ggchord_ribbons`、`filter_restriction_sites`、`find_common_features`、`find_ggchord_loci`、`find_primer_bindings`、`find_restriction_sites`、`focus_ggchord_data`、`focus_ggchord_track`、`get_chord_layout`、`link_ggchord_clusters`、`map_ggchord_position`、`merge_ggchord_ribbons`、`optimize_ggchord_layout`、`prepare_ggchord_data`、`project_ggchord_features`、`project_ggchord_links`、`segment_ggchord_data`、`validate_ggchord_data`、`validate_ggchord_hierarchy`、`view_ggchord`。

- 逐一测返回值是纯数据、报告、图对象还是布局结构；不得暗中改变原对象；区分“建议布局”与“执行布局”。
- 特别审查 `filter`/`deduplicate`/`merge`/`bundle`/`stat` 之间的数据丢失、来源集合、方向保持和聚合统计的不同语义。
- 验证 `focus_*` / `segment_*` / `map_*` 处理源坐标、显示坐标、窗口和断点映射时的可逆性。

### H. 其他自有 API 与第三方重导出

自有：`primer_catalog`。

第三方 22：`PositionIdentity`、`aes`、`after_scale`、`after_stat`、`annotate`、`arrow`、`element_blank`、`element_line`、`element_rect`、`element_text`、`expansion`、`ggsave`、`ggtitle`、`guides`、`labs`、`last_plot`、`margin`、`position_identity`、`rel`、`theme`、`unit`、`waiver`。

- 必须逐个说明保留的便利性、移除后的用户代码、是否可以直接使用 `ggplot2::` / `grid::`、与 IDE 补全及其他包冲突的风险。

### I. S3、ggproto、隐式注册与返回契约

S3：`+.ggchord`、`ggplot_build.ggchord`、`ggplot_gtable.ggchord_built`、`print.ggchord_clean`、`print.ggchord_diagnostics`、`print.ggchord_validation`、`print.ggchord_validation_summary`、`summary.ggchord_validation`、`as.data.frame.ggchord_validation`。

- 不止罗列函数：检查 S3 dispatch、ggproto 图层类/Position/Stat/Geom/Guide、`.onLoad` theme registration、未显式注册的动态对象、对象结构稳定性。
- `ggplot_build`/`ggplot_gtable` 的转换、缓存和设备行为是发布阻断级风险；输出必须在重复构建和多图并行时正确。

## 5. 参数级审计的强制数据表结构

Codex 必须生成 **一行一参数** 的机读 `api_arguments.csv`，而不是一行一个函数的粗略清单。字段至少包含：

`function, source_file, is_exported, category, parameter, default_expr, type_contract, allowed_values, unit, null_semantics, na_semantics, applies_to_coords, applies_to_geoms, validation_location, forwarded_to, behavior_source, actual_usage, duplicated_with, deprecation_status, proposed_name, proposed_default, decision, rationale, test_case, migration_example, confidence`。

审计还需生成 `api_functions.csv`、`api_aesthetics.csv`、`api_return_contracts.csv`、`api_internal_dependencies.csv`、`api_reexports.csv` 和 `api_compatibility_matrix.csv`。同一类参数（如 `offset`、`gap`、`length`、`ratio`、`size`、`width`、`colour`、`side`、`orientation`、`validate`、`bounds`）应横向比较，而不是只在各文件里局部判断。

**特别检查**：`mapping` / `data` 组合、`inherit.aes`、`show.legend`、`...` 的传递与拼写错误；固定值和映射值冲突；不同 coord 下是否悄悄无效；`NULL` 是自动还是禁用；语义类似参数有不同默认；绝对单位 vs 相对布局单位；用户误填时是否报错；多图层组合下参数是否互相覆盖；新图层是否无需修改内核巨大 `switch` 才能加入。

## 6. 目标 API 的建议结构（草案，非自动执行的定稿）

### 6.1 三个主要外部入口

1. `ggchord(...) + ... + coord_chord()`：多序列弦形同源/事件连接。
2. `ggchord(...) + ... + coord_circular()`：单序列圆形基因组/质粒及密集注释。
3. `ggchord(...) + ... + coord_collinear()`：直线多序列共线性及分段比较。

尽量让序列/feature/link/track 的数据契约与标签策略跨模式共享；对于确有不同用户心智模型的专用便利 geom，可保留简短的 wrapper。禁止因为三个坐标都支持某些 layer 就强制在不合理场景里产生伪效果。

### 6.2 将接口分层呈现而不是一味减少导出数

- **基础层（推荐最先学）**：构造器、三坐标、`geom_seq()`、gene/feature/link、标签、标准美学与 ggplot2 常用函数。
- **常用进阶层**：Position、Seq axis/scale、track、研究中常见数据导入、`focus_*()`。
- **专业层**：限制性酶切、引物、CDS 投影、stat ribbon、复杂分段/聚合、层输出。
- **内部层**：路径搜索、collision registry、缓存、候选布局、设备注入、多阶段局部状态，不列为稳定公开 API。

用户应通过网站 reference index 一眼找到“如何画基因？”“如何让标签不碰撞？”“如何展示结构变异或同源区间？”“如何加定量轨道？”；不要求先学习包的内部对象名。

### 6.3 新 API 的进入门槛

必须满足以下 7 项才允许加入：

1. 有一个真实、可复现的科研场景，且现有 API 不能自然组合实现；
2. 能一句话解释“做什么”，并有直观名字；
3. 没有语义等价的旧入口，或附完整合并替代方案；
4. 有明确输入/输出、坐标约定及错误行为；
5. 不破坏已有三坐标模式与图层独立性；
6. 有最小可运行示例与针对性回归验证；
7. 有将来继续扩展而不新增专有特例参数的机制。

潜在新功能（诸如更多定量 geom、多 contig 组合、复杂结构事件表达、分组布局）只纳入 **待研究池**，先做现有 API 的 gap analysis；不因本文件提及就自动实施。

## 7. 分版本路线：从 v0.16.0 开始重构

> 版本号是建议的工作包编号，可根据改动量合并或拆分。绝不能等到 v1.0.0 才开始改变 API；**文档与网站的全面改版也必须在 v1.0.0 发布之前完成并验收**。

### v0.16.0 — 架构契约审计 + 基础 API 清理（开始重构）

**必须完成**：清点全部公开/内部函数与所有参数、建立调用/参数流/返回结构清单；明确用户心智模型和数据契约；修复已确认错误（例如 F01/F02）；梳理 deprecated 参数、直接别名、重复入参、不同坐标下差异；选定第一批低风险改名、合并、精简；公开迁移表。新功能仅限能修补关键抽象缺口者。

**不做**：大规模重写所有布局算法、不加第四种绘图大类、不重建整个网站。**交付**：接口基线 CSV、审计报告、第一批迁移记录、可运行最小例子、针对性回归、版本设计文件。

### v0.17.0 — 绘图语法/可发现性统一

**目标**：`geom`/`stat`/`position`/`coord`/`scale`/`guide`/`theme` 各层归属与参数命名；统一 label 和 link/ribbon 语义；确定固定/映射参数优先级；审查全部便利 wrapper、颜色拼写别名和重导出。若真实场景需要，在这一阶段引入符合统一架构的新 geom/Position，不能绕开新规则。

**验收**：三坐标 × 核心图层的兼容矩阵，规范化绘图示例，可迁移的旧调用，图层顺序/重叠不引起不合理状态污染。

### v0.18.0 — 数据 API、分析工具和生物学契约

**目标**：统一 `read_*`、`as_ggchord_*`、validate/clean/prepare、focus/segment/project、alignment/link、cluster/primer/restriction/common-feature 的语义；明确元数据、来源与错误模型。必要时引入稳健的分组/多 contig 或研究场景补充，但新接口必须先通过语义审查。

**验收**：每个读取器与适配器的坐标不变式；跨原点、未知链、反向、零宽、混合分段、缺失注释等边界；失败必须可追溯且不可假装成功。

### v0.19.0 — 内部引擎与稳定输出契约重构

**目标**：把巨型布局函数与 `eval(quote())` 隐式状态转换为明确的输入/输出布局阶段；统一图层注册、标签避让、link routing、track、设备尺度和渲染缓存的所有权；确定 `export_ggchord_layout` 及诊断的正式结构。

**实施约束**：先构建原版与候选版的可比渲染语义测试，再小步抽取纯函数；每次变化不得同时大改几何结果与用户 API。确认性能与内存退化后才优化；不为追求代码行数减少而重写已经正确的算法。

### v0.20.0（可选）— API 候选冻结与迁移整合

如前几阶段出现较多破坏性变更，设置独立稳定化版本：统一错误信息和废弃策略、减少重复重导出、冻结默认参数/对象字段、修复三种绘图的组合缺口，形成 1.0 的候选接口。若工作量小可以并入 v0.19.0；无须强行增加版本。

### v0.98.0 / v0.99.0 — 文档、网站、发布候选验收

**文档必须在 v1.0.0 前完成**。稳定 API 候选确定后，先完善最小教程、逐主题进阶教程、双语 reference、迁移指南、FAQ 和可执行代码；再重建 README 与 pkgdown 风格网站，生成所有示例图片，并做用户入门测试。v0.99.0 作为 release candidate：兼容性、视觉、性能、README、网站、双语内容与迁移表均验收通过，不再主动扩张范围。

### v1.0.0 — 仅正式稳定发布

不在本版本开始架构重构、首次重写网站或大量新增公开 API。只处理 RC 发现的 release blockers；冻结 API 与参数语义、发布权威参考文档、对应版本网站/源码/tag，明确发布后兼容政策。若还存在明显语义冲突或文档落后，延后 v1.0.0 而不是带债发布。

## 8. 各版本任务的强制执行步骤（Codex 作业指导）

**步骤 1：读取** `AGENTS.md` → 总路线图 → 当前版本详细设计文档 → 实际代码 → 当前 `NAMESPACE` 与测试 → 相关 `NEWS` 历史。

**步骤 2：基线冻结**：记录 `git rev-parse HEAD`、`DESCRIPTION` 版本、导出数、源码数；生成配套 `audit_ggchord_api.R` 的清单。所有报告保留所依据 commit SHA。

**步骤 3：逐 API 证据卡片**：每个函数/参数记录用户需求、现状、调用处、所有 coord 中实际表现、源码/运行证据、是否重复、推荐方案、至少一个替代方案、用户可见迁移示例。

**步骤 4：优先级分类**：
- `P0`: 生物学坐标、数据损坏、静默误导、重大重复 build 状态污染等正确性风险；
- `P1`: API 语义冲突、明确错误、阻断功能组合、公开设计核心缺陷；
- `P2`: 命名、默认、别名、过度暴露、冗余便利；
- `P3`: 美观性、文档补全、极低频增强。

**步骤 5：变更决策**：`KEEP / FIX / RENAME / MERGE / SPLIT / INTERNALIZE / REMOVE / ADD / DEFER`。对 H/D 不得无用户审批即移除/改名；语义变更、默认外观变更和不可逆输入变换必须特别评估。

**步骤 6：实施策略**：以小 PR 或小提交拆分，先补充最小反例，再改源码，再迁移示例/帮助，再核查旧接口；每批返回风险报告与修改矩阵。允许新 API 伴随重构引入，但必须先过第 6.3 节准入门槛。

**步骤 7：针对性验证**：普通迭代按 `AGENTS.md` 只运行必要的轻量测试/样例；正式发布授权后才跑完整打包/CRAN check/设备矩阵。视觉比较在相同 seed、尺寸、dpi、字体、ggplot2/R 版本下进行。任何预期视觉差异要附理由。

**步骤 8：迁移与文档**：`NEWS.md` 写破坏性变化，维护 `MIGRATION.md`（旧调用 → 新调用 → 变动原因 → 删除版本）；稳定前更新关键 Rd 示例避免误导，但网站/完整双语教程的大重建安排在 v0.98/v0.99，并在 v1.0.0 发布前验收。

## 9. 测试与验收矩阵

| 维度 | 必需代表场景 | 失败判据 |
|---|---|---|
| 三坐标 | chord / circular / collinear 基本 build + 输出 | 原图失真、越界、非法伸缩或异常图层遗漏 |
| 数据语义 | 1-based、0-based 导入、单碱基、反向、缺失、跨原点、分段 | 坐标静默变化、方向被推断或来源丢失 |
| 图层组合 | 多个同类 geom、不同 data/mapping、手动 scale | 同名图层覆盖、映射失效、图例错误 |
| label | 密集 feature、长名、同半径不同尺寸、多种 callout | 严重不合理遮挡、跨 coord 不一致、文字裁切 |
| link | 正反连接、非相邻连接、local orientation、线与 ribbon | 相对方向错判、锚点丢失、缺失数据被伪造 |
| track | 多 ID、多图层同 ID、缺失值、bar baseline、line closed | 数值域被物理布局污染、跨窗口误连 |
| 数据读入 | BLAST/GFF/BED/PAF/GenBank/FASTA/Bioc | 错误格式静默成功、非法推断、来源缺失 |
| preview/export | PNG/PDF/SVG、不同画布与 dpi、保存后重复 build | 缓存串扰、尺寸改变引发非预期坐标误差 |
| 公共迁移 | 所有重命名/删除 API 的成对前后代码 | 文档与实现不一致，用户缺少替代路径 |

验收不是逐参数无限膨胀 testthat。应采用**少量公开契约的表驱动案例 + 高风险缺陷的最小回归**，符合现有 `AGENTS.md` 的测试原则。内部重构阶段可使用短期诊断脚本，但不要把私有布局状态永久固定为公开接口。

## 10. 文档、README 与网站：API 稳定后、1.0 前完成

### 10.1 官网内容架构

参考 ggplot2/tidyverse 与 gggenomes 的**信息组织和使用者工作流**，不复制其品牌资源、文案或代码。建议结构：

1. **首页**：一句话定位 + 三类图的代表图片 + 最小示例 + 安装入口 + Get started。
2. **Get started（15 分钟上手）**：最少数据列 → `ggchord()` → 添加 `geom_seq()` → 第一条链接/feature → 切换三种 coord → 导出图片。
3. **Data model**：四类数据，坐标与方向，表格/Bioc 输入，常见错误。
4. **三套场景教程**：Chord / Circular plasmid / Collinear，各按“基础 → 美化 → 高级”逐层递进。
5. **图层与美学指南**：geom、stat、position、scale、guide、theme、coord 之间的职责。
6. **进阶专题**：track、标签避让、restriction/primer、CDS projection、focus/segment、导出诊断。
7. **Reference**：按用户任务分类，所有函数/参数及默认值实时从最终 R 源码生成。
8. **FAQ & troubleshooting**：数据读不进、坐标不对、标签拥挤、导出字体和尺寸差异。
9. **Migration & changelog**：从 0.x 升级的代码迁移、改名、默认变更与删除。
10. **Gallery**：真实科研案例、完整输入与生成代码，不用只有漂亮但无法复现的静态截图。

### 10.2 README 的边界

README 应短、准确、好读。推荐：定位 + 1 个最小示例 + 3 类图的各 1 张图片/简述 + 安装 + 网站入口 + 开发状态 + 引用/许可。**不要**在 README 重复全部函数列表、参数说明、详尽教程和复杂比较表。

### 10.3 教程教学原则

由一个新手无需提前理解 ggproto、staged aesthetics 或环形坐标几何也能完成的例子开场；每节只引入必要的新概念；代码必须可复制运行；每个教程使用完整且可重建的输入、正确的文件路径、明确预期输出；从概念引向使用，再解释进阶选择；最后提供“为什么该用 A 而不是 B”的决策指南。

### 10.4 双语与版本一致性

英语与现代汉语对应页面一对一，代码和图片共用生成源；术语表统一 feature、gene、link、ribbon、strand、orientation、track、aesthetic、scale 的翻译或保留规则。网站导航、search、暗/亮主题、移动端、色彩可访问性、标题层级、API reference 可搜索性应达到长期维护要求。对 ggchord 的物理尺寸敏感标签，图像必须通过固定 export-size 渲染流程生成，不使用 RStudio 的随机面板尺寸作为范例。

### 10.5 官网与文档可核查参考

- ggplot2 官方扩展文档：https://ggplot2.tidyverse.org/articles/extending-ggplot2.html
- gggenomes 网站与教程：https://thackl.github.io/gggenomes/
- tidyverse style guide：https://style.tidyverse.org/
- tidyverse 文档规范：https://style.tidyverse.org/documentation.html
- pkgdown 网站：https://pkgdown.r-lib.org/
- Bioconductor 基因组对象：https://bioconductor.org/packages/GenomicRanges

## 11. Codex 实际交付物及目录建议

每个版本在其分支下维护可追溯产物（具体路径可依项目既有规范整合，不要在根目录无限增加文件）：

```
design-roadmap/
  API-ARCHITECTURE-PRINCIPLES.md         # 设计原则和公开边界
  API-AUDIT-DECISIONS.md                 # 各函数与参数的决策台账
  API-MIGRATION.md                       # 旧 → 新 API 的可执行迁移示例
  versions/
    v0.16.0.md                            # 新范围与验收
    v0.17.0.md
    v0.18.0.md
    v0.19.0.md
    v0.99.0.md
  audit-generated/                       # 脚本生成、可重新计算的台账
    api_functions.csv
    api_arguments.csv
    api_reexports.csv
    api_internal_dependencies.csv
    api_return_contracts.csv
    api_compatibility_matrix.csv
```

**每次交付必须回答**：修改了哪些 API？改变了哪些默认？是否破坏旧调用？有哪些替代方式？哪几个真实案例验证？目前还有哪些 H/D 尚未确定？用户下一步需要作出的决策是什么？

## 12. 给 Codex 的启动指令（可直接复制）

> 你正在维护 `DangJem/ggchord`。**必须先通读仓库 `AGENTS.md`、当前 `design-roadmap/`、本审计规范，并以最新 `main` 的实际代码作为事实来源。** 项目已经完成 chord、circular、collinear 三种坐标大框架，不再添加新的大类。所有重构从 **v0.16.0** 开始，并分多个版本完成；v1.0.0 只是 API、文档、网站、示例、验收全部定型后的正式稳定版。
>
> **第一项任务是“审计与设计”，不是立即改写全部源码。** 在新的版本分支上清点全部导出、重导出、S3、ggproto、内部函数及逐一参数；使用提供的 `audit_ggchord_api.R` 生成完整台账，手工检查动态注册、`eval(quote())`、非标准求值、别名、跨坐标行为和真实工作流。对每个问题给出源码证据、严重程度、KEEP/FIX/RENAME/MERGE/SPLIT/REMOVE/INTERNALIZE/ADD/DEFER 建议与替代设计。未经证实的判断必须标 H/D；任何破坏性 API 方案先提出供用户确认。
>
> 先提交新的 v0.16.0 分阶段重构计划、API 总目录、问题清单和第一批**低风险**修复/整理候选，等待确认；不要直接大规模重写、自动 push、创建 GitHub Release 或跑完整发布验收。每批修改应含最小可复现案例和可执行迁移指南。文档与网站的完整双语重建必须在 v1.0.0 正式发布前完成，README 保持简洁，详细内容放在网站与专题文档。

## 13. 审计局限与复核要求

本文件基于已读取的 GitHub `main` 源码定义、公开接口、关键实现和项目规范，对已定位的问题提供初步证据；**未在本次会话中重新安装/执行整个 R 包、测试所有参数组合、构建完整动态调用图或穷举确认每个内部函数无引用**。因此，`F01`、`F02` 等有具体源码证据，其他 `H/D` 均需 Codex 取得完整 checkout 后用静态 + 运行验证再定案。尤其不允许把“候选无用/冗余”直接等同于“安全删除”。

审计完成的客观定义：**所有 130 个导出名称、9 个 S3 注册点、全部内部函数及逐个形参有台账；所有公开 API 有明确去留决策；所有修改有完整迁移方案、最小验证和主场景可用性评估；所有未知事项有公开记录。**

---

**End of audit handoff.**
