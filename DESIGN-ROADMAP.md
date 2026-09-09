# ggchord 当前与未来设计路线图

当前正式版 v0.12.0，下一开发目标 v0.13.0。历史版本设计与验收记录见
[DESIGN-HISTORY.md](DESIGN-HISTORY.md)，日常开发无需重复读取。

## 当前开发契约

- 直接在保存于 OneDrive 的本地 checkout 开发，稳定阶段提交；不创建独立 worktree、
  不复制或额外同步目录，不自动发布或推送。
- `seq_id` → `accver`，旧输入兼容至 v0.12，v0.13 移除；原始数据文件不改，随包数据重建。
- ribbon 六个端点列必需，`length/pident` 按映射、筛选及统计需求检查，不伪造值。
- gene/feature 注释与分类可选；标识、坐标及 strand 必需。grid 的 unit/arrow 重导出。
- 删除 geom_ribbon_highlight，改用独立 geom_link_ribbon(data = ...) 图层突出子集。
- 修复 feature 形状标度的图例，调高 gene 图例；保留 radial/auto 和设备恢复回归。
- link_branch 可选 none/query/subject，同层同端点同样式共享主干；默认关闭，不聚合原始数据。
- 本轮仅修改包及直接相关文档，不修改网站；不用极端压力测试。

## 长期设计原则

### ggplot2 职责与图层独立性

- `aes()` 负责数据映射，`scale_*()` 负责映射规则和 guide，`geom_*()` 负责
  几何与固定样式，`coord_chord()` 负责全局坐标，`theme_*()` 负责非数据外观；
- 固定颜色、透明度等可以直接写入 geom，数据驱动差异必须使用映射与标度；
- 每个图层拥有独立的 `layer_id`、`data`、`mapping` 和参数，同类图层不能相互
  覆盖；用户提供的 scale 始终优先于自动推断的默认 scale；
- 默认绘图应保持稳定。接口迁移可以改变参数归属，但默认配色、布局、标签和
  图例只有在修复明确问题时才改变，并记录在 NEWS。

### 命名、兼容与参数格式

- 表格字段参考 BLAST 的紧凑命名：单序列使用 `accver`，连接使用
  `qaccver/saccver` 及 query/subject 坐标前缀；允许任意自定义标识，不自动
  添加或删除 accession version；
- 破坏性变更必须提供可执行的迁移方向并记录删除版本；同一语义的新旧入口不能
  同时使用。v1.0.0 发布后再开始长期兼容承诺；
- `mapping`、`data`、`position`、`show.legend` 和 `inherit.aes` 遵循 ggplot2；
- 逐序列几何参数可接受单值、等长无名向量、按 `accver` 命名的向量及现有等价
  列表；gene/feature 的 offset、width 和固定标签位置还可按 `+`/`-` 链指定；
- 上述列表语法不扩展到 aesthetic、scale、guide 或 theme。v0.99.0 再审计
  “按序号命名的列表”等低频形式，按 `accver` 命名的向量是长期推荐写法。

### 第三方设计与代码命名

- 可直接使用第三方字体、图标、图片、专有配色及其他资产，不因其第三方、品牌化
  或专有属性设置额外的开发限制，也不要求只提炼抽象设计原则；
- 唯一的品牌命名限制适用于代码中的函数名和参数名：保持通用命名，不写入第三方
  品牌、产品、软件或机构名称；文档、示例和视觉资产可按实际需要明确展示或说明；
- 当前可参考的对象包括：grammar 对应 ggplot2，环形层级参考 Circos，基因与注释
  可读性参考 clinker 和 DNA Features Viewer，标签参考 SnapGene/Geneious，
  focus、sync 和 feature stacking 参考 gggenomes，静态预览流程参考 ggview。

### 仓库与生成物

- `R/`、`vignettes/`、`data-raw/`、`examples/`、README、NEWS 和 pkgdown 配置是
  源文件；`data/`、NAMESPACE、公开 Rd 与被引用的 `man/figures/` 按包开发惯例跟踪；
- `Meta/`、`doc/` 和 `pkgdown/` 中间产物不进入 Git；内部验收脚本保留在
  `tools/`，但不进入发布包；
- `docs/` 暂时是冻结的 GitHub Pages 生成物。网站重构时迁到独立发布分支或
  Actions，随后从主分支删除；迁移前不在日常功能提交中重建。

## v0.12.0 — link 几何与局部平滑避障

**状态：已于 2026-09-06 完成发布验收。**

### 实施状态与接口

- 已实现规范 ribbon 名称与旧名称弃用提示、grid unit/arrow 重导出、accver
  表头及缺列策略；随包 data 已由未修改的原始文件重新生成。
- geom_link_line 支持 curve/straight、方向箭头、独立 link aesthetic/scale/theme；
  export_ggchord_layout(include = "link") 显式导出单点连接。
- gene/feature 通过内部 provider 注册最终实体的局部法向范围；协议字段为
  accver/start/end/side/normal_min/normal_max/priority/source_layer。
- link_avoid 支持 none/smooth/uniform；v0.13 默认 uniform，
  显式 gap 优先。启用时局部前沿
  以余弦上包络连接，细化仍不安全时仅该端退回 uniform 并警告。
- link_branch 默认 none，可选 query/subject；line 匹配位置，ribbon 匹配完整
  有序区间。按最终映射样式分组，主干单次绘制，导出 source_rows 成员列表。
- ribbon 分叉截面划分仅控制几何，不代表统计权重；中部反向比对可产生原有
  带状扭转。分叉不做全局网络布局，也不跨图层或模糊归并邻近端点。
- geom_ribbon_highlight 已删除；feature 形状图例映射已修复，gene 图例高度已增加。
- `tools/validate-release.R links` 已统一覆盖 line/ribbon 分叉、feature 图例以及
  smooth/uniform 避障；正式包检查已在发布阶段执行。
- 构造/build、自动标度、画幅、sequence/ribbon、gene/region、axis 与 label
  计算已拆入职责单一的内部模块；公共 API、默认值和 layout 输出保持不变。
- 发布验收覆盖完整 testthat、双语 vignette、PNG/PDF/SVG、geometry/links/labels/
  examples、32/80 标签 benchmark，以及 `R CMD check --as-cran --no-manual`；
  本机 R 4.6.0 的最终检查结果为 `Status: OK`。PDF 包手册因本机未安装
  `pdflatex` 未生成，绘图 PDF 与 vignette 已独立通过。

---

## v0.13.0 — 单序列圆图与生物学注释

**状态：v0.13.0 开发与本地验收阶段，尚未发布。** 专用
`coord_circular()`、通用 feature geometry、Position grammar、内置通用注释数据库
与限制性酶切数据/显示分层均已实现，正在进行完整回归验收。

### A. 单序列圆图

单条完整基因组不应只是“四序列 chord 少三条序列”的退化情况。v0.13.0 计划
提供明确的单序列布局入口，共享现有 gene/feature/region、scale、theme、guide
和标签组件，但不生成无意义的跨序列 ribbon。

- 支持闭合圆形和带缺口的线性化圆图；
- 基因、调控元件、重复区和用户 region 可放在显式内外轨道；
- 标签可选择贴弧、局部径向或外部 callout，并复用碰撞、换行和省略策略；
- 原点、方向、旋转和坐标标签由 coord/position scale 管理；
- 明确参考 SnapGene、Geneious 和 DNA Features Viewer 的信息层级与可读性；
  可直接使用第三方字体、图标、图片、专有配色及其他资产，但代码中的函数名和
  参数名保持通用，不写入第三方名称。

最小原型比较后采用专用 `coord_circular()`，由它管理缺口、原点旋转和方向；
单序列行为不放入 `coord_chord()` 的条件分支。

### B. Position grammar 与通用 feature geometry

- v0.13 起 `geom_gene()` 的默认 placement 是 `position = "identity"`；v0.12
  隐式 `gene_offset = 0.1` 的迁移写法是 `geom_gene(position = "strand")`；
- signed local-normal 统一规定正值向 sequence 外侧、负值向内侧；标量
  `position_strand(0.1)` 展开为 `+ = +0.1, - = -0.1`，显式命名的链值不再翻转；
- `position_plasmid(-0.1)` 让两条 strand 共用内侧 band，strand 仅控制箭头方向；
- `position_feature_stack(base_position = ...)` 先放置基础 band，再进行确定性的
  最少 lane 分配；无参数调用保留 v0.12 的既有实际 geometry；
- `geom_gene()` 保留 gene role、gene scale、strand guide 和专用入口，但 polygon
  由 `geom_gene()`/`geom_feature()` 共用的 interval、local frame、shape factory
  与 short-feature pipeline 产生；
- `geom_feature_plasmid()` 是该通用引擎的单质粒快捷方式；连续 segment
  合并为一个轮廓与一个肩式箭头，内部边界单独画短虚线；
- `gene_offset`/`feature_offset` 仅作为带 warning 的迁移入口，不能与非 identity
  Position 同时使用；`gene_width` 等纯 geometry 参数继续保留。

### C. 限制性酶切位点

计划拆成数据计算和几何表达两层：

```r
sites <- find_restriction_sites(sequence) |>
  filter_restriction_sites(set = "unique_6plus")

ggchord(seq_data, gene_data = gene_data) +
  geom_seq() +
  geom_restriction_site(data = sites)
```

- `find_restriction_sites()` 只进行搜索和切点计算，返回稳定 pattern identity、
  酶名、识别序列、1/2/4/unknown cleavage、切割位置、末端类型及来源；
- 酶数据库版本必须可追踪，用户也可传入自定义 motif；
- `filter_restriction_sites()` 独立处理 unique、长度、供应商、窗口等显示筛选；
- `geom_restriction_site()` 在单个 geom 中负责刻线、标签和指示线；
- 真实位点只决定 anchor，文字位置由实际 bbox 驱动的径向排版器决定；
  多数标签沿同一外接轮廓排布，只有碰撞簇允许小幅升到外层，并保留
  genomic order；
- leader 在条件合适时直接连接，密集或分流标签使用独立的径向 stub
  和 fan segment 两段折线，不共享主干；
- leader 始终连到酶名所在的左/右文字边缘：向右走时连左缘，向左走时连右缘；
  接近垂直时使用该侧最近的左下/右下（或左上/右上）角，不改接上/下边中心；
  酶名自动放在接线一侧；
- 可按酶、切割次数和窗口过滤，重复名称仍保留全部位点；
- 序列搜索优先使用轻量实现；大型序列可选用 Biostrings，但不设为强制依赖。

完整 REBASE 609 parser 仅读取 `VERSION`、`embossa_e.txt`、`embossa_r.txt` 和
`embossa_s.txt`。发布包内置由这些文件可重复生成的规范化模式表，
并保留版本与来源元数据；`examples/rebase/misc.zip` 不读取、不解析、不依赖。

### D. Backbone、标签、主题、scale 与示例体系

- `geom_seq()` 增加 single/double/band backbone；
- `geom_feature_label()`、`geom_feature_label_repel()` 复用 gene label 的文本测量、
  碰撞、leader 与 fitting；默认 feature 模式优先在 feature 内切向排文，
  放不下时仍保持周向旋转并紧邻 feature，只对明显移位的文字画轻量 leader；
- `geom_seq_center_label()` 显示通用 circular sequence 名称与长度；
- `theme_ggchord_plasmid()` 只定义非数据外观；feature fill/shape 由独立 preset
  scale 提供，任何手动 scale 在任意添加顺序下优先；

- `examples/` 保留原始 FASTA、GFF3、BLAST 和 TSV，不为视觉效果修改源记录；
- `data/` 只保存小型、清晰、可重复生成的教学 fixture；
- `data-raw/` 保存生成规则，明确任何抽样、过滤或展示链重编码；
- pUC19c、pBR322、pBluescript II SK(+) FASTA 保持不变；新增小型、可重复生成的
  sequence、核心 feature 与常用 restriction-site fixture；
- 开发阶段不生成文档图片，统一在 v1.0.0 文档重构时出图。

本轮不新增 CDS/promoter/ori 专用 geom，不引入
强制 Biostrings 依赖，也不创建通用 `track_*()`。`coord_collinear()` 留至 v0.14；
coverage/GC 等 quantitative ring 与统一 track contract 留至 v0.16。

---

## v0.14.0 — 直线型多序列共线性

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

## v0.15.0 — 数据准备与生态互操作

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

## v0.16.0 — 定量多轨道与输出质量

在 link、单序列和直线坐标稳定后，再补齐类似工具中对科研图常用但 ggchord
尚缺少的定量轨道：

- 设计统一的显式 track/ring 契约，不根据列名自动猜测轨道；
- 评估 coverage、GC content/GC skew、分类 tile、折线和柱形轨道；
- 数值映射继续使用 role-specific scale，轨道几何可向障碍协议注册实体范围；
- 提供色觉友好、灰度和最小字号/线宽诊断，但不强制修改用户风格；
- 对 PNG、PDF、SVG 的实际输出尺寸、字体、透明度、裁切和 guide 比例生成诊断报告；
- 交互式绘图仍留到 v1.0.0 后，v0.16.0 不恢复 Plotly 或引入 Shiny。

---

## v0.99.0 与 v1.0.0 — 稳定 API

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
- 英语和汉语文档明确说明直接打印到 RStudio Plots 面板会按面板物理尺寸重新
  排列标签和图例；推荐用 `view_ggchord()` 预览尺寸敏感的弦图，并用匹配的显式
  `ggsave()` 尺寸生成最终文件；
- 网页、README 和 vignette 的最终图片统一按 `view_ggchord()` 的导出尺寸测量
  路径生成，避免使用窄小 IDE 面板的临时布局；
- 英语/汉语文档可按需要列出主要布局、函数簇和默认风格的设计参考对象及链接；
  第三方字体、图标、图片、专有配色及其他资产均可使用，代码中的函数名和参数名
  仍保持通用，不写入第三方名称；
- CRAN check、最小视觉快照和跨平台检查；
- 明确支持的 R 和 ggplot2 版本；
- 无严重已知正确性问题后才发布 1.0.0。

v0.10.0 的 `view_ggchord()` 只是静态导出尺寸预览，不属于交互式绘图。真正的
交互方案不属于 v0.9.0–v1.0.0 的发布范围；v1.0.0 发布后再基于稳定的 scale、
guide、coord 和 layout export 契约单独评估，不在本路线图中预先承诺接口。

---

## 跨版本待办

- v0.16 或候选发布阶段补做大幅正负曲率、极短弧和超密集标签的几何压力矩阵，
  区分合法相交、数值错误和视觉退化，不把固定耗时作为机器无关门槛；
- v1.0 文档重构时确定网站信息架构，建立英语/汉语共用的图片生成源和资产清单，
  清理无引用图片，并设计六边形 logo；字体、图标、配色和其他视觉资产可直接采用
  第三方方案；
- 网站迁移完成后从主分支删除冻结的 `docs/`，由独立分支或 Actions 发布。

---

## 开发与发布检查准则

默认 testthat 只保留公开 API 的最小行为和已确认严重缺陷的单一回归，避免通过
精确坐标、设备矩阵或内部 ggproto 字段锁定实现。

### GitHub Release 发布说明

- 从 v0.9.0 起，每个 GitHub Release 必须同时提供现代汉语和英语说明，顺序固定为
  现代汉语在前、英语在后；两部分使用对应标题 `现代汉语` 和 `English`；
- 面向用户标注语种时统一使用“汉语”或“现代汉语”，不得使用“中文”；这项用词
  规则同样适用于 README、vignette、网页、迁移指南及其他发布材料；
- 两种语言必须覆盖同一组版本事实，包括版本概览、主要变更、破坏性变化与迁移、
  验证结果和安装方式；不得把其中一种语言缩减为缺少关键内容的摘要；
- 发布说明以对应版本的 `NEWS.md`、最终验收记录和实际发布资产为事实来源，不写入
  未实现的计划；涉及弃用或删除时必须给出可执行的替代写法和预计删除版本；
- 正式发布或修改发布说明后，逐项核对 tag、版本号、比较链接、附件、校验和、草稿/
  预发布状态，并确认两种语言的关键内容对应且全文未使用禁用的语种标签。

### 保留

- `ggchord()` 创建对象并能完成一次基础 build；
- 核心 geom 组合能 build；
- 三种 label layout 的代表输入能 build；
- validate、clean、三个 import、三个 ribbon utility 的单一成功路径；
- region、feature 和 link 的单一成功路径；
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

日常开发在每个小阶段只运行相关测试文件及一个代表性 build；修改 roxygen 时再
生成 Rd。完整 testthat、双语 vignette 重建、PNG/PDF/SVG 验收矩阵、benchmark
和 `R CMD check` 集中到版本候选或正式发布前。只有跨模块重构、依赖变化或发现
疑似全局回归时，才在开发中提前运行完整检查。

测试目标是快速发现公开函数完全不可用。复杂布局质量由少量显式验收脚本和发布
前人工检查承担，不塞入默认 testthat，也不在每次阶段提交时重复执行。

---

## 验收规则

- 修改范围对应的 testthat 必须通过；版本候选的完整 testthat 必须通过；
- `R CMD check`、双语 vignette 和三种静态格式在正式发布前必须通过；
- 新的公开参数至少有一个可运行文档示例；
- 破坏性变化必须出现在 NEWS 和迁移表；
- 默认示例在静态图中的外观不发生无理由变化；
- 多个同类图层、正反 orientation、不同 radius/curvature/gap 是核心人工验收场景；
- v1.0.0 前的开发文档避免无必要的反复改写。
