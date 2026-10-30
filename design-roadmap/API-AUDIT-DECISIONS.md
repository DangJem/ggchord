# API 审计证据与决策台账

**状态：v0.16.0 第一轮源码复核；不是全部参数或动态行为审计的完成记录。**

审计起始提交为 `69b907a6ddffea60a73ed1d4b1392c8cc5f7da22`，其 R 源码与附件所用
`bdbb3326847fcbb71928eed41647c209bd545871` 一致。下列行号固定指向这个修复前基线；
后续改动应按函数名定位并在实际实施记录中注明新的源码状态。
范围与版本编号见 [分阶段重构计划](API-REFACTOR-PLAN.md)。

附件的 30 项记录保留原编号、优先级及 C/H/D：C 是报告声称已经源码定位的事实，
H 是待验证假设，D 是设计取舍。这里另列本地静态结论、缺少的运行证据和建议；
“事实吻合”不等于“缺陷已复现”，建议 FIX/KEEP/DEFER 也不等于用户已批准改名或删除。
本表不能替代逐导出与逐参数清单，未完成的语义栏仍需逐项审阅。

## F01–F10：清洗、布局与入口职责

### F01 — P1/C：清洗输入被当成未提供

- **静态证据**：`R/clean.R:74–127` 的 `clean_ggchord_data()` 先归一化 gene ID，
  随后以 `is.data.frame()` 决定 ribbon/gene 输出，否则赋 NULL；没有对这两项的
  非 NULL 非表格类型建立对称错误。类型不支持与用户未提供因此可能不可区分。
- **建议与缺证**：FIX；优先取得不支持输入的最小复现，并保留 NULL 与 data.frame/
  tibble 路径。不将 `ggchord()` 的 Bioconductor 自动转换推定为清洗器已有契约。
  本轮已复现字符串输入得到 NULL 且空报告，并已修复为具体的表类型错误；局部直接
  断言核对错误、NULL 和 tibble 路径通过。没有运行完整 testthat。
- **阶段**：v0.16.0。错误信息应指出具体表及转换方向，不静默损失输入。

### F02 — P1/C：反向 ribbon 排序报告的新值错误

- **静态证据**：`R/clean.R:297–312` 的 ribbon sort 分支，在交换端点之前用同一
  `coords` 写入 original/new value。`R/clean.R:476–481` 的 gene 分支已分别写
  `st/en` 与 `en/st`，不能将同一错误笼统归给全部区间。
- **建议与缺证**：FIX ribbon 报告；排序后的新值从交换后的坐标构造，原相对方向
  和源行不变。本轮已复现并修复，局部直接断言核对输出坐标、direction、原输入和
  original/new 值通过，不以精确布局浮点作为断言。
- **阶段**：v0.16.0。不改变默认 `reversed_interval` 或 gene 方向语义。

### F03 — P1/H：clean/prepare 的反向端点默认差异

- **静态证据**：`R/clean.R:80` 默认 sort；`R/data-prepare.R:72–74,139–145`
  明确将 prepare 的清洗缺省设为 keep，并拒绝会使后续 ribbon 步骤丢失端点方向的 sort。
  这不是未注明理由的任意差异，已有正确性保护与帮助说明。
- **建议与缺证**：DEFER 默认值改变；先记录两入口的不同用途，再比较“保留差异并强化
  契约”和“统一为显式策略”的迁移成本。缺少全部 pipeline 方向保持的运行证据。
- **阶段**：v0.16.0 定义取舍，v0.16.2 整理工作流；不能为名称一致删除方向保护。

### F04 — P1/C：共享环境中的布局阶段

- **静态证据**：`R/layout.R:54–138` 的 `compute_chord_layout()` 具显式长签名，
  以自身 `environment()` 顺次 eval 七个延迟阶段；`R/layout-sequence.R:2` 的阶段
  为 quote 表达式。共享状态确实存在，但结构本身不证明数值或绘图错误。
- **建议与缺证**：DEFER 内部抽取到 v0.16.3；先列每阶段读写、顺序和返回依赖，
  以相同输入的公开数据保留与几何语义比较支持逐步改写，不直接替换全部引擎。

### F05 — P1/H：按类型汇总的单次布局与图层独立性

- **静态证据**：`R/ggchord-build.R:9–70` 的 single 函数以 switch 覆写局部参数；
  但外层 `compute_chord_geometry()` 在 `:1035–1105` 分配 layer ID、按 ID 分组并
  处理依赖，`:1474–1495` 重建图层时保留独立输入等字段。不能只看内层认定同类层覆盖。
- **建议与缺证**：DEFER 缺陷定性；需同类多层、不同 data/mapping、交换顺序、复制图
  和重复 build 的最小运行证据。实际发现的数据覆盖按严重度提前修复。
- **阶段**：v0.16.1 核对公开行为，v0.16.3 整理状态所有权。

### F06 — P2/D：theme 参数规模与动态预设

- **静态证据**：`R/theme.R:53–122` 显式列出标准与角色主题参数；`:126–147`
  通过复制 `theme_ggchord` 并替换 body 生成四个预设，保留共同 formals。
  `R/zzz.R:7–50` 注册自定义主题元素；数量多不自动说明冗余。
- **建议与缺证**：DEFER 公开拆分；先建立标准元素、真正自定义元素和角色图例矩阵，
  比较保留显式签名/改进分组文档与调整入口两方案。需核对继承和预设实际组合行为。
- **阶段**：v0.16.1；预设默认视觉应保持稳定。

### F07 — P2/D：getter 与稳定 export 的公开边界

- **静态证据**：`R/zzz.R:115–142` getter 返回 plot-owned 原布局；
  `R/layout-export.R:64–83,197–239` export 选择组件并附元数据、诊断等。
  其 `:18–52` 公开返回说明已经列出 registry 的 band/slot/leader 等字段。
- **建议与缺证**：DEFER 删除或内部化；先逐字段比较两个出口及真实调用需求。不能把
  已有公开返回说明中的字段直接改称“内部”后删除；收窄需迁移与删除版本。
- **阶段**：v0.16.3；可选 v0.16.4 收束迁移。

### F08 — P2/C：plasmid 是便利封装

- **静态证据**：`R/geom-feature.R:445–479` 的 `geom_feature_plasmid()` 构造默认
  stacking/plasmid position 后调用 `geom_feature()`，并不是独立复制整个绘图引擎。
- **建议与缺证**：KEEP 现有入口作为暂行建议；共享实现事实吻合，不应以重复名称
  为由删除。仍须核对参数透传和默认视觉，特别是 F11 的固定形状透传候选问题。
- **阶段**：v0.16.1，是否长期保留便利入口仍属于公开设计取舍。

### F09 — P2/D：gene 与 feature 的语义分离

- **静态证据**：`R/geom-feature.R:89–173,379–410` 处理方向性、形状与复合区间，
  且复用 gene 类型布局通道；`R/geom-gene.R:48–155` 提供方向基因入口。
  `R/data-convert.R:263–267` 明确未知链归 feature，不能据共享布局断言输入域相同。
- **建议与缺证**：DEFER MERGE/REMOVE；比较保留语义便利入口并共享实现与统一 geom
  后显式角色两方案。缺少定向 gene、未知链 feature、复合区间和两者标签的成对案例。
- **阶段**：v0.16.1；不得静默将未知链编成正链。

### F10 — P2/C：offset 已弃用，但删除阶段须协调

- **静态证据**：`R/geom-gene.R:21–23,65–72` 与 `R/geom-feature.R:219–222,287–295`
  已发弃用警告，且拒绝与非 identity Position 同时提供。帮助原已计划 v1.0.0 删除。
- **建议与缺证**：KEEP 现有警告和迁移方向，DEFER 删除；将
  `gene_offset = x` / `feature_offset = x` 对应到 `position_strand(offset = x)`，
  再核对旧示例与逐序列格式。删除前移到冻结阶段须记录对原删除承诺的修订。
- **阶段**：v0.16.0 建迁移记录，v0.16.1/v0.16.4 确认处置。

## F11–F20：参数、输入与坐标

### F11 — P2/C：固定 shape 双入口与 wrapper 透传

- **静态证据**：`R/geom-feature.R:256–282` 提供 `feature_shape` 与简写 `shape`，
  使用 `missing(feature_shape)` 拒绝双入口。映射仍是 `aes(feature_shape = ...)`。
  `:445–466` 的 wrapper 总是显式传两参数，调用者只给 wrapper 的 shape 也可能误判。
- **建议与缺证**：DEFER 首选名及别名删除；wrapper 误判单独作为 FIX 候选先最小
  复现。本轮已确认并最小修复：仅调用者提供 shape 时在 wrapper 解析别名，再向
  核心传递唯一固定值。`shape = "block"` 的圆图 build、用户手动 scale 以及明确
  双参数仍报错的直接核对通过。首选名称尚未决定，修复不改变两入口的现有设计。
- **阶段**：v0.16.0 核查关联错误，v0.16.1 决定命名。

### F12 — P2/H：标签命名相似但缺省并非完全相同

- **静态证据**：`R/geom-feature-label.R:26–48` 将 `label_*` 翻译为 gene 内部参数，
  固定 feature 标签缺省方向为 feature；`R/geom-gene.R:227–236` 缺省为 horizontal。
  repel 的 feature/radial/auto/callout 选择及 shrink 路径见前者 `:118–192`，
  gene 的 radial 缺省和 segment 参数见后者 `:396–409`。
- **建议与缺证**：DEFER RENAME/MERGE；先做含默认与适用坐标的语义对照，运行代表
  固定/避让标签例子。不能只删前缀而意外改变默认方向和布局。
- **阶段**：v0.16.1。

### F13 — P2/H：line/ribbon 相近几何参数有不同默认

- **静态证据**：`R/geom-ribbon.R:34–50,94–98` 的 ribbon 缺省 uniform 避障；
  `R/geom-link-line.R:8–16,29–39,119–134` 的 line 缺省 none、自动 gap 回退 .035，
  straight 还拒绝控制点或 branch。共同 branch_fraction 已共享校验函数。
- **建议与缺证**：DEFER 统一命名；先比较端点间隔/边缘控制及正反方向的含义，
  保留必要几何差异。缺少跨半径、曲率、避障与 branch 的代表性等义比较。
- **阶段**：v0.16.1；不据相近参数名合并两个 geom。

### F14 — P2/D：color/colour 是直接别名

- **静态证据**：`R/scales.R:17,132`、`R/link-scales.R:15,30`、
  `R/guides.R:132` 是直接绑定别名，不是五套独立行为。
- **建议与缺证**：DEFER 删除；比较统一文档首选 colour 同时保留 color，与分阶段
  弃用 color 两方案。需要用户负担、现有示例和下游使用证据；不默认认定别名有害。
- **阶段**：v0.16.1；任何删除有替代调用与删除版本。

### F15 — P2/D：22 个第三方重导出

- **静态证据**：`R/reexports.R:13–137` 列出 20 个 ggplot2 和 2 个 grid 符号，
  包括 `PositionIdentity`、`last_plot`、`unit`、`arrow`；NAMESPACE 同时导入并导出。
- **建议与缺证**：DEFER 去留；逐项记录便利用例、替代命名空间及迁移负担，不能把
  包内引用少等同为用户不用。可比较保留完整便利集合与仅保留入门必需集合。
- **阶段**：v0.16.1，稳定候选时核对最终清单。

### F16 — P2/H：file/files 已共用解析器，但存在双入口差异

- **静态证据**：六个 reader 采用 `ggchord_resolve_files()`；`R/data-import.R:10–53`
  已接受路径向量和 glob，两者同时提供时先解析 files、再将 file 原文拼到前面。
  `file` 的单入口实际也可为向量，帮助常称单文件；两者并存时的解析不完全对称。
- **建议与缺证**：DEFER 统一接口；先验证双入口/glob/重复路径和来源顺序，不将
  共有 resolver 误写成每个 reader 各自重复实现。保留双入口或选一入口均需迁移。
- **阶段**：v0.16.2。

### F17 — P2/H：GenBank 的外部解析边界已显式化

- **静态证据**：`R/data-genbank.R:5–64` 保存 Python 解析脚本，`:146–180`
  校验解释器、可选 jsonlite、临时文件并 system2 执行；帮助 `:67–95` 已说明不安装/
  下载、内存读取、坐标与 unsupported 政策。缺少成熟解析器不是当前事实。
- **建议与缺证**：KEEP 显式可选依赖边界，DEFER 改解析器；需环境缺失、非 ASCII
  路径、离线、解析器版本、模糊/远程/跨原点位置及进程失败的定向证据。
- **阶段**：v0.16.2。不因有外部进程就断言执行风险已发生或删除此入口。

### F18 — P2/H：track 已有冲突、数值域与设备保护

- **静态证据**：`R/geom-track.R:46–65,138–173` 校验 track Position、物理单位、
  baseline/closed；`R/layout-track.R:132–178` 拒绝不合条件的 closed、同 ID 位置
  声明冲突，clone 用户 scale，并将同 ID 值和 baseline 共同训练；`:291–325` 处理 bounds。
- **建议与缺证**：DEFER 缺陷定性；报告要求检查的多项保护已经实现，应验证是否
  正确组合，而不是重建一套。缺少同 ID 多图层、变换/OOB、缺测分段和实际设备度量
  在组合下的运行证据。
- **阶段**：v0.16.3；严重数值或坐标问题发现即提前修复。

### F19 — P2/H：三个 coord 的共有参数已有共同校验

- **静态证据**：`R/coord-chord.R:59–93` 定义 ratio、fit、expand、clip 和 limits；
  `R/coord-circular.R:20–40`、`R/coord-collinear.R:91–118` 委托 chord 校验共有部分，
  再添加专属方向/窗口规则。`R/ggchord-frame.R:3–54` 按 fit 和模式重建坐标。
- **建议与缺证**：KEEP 共同校验的暂行结构，DEFER 语义改动；生成参数矩阵并核对
  labels/geometry/manual、用户 limits、实际 aspect 和 clip 的代表输出。
- **阶段**：v0.16.1；不能把共用参数“完全没有统一”当作已确认缺陷。
- **本轮补充证据**：当前 ggplot2 4.0.3 的 CoordCartesian 已实现读取 ratio 的
  aspect 方法。一次手动限域的成对 build/gtable 核对确认 ratio 1/2 分别产生
  不同面板高宽约束，排除“只保存字段、完全不生效”的假设；不据此保证其他依赖版本。
  rotation 顺逆时针文案与标准正角矩阵的口径，以及重建 Coord 时 segments 字段
  的可见性，保留为待核查事项；未取得用户可见错误证据前不改几何或复制状态逻辑。

### F20 — P2/H：Bioconductor 适配已有保真约定

- **静态证据**：`R/data-convert.R:56–68,98–124` 保留元数据和 Seqinfo；`:263–267,
  301–326` 明确保留未知链来源且输出 none；`:331–347,363–402` 保持两端存储顺序，
  不补造 identity/length。`R/data-relations.R:266–274` 只在方向已知时计算相对方向。
- **建议与缺证**：KEEP 已有明确语义，DEFER 缺陷定性；需真实类对象验证 metadata
  冲突、缺少长度、seqlevels、未知链及两端顺序，不能凭注释即宣布全部保真通过。
- **阶段**：v0.16.2；可选依赖保持可选。

## F21–F30：注册、输出、状态与专业能力

### F21 — P2/H：annotation registry 与标签算法共享边界

- **静态证据**：`R/layout-circular-registry.R:3–16,119,387` 有圆图资源与外侧注册表；
  `R/layout-collinear-annotation.R:409` 另有线性标签共享排布；
  `R/helpers-label-routing.R:151,649,893` 提供不同环形/弯曲标签处理。
- **建议与缺证**：DEFER 引擎合并；先比较算法输入域、单位和几何约束，再决定复用。
  registry 的一部分已在 export 帮助中公开（F07），须逐字段分级；候选布局状态
  不因可查看就成为新的稳定承诺。缺少跨模式等义布局与碰撞的运行比较。
- **阶段**：v0.16.3。

### F22 — P2/H：export、preview、diagnose 已有不同职责

- **静态证据**：`R/layout-export.R:64–68,197,231` 输出受选几何、render_context
  和 data 单位元数据；`R/viewer.R:175–190` 返回临时图片路径、接受目标尺寸；
  `R/diagnostics.R:153–167` 返回测量与建议，明确 measured 不是全部质量通过。
- **建议与缺证**：KEEP 职责分离，DEFER 返回字段冻结；核对三者在同一目标尺寸的
  单位、原始来源、对象不变性和估计/实测状态。不能把预览当交互图或诊断当无条件保证。
- **阶段**：v0.16.3。

### F23 — P3/D：标题范围小于实际能力

- **静态证据**：`DESCRIPTION:2,5–11` 标题仍只说 alignment chord，描述已包括
  collinear、single circular、feature 和 quantitative track。
- **建议与缺证**：DEFER 文案到 v0.16.5；比较保留包名并更新标题/首页定位与维持
  现标题只强化副标题两方案，不改变源码函数名。最终措辞由文档设计确认。

### F24 — P3/H：迁移错误分支是显式防误用机制

- **静态证据**：`R/zzz.R:81–92` 的 `ggchord_reject_retired()` 读取 dots 和原调用，
  以替代写法报错；`R/ggchord.R:114–125` 另拒绝未知参数；
  `R/geom-gene.R:416–450` 特别检查 R 部分匹配和已删除取值。
- **建议与缺证**：KEEP 有用错误方向，DEFER 清理；先列每项旧参数、已删除版本、
  替代和触发路径。不能以“兼容分支太多”删除会防止静默误配的保护。
- **阶段**：v0.16.0 建清单，v0.16.4 收敛；需代表性拼写/部分匹配运行证据。

### F25 — P2/H：内部无用函数尚无充分证据

- **静态证据**：`R/layout-sequence.R:2` quote 阶段、`R/aaa-geom.R:6–28`
  运行时 ggproto 回调、`R/theme.R:126–147` 动态函数以及 `R/zzz.R:7` 注册和
  `NAMESPACE:3–11` S3 均说明词法直接调用图不能覆盖全部使用。
- **建议与缺证**：DEFER REMOVE/INTERNALIZE；必须联合直接调用、对象字段、注册、
  字符串引用、延迟表达式、文档/示例及潜在下游。当前未认定任何内部函数“确定无用”。
- **阶段**：v0.16.0 盘点，v0.16.3 按证据处理。

### F26 — P2/H：有临时全局上下文，不能等同全局布局缓存

- **静态证据**：`R/render-context.R:3–18` 确有命名空间级 environment 保存当前
  context，`ggchord_with_render_context()` 用 on.exit 恢复旧值；`:21–25`
  为渲染副本创建新 ref 并 clone scales，`:95–136` 作目标设备测量与有限重建。
  `R/zzz.R:1–3` 的“没有全局 rendering state”注释比实际实现更绝对。
- **建议与缺证**：DEFER 泄漏定性；分别记录临时动态上下文与图对象持有的缓存。
  需嵌套、异常退出、多图、复制、尺寸变化和重复渲染验证，不能宣称“完全无全局状态”，
  也不能仅因这个 environment 就认定已发生串扰。
- **阶段**：v0.16.3；可先修正文档性绝对表述，行为改变需对应证据。

### F27 — P2/H：Stat 在布局前计算，再暴露 staged 值

- **静态证据**：`R/stat-ribbon.R:4–16` 的 compute_panel 透出私有美学中的计算列；
  `:88–115` 从 geom 创建层再挂统计和 private mapping，实际 bundle/density 数据
  在布局前准备。`:55–82` 的 density 保持源行并另附报告。
- **建议与缺证**：DEFER Stat 重构；验证 after_stat/after_scale、固定和映射优先级、
  grouping 及源坐标使用。结构不同于常见 Stat 不能直接证明 staged mapping 不兼容。
- **阶段**：v0.16.1。
- **有限运行证据**：本轮一次 chord 组合包含 point link、原始 ribbon、bundle 和
  density，手动 link_colour scale 保留；bundle 为两组且 bundle_n 为 1/2，density
  保留三组，after_stat(density) 基础映射得到有限值。不外推到 collinear、branch、
  全部计算字段、after_scale 或完整 staged mapping 契约。

### F28 — P2/H：dots 并非统一吞掉未知参数

- **静态证据**：`R/ggchord.R:114–125` 明确拒绝剩余 dots；
  `R/geom-gene.R:61–100,243–257` 做迁移检查后向 layer 透传；
  `R/geom-track.R:110–139` 经共享 layer helper 传 dots；theme 有独立通道。
- **建议与缺证**：DEFER 统一错误策略；逐入口记录合法透传目标和负责报错的层，
  检查拼写、匿名参数、重复参数及部分匹配，不能对构造器和标准 ggplot 层使用同一禁用策略。
- **阶段**：v0.16.1；已知迁移错误清单关联 F24。

### F29 — P2/H：比例尺图元与位置映射已是两种能力

- **静态证据**：`R/geom-seq-scale.R:34–68,76–113` 的 geom 是独立 bp 参考条，
  明确只支持 collinear；`R/scales.R:324–350` 返回 seq_position 连续 scale。
- **建议与缺证**：DEFER RENAME，优先比较文档任务分类/显著说明与新名称迁移两方案。
  需要新手使用或文档可发现性证据，不能把两者合并成同一 API 或伪装三坐标都支持。
- **阶段**：v0.16.1；v0.16.5 完成用户入口组织。

### F30 — P2/H：识别算法和来源字段已经存在，置信语义仍需审计

- **静态证据**：`R/common-features.R:321–329,676–706` 区分命中方法并返回来源与
  启发式 confidence；`R/primers.R:154–180,207–227` 区分 exact/three_prime，
  mismatch 固定零，TM 标明方法；`R/restriction-sites.R:54–63,142,234–256`
  保存来源行/版本，解释识别、切割与甲基化显示边界。
- **建议与缺证**：KEEP 搜索与绘图分开，DEFER 专业接口调整；核对数据库构建来源/
  许可、命中置信词、未知值、环形边界及实际规模性能。源码字段存在不能代替来源许可
  审阅，算法 high 不能理解成实验已验证或概率置信区间。
- **阶段**：v0.16.2；任何新增分析能力须另过真实需求准入。

## 公开入口与逐参数核查进展

已为四个入口的 29 个参数建立独立人工记录，放在 `audit-reviewed/` 并合并到生成
参数表。标记 SOURCE_REVIEWED，表示已阅读相关源码；不将 NOT_RUN 场景写成通过。

| 入口 | 参数行 | 当前建议 | 返回与错误边界 |
| --- | --- | --- | --- |
| ggchord | 6 | KEEP 构造职责与现有入口 | 返回 ggchord 图对象；validate=none 仍拒绝结构上不安全的输入；缓存 ref 不作为新公开字段承诺 |
| clean_ggchord_data | 9 | FIX F01/F02/F31，其他保留现状或 DEFER F03 | 返回 ggchord_clean 的三表和报告；不支持对象须先显式适配；keep 缺失 ID 不等于验证通过 |
| prepare_ggchord_data | 7 | KEEP 显式编排，F03 默认差异待取舍 | 返回三表、独立集合和 report；适配结构错误停止，生物学参照/边界问题以报告表达；不自动加图层 |
| validate_ggchord_data | 7 | KEEP 结构化检查职责 | 返回 ggchord_validation；strict 在收集问题后决定报错，不是逐行 fail-fast；四个逻辑开关的非法值立即报参数错误 |

记录逐项区分 NULL、NA、缺参、单位、嵌套参数、验证位置和透传路径。例如 prepare
顶层 mapping/steps 的 NULL 被拒绝，而 features/links 的 NULL 表示空集合；不能
为了“一致”在本轮直接统一。其余公开参数和内部参数仍须继续复核。

后续新增三坐标的 28 行、四种自有 Position 的 10 行，以及 feature/preset 的 data
两行，共 69 行 SOURCE_REVIEWED，涉及 13 个入口；这里的“涉及”不表示每个入口的
全部参数均已复核。feature 数据两行记录 F32 修复，其他 core geom 与 link/stat
参数的完整人工台账仍未完成，不能把有限运行例子当成逐参数审计已结束。

Position 的一次 collinear 公开导出确认：两条序列显示方向相反时，+/- 链 feature
仍分别位于各自骨架上/下。stack 的 legacy side 历史名称与 track 的法线 side
不是同一个契约，本轮只登记差异，不静默翻转几何或删除便利入口。

## 后续逐参数核查新增事项

### F31 — P2/C：unknown_id 的缺失/空 ID 合同冲突

- **证据基线**：`cf0dadd0da7b3125fd8a8ae14a0436958b837436`。
  `R/clean.R:22–24` 的公开参数说明将缺失、空和未知 ID 都纳入策略；实际在
  `:171–175`、`:379–381` 对缺失/空 ID 先报错，仅非空未知 ID 进入策略分支。
- **直接复现**：同一合法序列参考和 ribbon 输入，非空未知 subject ID 在
  `unknown_id = "drop"` 下被移除并有报告；改成 NA 则立即报 ID 缺失错误，
  不能按文档宣称该策略已覆盖缺失 ID。
- **用户决定与实现**：2026-10-30 用户明确选择 A，即为 ribbon/gene 的缺失/空 ID
  实现既有 drop/error/keep 策略并保留报告；B（保留前置错误并收窄文档）未采用。
  已于 `2c90371ce2a27ee4083c07a396150e78c46040b0` 实现，默认 drop 不变，缺失/空
  ID 的行为从前置报错改为移行并报告；error 显式停止，keep 保留且跳过无参考行。
- **实现证据**：新基线 `R/clean.R:148–153,188–219,391–422` 保留原始缺失掩码、
  报告文字值及源行；keep 同时检查原始缺失，防止 NaN 误配合法字面 ID `"NaN"`。
  两表的三策略、双端缺失、原输入保留和结构检查均已有局部直接断言，回归集中在
  单一 F31 案例中；未运行包级 testthat 或发布验收。
- **边界与阶段**：v0.16.0。seq_data 本身作为坐标参考仍必须具有非空、唯一的 ID；
  不猜测或生成 ID，keep 不保证可绘图。prepare 仍先严格适配，需要先独立清洗。

### F32 — P1/C：feature 的 factor 坐标被静默转为等级码

- **修复前基线**：`e44f1f76a3041e1ebe271b3dc30024bdcce4f50a` 的
  `R/geom-feature.R::ggchord_feature_data()` 直接 as.numeric(start/end)。公开
  geom_feature 的 collinear 布局导出证实：factor 标签 100/300 变成 1/2，并真的
  使用错误坐标绘制；不是仅函数名、类型风格或未使用参数的假设。
- **决定与修复**：FIX，提交 `4adab67f5ea680eb825a50f226dff119837e0923`。仅将
  顶层 factor start/end 先按数值标签转换，再交给既有展开/规范化路径；既有 numeric
  和 numeric-text 路径不变。修复恢复真实坐标，不新增函数或改变默认布局。
- **直接验证**：修复后公开 export 返回 start=100、end=300，原始 factor 输入保持
  不变。新增单一公开数据保留回归，不以私有 ggproto 字段或精确几何浮点坐标锁定实现。
- **文档与边界**：feature 和 plasmid preset 的 data 帮助、NEWS 与迁移已同步。
  子 segment 表自身仍遵循既有数值坐标合同；不声称所有 factor/复合区间格式均已验收。

## 待用户选择的跨模块方案

以下均为尚待确认方案，不产生立即改名、合并、删除或默认变化。F31 已选择并实现，见上节。
每项实际决定还需补上旧/新可运行调用、受影响用例、错误提示和删除版本。

| 事项 | 现状 | 方案 A | 方案 B | 当前建议与下一步 |
| --- | --- | --- | --- | --- |
| getter/export（F07） | 一个返回内部布局，一个受选输出，但 export 也公开说明了多种 registry 字段 | export 作为稳定推荐入口，getter 暂留并明确调试性质；收窄字段走迁移 | 两者都保留明确分层的公开契约，分别稳定必要字段 | 先列字段与用途，v0.16.3 决策；不能直接删除 getter |
| gene/feature（F09） | 有向基因与通用/未知链 feature，部分内部共享 | 保留语义便利入口，进一步共享内部实现 | 单一通用 geom 加显式角色，配完整迁移 | 优先比较真实生物学输入；尚不推荐公开合并 |
| 固定形状（F11） | shape 是 feature_shape 的简写，映射名仍 feature_shape | feature_shape 为首选，shape 保留或分阶段弃用 | shape 为固定参数首选，feature_shape 保留为映射名且迁移固定调用 | wrapper 透传错误已独立修复；选名不与该修复绑定 |
| 标签命名（F12） | label_* 与 gene_label_* 并存，默认方向/布局不同 | 各 geom 统一短 label_*，保留角色美学，迁移旧参数 | 保留角色前缀，提供清晰横向说明并共享实现 | 比较默认与跨坐标后在 v0.16.1 选择，不能仅文本替换 |
| color/colour（F14） | 五个直接别名 | colour 为文档首选，color 长期作为便利别名 | colour 唯一推荐，color 经警告与迁移后删除 | 两方案均可，需使用负担证据；不预设“别名无用” |
| 重导出（F15） | 22 个便利符号 | 保留便利集合，reference 明确来源 | 仅保留入门必需项，其他改显式 ggplot2/grid 调用 | 逐项建立用途、来源和迁移，不做数量驱动裁剪 |
| clean 默认（F03） | 独立 sort，pipeline keep 且有方向保护 | 保留任务相关默认并强化文档/报告 | 统一显式策略或统一默认，给出现有调用影响 | 先验证方向依赖，不静默统一为 sort |
| file/files（F16） | 共享 resolver、允许同时提供，两路径解析不对称 | 保留兼容双入口并统一解析/来源规则 | 选择一个支持路径向量/glob 的入口并迁移另一项 | 先最小验证实际差异，v0.16.2 选择接口 |
| offset 删除（F10） | 帮助原承诺 v1.0.0 删除，Position 已可替代 | 保留至原定删除点，在此前完成迁移 | 在候选冻结前删除，明确修订原承诺并给迁移版本 | 具体删除时间仍待确认，不能把总体路线当作单项删除授权 |

## 本轮交付的边界

已逐项定位 F01–F30 的相关源码事实，纠正“全部同类图层一定覆盖”“所有 registry
都从未公开”“track 没有冲突检查”“完全没有全局渲染上下文”等可能的过度概括。
这不意味着已动态确认所有假设或已完成全部函数/参数审计。

下一步继续以机读全量清单补上其余入口和参数证据，优先解决影响实施安全的未知；
F01/F02/F11/F31/F32 的局部复现与修复已单独取得直接验证，其他项仍按各卡片缺证范围推进。
未运行包级 testthat、打包/check、vignette
或设备矩阵，未根据本表执行批量 API 改名删除。
