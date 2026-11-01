# v0.16.6 RC 本地开发完成记录

2026-11-02（Asia/Shanghai）。**前阶段六个工作包已完成；其本地main整合随后按用户
要求撤回，当前继续在 `codex/v0.16.6` 开发。没有正式发布、push、部署网站或创建
v1.0.0分支。** 本记录保留前阶段证据，新增要求见版本文档。开发包标识 `0.16.6.9000`，正式版及安装
tag仍为v0.16.5。[合同审计](API-RC-CONTRACTS-v0.16.6.md)与本版rc表记录具体边界。

## 六个工作包

| 工作包 | 实际交付 | 产物 |
| --- | --- | --- |
| 1 当前覆盖 | 130exports、1,236自有形参、49美学角色、9S3、130返回说明、200迁移项均有关联合同；当前签名无漂移 | rc-exports/parameters/aesthetics/s3/returns/migrations/source-hashes；tools/audit-rc.R |
| 2 契约取舍 | getter/export、gene/feature/plasmid、点/区间link、5别名、22重导出保留理由；27条灵活参数边界；默认/来源/返回字段复核 | API-RC-CONTRACTS；rc-parameter-formats；rc-data-fields（34入口）、rc-engine-fields（16组字段） |
| 3 行为修复 | 六组有限行为；F54内侧引物不再被移到外侧；F55内侧酶切registry纠正；两个最小回归体直接验证 | rc-behavior/regressions；R/ggchord-build.R、layout-external-annotation.R、layout-perimeter.R、layout-circular-registry.R；现有两份回归文件 |
| 4 支持基准 | 声明/实际/未知分开；四种真实输入的小型build基准；两个PNG场景与显式场景清单 | rc-support/benchmark/scenes/compatibility/tools-review；tools/benchmark-rc.R、render-rc.R、rc-fixtures.R、audit-rc-support.R |
| 5 文档候选 | NEWS、相关roxygen/Rd及汉语参考摘要、双语参数/安装/质粒/迁移/版本页、双语Rmd源同步 | NEWS、R/geom-primer.R、R/layout-export.R、man、website、MIGRATION两份、vignettes两份 |
| 6 完成核对 | 前阶段问题/限制、提交和验证对应；曾本地fast-forward main，用户随后撤回该整合；全部成果留在版本分支 | 本文件、versions/v0.16.6.md、DESIGN-ROADMAP；main恢复d697345，版本分支继续开发 |

当前帮助相对候选有149处参数文案、11处返回文案差异，包含本版侧别修复的澄清；
逐项保留，不冒称149个新行为。69个原支持合同中，17个geom/coord对被本次有限场景
实际覆盖，其他继续NOT_RUN_IN_RC。历史表没有改写，NOT_RUN没有批量升级为PASS。

## 修复、保留与限制

- F54：真实pUC19单引物及六个universal引物均复现inside被移到outside。构建局部
  几何现在保留请求侧别，外侧任务/slot不再接收内侧引物，内侧登记实测框；
  outside默认、锚点和公开签名不变。新增一个最小公开回归。
- F55：inside restriction几何原本在内侧，registry却标outer；现登记inner，保留
  原内侧几何，不分配exterior slot。复用酶切案例增加一个侧别断言。
- F52/F53、稀疏混合排序、@export及ASCII帮助显示副本为v0.16.5已发布事实，合同
  审计已逐项追踪，不冒称本版首次实现；首发/重发标签区别保留。
- 所有既有有效向量/列表形式、5个拼写别名和22个重导出保留，推荐命名向量；
  没有新增弃用、删除、迁移期限或需用户选择的破坏性方案。
- RC-R04：完整显式inside圆图中部分长引物文字与feature区域交叠，侧别/内容正确。
  默认max_overlaps=Inf不提供跨feature零交叠保证；默认outside更有空间，缩短显示
  文字后仍需看最终尺寸。此为已记录的非阻断视觉限制，inside未标为全视觉PASS。
  进一步全局内侧避障是待确认优化，不承诺下版，也不作为本版已实现功能。
- 原生aa坐标仍另案。CDS投影是显式坐标转换；同源例子核对真实CDS DNA，显示连接
  不冒充新比对。源strand、显示orientation和link direction保持分开。

## 实际验证

本机R4.6.0、ggplot2 4.0.3、arm64 macOS（Darwin27）。仅运行针对性开发检查；
没有完整testthat/devtools::test、R CMD build/check、vignette重建、发布设备矩阵或
远端CI。两个新增回归体使用直接断言执行，没有加载testthat运行。

| 命令 | 结果和边界 |
| --- | --- |
| `Rscript tools/audit-rc.R` | 当前覆盖和源码字节指纹生成；不是逐参数执行；历史candidate/engine/data/grammar CSV未改 |
| `Rscript tools/freeze-api.R check` | 130候选导出签名、5别名、22上游identity、9S3保持；未运行snapshot |
| `Rscript tools/validate-rc.R` | layers/circular/linear/labels/output/data六项PASS_FINITE；修正报告闭包读取初始结果的问题后，由实际result写文件，没有手动改PASS |
| `Rscript tools/validate-rc-formats.R` | 一个真实数据命名半径重排、嵌套链width/offset组合通过；未穷举等价形式 |
| `Rscript -`（新增F54/F55回归体直接断言） | 默认主题下inside引物源行/侧别、inside酶切资源数量/侧别通过，见rc-regressions；长期案例在现有测试文件 |
| `Rscript tools/audit-rc-support.R` | 声明、安装版本与未知组合分开，安装存在不冒充行为通过 |
| `Rscript tools/benchmark-rc.R 1` | 四种真实输入build各一次，0.251–1.651秒；对象约28.5–29.0MB，非峰值内存，无跨机器门槛/加速声明 |
| `Rscript tools/render-rc.R circular` | 9×8in、120dpi、sans PNG实际渲染并检查；完整feature+6primers+3cuts，默认外侧内容/引线可读 |
| `Rscript tools/render-rc.R inside` | 同尺寸PNG实际渲染；内容/侧别保留，RC-R04交叠限制如实登记 |
| `python3 tools/sync-tutorials.py` | 只同步Markdown/Rmd源，不knit/建站；README两份无差异，未恢复reallinear展示 |
| `python3 tools/check-documentation.py --source-only` | 27语言对、共享代码/图片、130exports/93Rd通过；865段汉语参考及源摘要对应 |
| `git diff --check` | 通过；历史CSV、冻结docs、README两份、原始输入及许可证保持无差异 |

真实HSP有反向端点，constructor出现1类预期validation warning；未排序/篡改输入，
也未将它当作设备失败或隐藏未知警告。图片、manifest和探查脚本在忽略目录
`.codex-tmp/rc-v0.16.6/`，PNG哈希/视觉状态另存rc-scenes，可显式重现。

## 支持限制与待授权发布检查

- 保留R>=4.1.0、ggplot2>=4.0.0声明；下限组合本轮NOT_RUN。Windows/Linux/其他
  macOS和其他字体未在本版执行。v0.16.5四平台及27设备文件只作历史证据。
- PNG仅两个场景、两个目标尺寸和一个preview；PDF仅私有度量设备，SVG/Cairo/
  Quartz输出未验收。测量成功不认证嵌入字体、所有碰撞或跨软件效果。
- XLSX/Bioconductor/GenBank外部parser及完整文档构建未执行。已安装Suggests标为
  INSTALLED_NOT_EXERCISED，缺失依赖不虚构运行。安装情况见rc-support。
- 密集场景最多40个真实非冗余单切位点，可显式选择；本轮未渲染该场景、未做
  极端压力或全排列。性能为单次小型观察，不作跨机器或旧版比较。
- 用户以后明确授权本版GitHub正式发布时，才执行完整testthat、标准build/check、
  英语/汉语vignette、PNG/PDF/SVG及跨平台验收，按同一最终候选准备双语Release、
  tag和资产校验。本地完成不等于发布门槛已通过。

## 提交与整合

- `f6d4b49`：来源会话启动计划审计，本轮起点。
- `409e216`：当前合同覆盖台账与开发版本标识。
- `b4c57aa`：F54/F55修复、单一必要回归及帮助。
- `28af660`：合同取舍、有限行为、支持与基准工具/证据。
- 最终文档提交：本记录、路线图、NEWS、相关双语源和同步的教程/迁移。

原OneDrive checkout直接操作，无worktree、仓库复制或子代理。全部实现、测试源及
文档在版本分支完成后，曾按旧规则本地fast-forward main至 `5c83015`；用户随后明确
取消整合。切回 `codex/v0.16.6` 后，main指针恢复至 `d697345`，版本分支保留全部
开发提交。新的长期规则是不经用户明确批准不整合main。没有push、Release、部署
或新版本分支。

2026-11-02追加工作：WORDLIST已更新并用显式工具作针对性文档拼写审阅；7个官方
绘图包的测试组织调研已完成，见 [记录](TESTING-RESEARCH-v0.16.6.md)。测试精简仅为
待确认草案，当前tests未删改，既有证据仍保留。当前在版本分支继续开发，不自动整合main。
