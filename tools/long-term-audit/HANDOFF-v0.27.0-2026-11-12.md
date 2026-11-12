# ggchord 新会话交接说明（2026-11-12）

本文件用于恢复普通软件工程开发，不是完成验收记录，也不保证新会话能解除服务端限制。
保留全部来源约束及拦截事实；不得通过隐瞒背景、伪造来源、改名或关闭防护规避安全检查。

## 用户目标与授权边界

- v0.26.0 已获授权整合本地 main；从整合后的 main 创建 v0.27.0，继续完成 v10 任务书其余开发。
- 按稳定阶段、单一范围及时本地提交；不要每完成一个工作包就停止。全部任务尚未完成。
- 后续样例使用**不涉及高风险病原体、来源可追溯的真实数据**。格式控制用例可以保留，但须明确其用途，不能冒充真实业务数据。
- 现有主题系统讨论纲要是设计输入，不是预先决定的 API 实施清单；须结合源码审计和真实用例论证。
- 继续使用原 OneDrive checkout；不创建 worktree、不复制仓库。不自动 push、部署、发布，未经新授权不把 v0.27.0 合入 main。
- 普通开发仅做直接相关的轻量验证。未获正式 GitHub 发布授权，不运行 testthat、devtools::test()、R CMD build/check、vignette 编译和完整发布设备矩阵。
- 冻结 docs/ 不重建；不覆盖用户无关改动。沟通使用汉语。

## 当前 Git 与版本快照（本轮实际核对）

仓库：`/Users/vutmvu/Library/CloudStorage/OneDrive-个人/R语言/ggchord`

| 项目 | 当前值 |
| --- | --- |
| 当前分支 | v0.27.0 |
| HEAD | 6c8eca594b4ffceb3681cd70619e649b8440efc9 |
| 本地 main | 56f76af3c5d1bdde8f0bd74f235c70857abbfa2d |
| DESCRIPTION Version | 0.27.0.9000 |
| 最近正式版（路线图记录） | v0.16.6 |
| 当前版状态 | 开发中，未整合 main、未发布 |

本轮未 fetch，不声称已核实远端当前状态。路线图和 V10-STATUS 含不同历史时点的 main/origin/main 表述，
不得用旧摘要覆盖以上实测本地引用，亦不得把本地 origin/main 引用当远端实时事实。

## 开始工作时读取

1. 根目录 AGENTS.md。
2. design-roadmap/DESIGN-ROADMAP.md。
3. design-roadmap/versions/v0.27.0.md **全文**，再读相关实现。
4. tools/long-term-audit/V10-STATUS.md、DEVELOPMENT-v0.27.0.md 和 THEME-DESIGN-v0.27.0.md。
5. 原始任务书：`/Users/vutmvu/Downloads/ggchord_长期开发路线图与统一架构_Codex任务书_v10_独立执行版.md`。

不能只靠本交接说明或状态表开发。V10-STATUS 后半部分保留 v0.25/v0.26 当时判断；M2、
部分 M3/Fusion 及适配器已有后续实现，需按最新记录与实际代码逐项核对，避免重复开发或误关任务。
`design-roadmap/` 当前被 Git 忽略，仍须维护本地权威设计文档，不强制加入版本控制。

## 已提交工作与证据边界

以下全部已在当前分支提交；“完成”限于所述实现与定向验证，不等于完整任务书、全部输入或发布验收完成。

| 提交 | 工作 |
| --- | --- |
| 6340485 | 原生 theme 元素、Track 数值轴与局部 guide 外观，主题职责审计；保留原 109 形参，未任意删 API |
| cf262ca | Fusion 八方位/角度/显式 box、Sankey/Radial 注册、嵌套 magnify、真实来源命名空间与重复展示许可 |
| e933efb | 显式 circular 来源与跨原点选区，保留精确源位置；不从圆形 Coord 推断生物学环状 |
| 6f47bf4 | 实测毫米障碍避让和显式共享显示通道；不合并源观测，失败明确报告，不宣称全局最优 |
| a1a7b79 | Fusion 外层原生多侧 Guide 联合预算、用户锁定及最终尺寸报告 |
| 6d19945 | 显式 side="auto" 外围环绕、有界确定性排布及不可容纳诊断 |
| 96ad72d | 原生组合标题/副标题/图注/tag、背景/边距与毫米转换；修复嵌套继承裁剪/透明绘制 |
| f55772f | 修复普通子图边距造成的无依据连接断带，保留前景和后层覆盖语义 |
| f6056e7 | 闭合 Track 接缝周期法线修复，真正越过圆心的折叠仍检测 |
| c2e5105 | read_bedgraph/read_bigwig；坐标、缺失、查询与来源边界；修复 BED 类合法 track/browser ID 被误丢 |
| 4c5d2da | read_vcf、现代 BCF 可选成熟解码、显式 SV 断点桥接；保留 INFO/基因型/置信区间，不推断支持度或 mate |
| 01e0804 | read_maf/read_chain 成熟解析桥接；真实字符/负链/gap/块来源，chain 不伪造字符或 identity |
| f235a53 | read_delta 成熟 show-coords 桥接；原始方向/长度/脚本/原生统计，不运行 aligner |
| 79d8530 | 共用可选解析进程边界，保留 venv 解释器入口，避免 symlink 解析后丢环境 |
| 6c8eca5 | read_signature_matches 离线 TSV/缓存 JSON；原生 AA segments、源 score 类型和完整对象，不填 gap 或增加命中数 |

相关 NEWS、roxygen/Rd、双语课程和 API 清单已随对应阶段更新。具体脚本、输入、结果与限制见
DEVELOPMENT-v0.27.0.md；代表图已按各阶段记录实际生成/查看，不能据此宣称发布设备矩阵通过。
最近记录的 formats 独立课程执行通过、reference 翻译检查 1134 段通过；本次交接没有重跑这些检查。

主题后续仍遵循：Theme 管非数据外观；Scale 管数据映射；Guide 管显示组织；Geom 管固定样式和形态；
Position/Coord/Data/Stat 各自承担布局/投影/处理。显式用户配置优先，同 Track ID 与不同 ID 的共享/隔离语义需保持。
预设切换不改变原始数据、序列方向、筛选或真实连接，不为每条 Track 增加一整套冗余参数。

## 交接前的未提交改动

实际核对时只有以下六个未跟踪文件，无已跟踪文件修改；本交接文件及本轮本地版本记录补充随后新增。
不要清理这些资产，也不要直接视为已完成的 reader。

| 文件（均在 inst/extdata/） | 内容 |
| --- | --- |
| FEATURE-RELATIONS-LICENSE.txt | 上游原始 BSD-2-Clause 许可 |
| feature-relation-sources.json | 固定源提交、原文件 URL/SHA256、抽样规则、原物理行号和产物哈希 |
| plant-feature-pairs.homology | 前 20 条完整原始植物关系记录，未改 ID/坐标/score |
| plant-pair-features.gff | 匹配上述端点 ID 的 28 条完整原始位置记录 |
| plant-feature-groups.txt | 3 个完整原生命名组，每组 4 个成员；原文件行 1167–1169 |
| plant-group-features.gff | 匹配上述成员的 12 条原始位置记录 |

来源为作者公开 MCScanX 仓库，固定提交 `0956cb8f900c152e2be8ba3196829e50ab454b94`。
完整 URL 和哈希以来源 JSON 为准；第三方名称可出现在来源文档，不写进新增函数/参数名。
该四列 .gff 是该工具的位置表，**不是标准 GFF3**，也没有完整染色体长度。
成对 ID 带 at0| 等前缀，组文件 ID 不带此前缀，各有匹配位置表；不能自动删前缀合并命名空间。
第三列 score 是原始关系分值，不自动称 evalue、identity、置信度或支持度。

临时准备脚本 `/tmp/ggchord-prepare-relations.py` 尚未版本化。原始下载位于 /tmp/ggchord-plant.homology、
/tmp/ggchord-plant.gff、/tmp/ggchord-plant.groups、/tmp/ggchord-plant-groups.gff；临时文件可能消失，不能成为唯一来源。
下一阶段须提供 data-raw 中可重复的生成规则，校验来源哈希并保留完整原始记录抽样依据。

## 紧接着的开发阶段（尚未实现）

第十五工作包：预计算 feature 对及同源组。当前只有真实输入调研/资产准备，没有 reader、帮助或验证闭环。
候选接口 read_feature_relations 仍待原型论证，不是已冻结承诺。

- 显式区分 pairs、named_groups、无名 groups；数字状原生 group ID 保留字符串，无名组仅有源组索引。
- 成员关系返回 membership，不能自动扩展为全体两两连接或推断同源/直系同源关系。
- 可选 feature_data 按唯一 feature_id 精确映射真实区间；不猜 gene/CDS 父子关系或修改版本/大小写。
- 歧义须明确失败；缺失默认错误，可显式 report 保留原记录。原始表、物理行号、文件范围与 score 保留。
- feature strand 与 alignment orientation 不同，不把前者当后者；没有方向证据时保留未知。
- 不从最大 feature 端点推断整条序列长度，不生成伪完整 seq_data，不声称已支持未经真实输出验证的 .collinearity。
- 完成原型/实现、来源生成规则、直接轻量核对、文档/NEWS/Rd/双语课程和 API 同步后再单独提交。

## 其余实质任务（必须回查任务书逐项闭环）

1. 适配器：完成上述关系/组；BAM/CRAM 按需、GFA/ODGI 探索性范围需有真实调研/桥接证据与明确取舍，不造空壳 reader。
2. 样例治理：后续主动选用符合新约束的真实来源。delta 历史 NUCMER 验证含不符合新约束的病原体样例，
   不再将其作为后续主动开发样例；评估用可追溯、符合约束的真实输入替换活动验证依赖，保留历史来源事实，不能只改名冒充新数据。
3. M3/Fusion F2/F3：已具备显式放置、嵌套、路由、多侧 Guide 和自动环绕，仍须补高密度/复杂组合质量证据及适用边界，不因有限定向通过而整体清零。
4. 制图质量：密集标签、引线穿字/框、裁切等指标；真实同输入同尺寸公平对照；Adaptive 更广视图与宽高比、显式约束及可解释失败。
5. Track 范围：已有 14 种几何及 identity/count/bin/summary 等；facets、smooth、boxplot/violin 等仍有限制。
   按 v10 确认实际要求，不能默认为任意第三方 Geom/Stat 全透传，不能擅自把已约定核心项延期。
6. API/兼容：集中职责、冗余与迁移审计；109 个主题形参目前保留，不以精简为由任意删除；破坏性变化须可执行迁移与 NEWS。
7. 文档：旧 lesson engine 和 178 张旧图实际重跑，不能只改指纹/stamp；真实主示例来源、双语页面配对与源站响应式检查。
   日常不得重建冻结 docs，也不得把课程运行自动升级为包级发布流程；若具体任务与规则冲突，先指出冲突。
8. 稳定化：按完整任务书正文/附录和最新实现建立真实完成证据；正式发布检查与远端 CI 等留待明确授权。

## 可选开发环境（临时，不保证新会话仍存在）

- R 4.6/macOS arm64；可选包库 `/tmp/ggchord-adapter-r-lib`，含 rtracklayer、VariantAnnotation、StructuralVariantAnnotation 等。
  定向脚本可显式设置 .libPaths，不能修改用户默认库或给运行时强加整个 tidyverse。
- `/tmp/ggchord-hts-python-venv/bin/python`：已记录 pysam 0.24.1、Biopython 1.85、numpy。
  保留解释器入口 symlink 路径，不能 normalize 成另一个 Python。
- `/tmp/ggchord-show-coords`：仅编译成熟解析组件，构建源 v4.0.0rc1；没有运行比对器。
- 既有脚本环境入口 BCF_PYTHON、PARSER_PYTHON、ALIGNMENT_PYTHON、DELTA_DECODER，按脚本实际参数核对。
- 文档源工具：generate-parameter-atlas.py、generate-api-curriculum.py、sync-tutorials.py、reference-translations.py。
  只在相关修改时运行；保持清单顺序，避免整文件无关排序和 roxygen 无关空白变动。

## 服务端拦截记录与未知事项

原会话 ID：`01a1ce55-c45c-7061-a6da-e63c39278c40`。用户可见提示为“无法显示此内容”，并提到生物研究风险及可信访问。
2026-11-12 用应用 read_thread 核实以下失败记录：

- 较早失败 turn：`01a1cf73-9871-7cf2-9afa-2dcd9ccdf561`。
- 最近连续两次失败 turn：`01a1cfcc-e46f-7dc0-b2ab-11113954f533`、`01a1cfd0-aa06-7a83-9a4e-3a90bc22386e`。
- 状态均 failed；错误开头原文：`Invalid prompt: we've limited access to this content for safety reasons.`
- 服务端附带说明涉及信息可能有益或有害；additionalDetails 为 null，没有给出具体分类条件。

这证明请求层发生安全限制，不是已观测到的 R/Git 执行错误，也不是工具自动审批拒绝。
不能确定单个关键词、某个样例或累积上下文是唯一原因；不能推断账户被封禁。
此前仅归因于物种名的说法证据不足，已向用户更正。最近这两次失败没有产生新的开发成果。
新会话仅用于准确交接和合法工程延续；若继续遇到限制，应保留原始错误、明确报告，不能反复改写隐藏背景绕过。
本文件没有修改任何安全设置，没有把限制写成 ggchord 的功能要求，也没有宣称问题已解决。

## 建议恢复顺序

先核对本地状态与上述快照，读完规定文档，保留这六个资产；再完成第十五工作包并按单一范围提交。
随后按依赖推进其余任务，每阶段同步实现、来源、轻量验证、当前版本设计与面向用户文档。
有冲突或关键输入缺失时具体说明，同时继续独立工作；不要在单个工作包后把全任务停止或标成全部完成。
