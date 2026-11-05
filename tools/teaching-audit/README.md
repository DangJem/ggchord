# 三模块教程重构审计

版本 **0.16.7.9000**，分支 **v0.16.7**。首轮审计基线 `43967187ca81c7e531d4bb9182c29b9563d2ab84`；代码修复提交 `fcf5f4e`，当前矩阵读取该提交之上的工作区；父 SHA 不表示未提交内容已属于该提交。具体文件 SHA256 见 `teaching-summary-v0.16.7.json`，本地站构建清单另保留源码哈希。此前阶段的待办和失败记录移至 [历史记录](EARLIER-STAGES.md)，不能当作本次最终状态。

## 七项交付

| 交付 | 文件及范围 |
| --- | --- |
| API/参数矩阵 | `teaching-api-v0.16.7.csv` 覆盖 130 导出，108 包自身入口、22 上游/对象；`teaching-parameters-v0.16.7.csv` 包含 1474 formal 的真实默认、类型规则、约束、模式、优先级、等级及证据 |
| 缺口及完成证据 | 下表与 `completion-checklist.json`；状态依据正文、执行和图像，不按函数提及清零 |
| 课程与旧页映射 | 下表；`website/navigation.json`、`pages.json` 为最终目录；原页面/锚点保留 |
| 双语源码 | `website/en/`、`website/zh/`，共 59 对页面；新增主讲课而非复制 Rd；轻量离线教程只同步源，不重建 |
| 对照目录 | `website/parameter-comparisons.json` 共 71 个对照输出，记录输入、基线、参数/默认、尺寸、依赖和状态；169 张现版本 PNG 均实际重新导出 |
| 验证 | `lesson-execution.json`、`version-check.json`、`layer-compatibility.csv`、`verification.json` |
| 变更与边界 | 本文件、NEWS.md；两项确切渲染错误已修复，无破坏性签名或依赖变更；未整合、推送、部署、发布 |

## P0/P1/P2 缺口处理与最终课程

| 级别 | 原缺口 | 最终实现与证据 |
| --- | --- | --- |
| P0 | 正式版安装与开发示例混用 | installation 明示本地源码版本、全新会话、正式版对应文档；三个首课各在 0.16.6/开发版独立进程渲染，共六次；实际探测新参数/蛋白对象 |
| P0 | 跨页准备与可选导入结果不明 | 三类首课显式 library；47 页按顺序执行准备和原生 R 构图；11 个新专题/完整案例另在新进程执行；Biopython 1.85 实际读取 L09137.2，不将跳过写成成功 |
| P0 | 模式按名字误判、帮助矛盾 | 21 geom × 3 坐标共 63 组，53 成功/10 按契约拒绝；compatibility 与对应 Rd 澄清；原始合成几何 fixture 不冒充生物数据 |
| P1 | 共线性布局/方向未成课 | collinear-layout：left/center/right、anchor、manual/auto、min_support、四种方向及冲突；views/homology 承接窗口、分段、源坐标 |
| P1 | 箭头、分支、点连接不足 | feature-geometry：真实 integrase 三种头部与短特征回退；chord-branches：显式标注共享端点示意、分支长度和直/曲点连接 |
| P1 | 圆图中心、骨架、注释混于大案例 | circular-layers：单/双/band、线间距/带宽、方向/开口、中心独立样式、通用形状和 Position；band 真实代码缺陷已修复 |
| P1 | 固定/自动标签欠缺参数对照 | circular-labels：方向、内外、旋转、feature/radial/auto/callout、换行/缩小/省略与外圈回退；labels/chord-labels 讲基因/序列及诊断 |
| P1 | 引物/酶切搜索与绘制混淆 | primers、restriction-sites：目录、真实搜索表、弧线与文字、内容/排序/侧向、引线与甲基化限制分别主讲 |
| P1 | 定量轨道缺 tile/组合 | tracks：真实 GC 窗口的 bar/line/tile、4/10 mm 色带、两条独立轨道/次序、边界裁剪、分母与模糊字符规则 |
| P1 | API/参数只有名称盘点 | 全部 130 条显式语义判断；parameter-atlas 1474 条参考合同；api-curriculum 提供用途→课→模式→参数→案例映射；task-index 40 项任务 |
| P2 | scale/guide/theme/统计辅助入口不全 | scale-families、ribbon-statistics、workflow-tools；所有非别名包内入口在展开课程有调用，拼写别名和上游明确分类；密度/权重/成员和预处理精确定义 |
| P2 | 完整案例只有大代码 | chord、collinear、plasmids、real-data 补研究问题、来源权利、筛选/变换、选择理由、尺寸与推断限制；case-reviews.json 按完整案例评 L4；13 图库为可选参考 |
| P2 | 旧链接、搜索与移动阅读 | source/full documentation checker 检查双语代码、875 译文、图像依赖、旧 URL、锚点、搜索；浏览器结果见 verification.json |

基础必学：三个 *-start → 对应 data/build → 基本 tuning → 对应完整案例。进阶可选：几何、布局、标签、分支和局部窗口。高级专题：投影、统计、多轨道、复杂引物/酶切与诊断。原 A01–A16、B01–B15、C01–C19、D01–D12 由这些职责承接，不机械新建每个单元：

| 单元 | 主讲位置 |
| --- | --- |
| A01–A07 | chord-start/data/build/tuning/geometry，feature-geometry，modes#coordinate-controls |
| A08–A14 | chord-branches/routing/labels，ribbons，ribbon-statistics，scale-families |
| A15–A16 | chord-import，chord，real-data |
| B01–B07 | collinear-start/data/build/tuning/layout，feature-geometry |
| B08–B12 | collinear-routing，views，homology，labels |
| B13–B15 | tracks，formats，collinear，real-data |
| C01–C08 | plasmid-start/data/build/tuning，circular-layers，annotations |
| C09–C14 | circular-labels，primers，restriction-sites，recognition |
| C15–C19 | tracks，formats，plasmids；plasmid-catalog 可选 |
| D01–D06 | modes，adjustments，styling，scale-families |
| D07–D12 | labels/circular-labels，annotations，tracks，views，output，faq，workflow-tools，compatibility |
| 任务/API 检索 | task-index，api-curriculum，parameter-atlas，reference |

## 等级和证据的解释

130 条 API 评为 L3 28 条、L2 80 条、上游/对象 L1 22 条。五个纯拼写别名复用 canonical 教学，不重复创建例子。L3 只指登记的关键参数行为，不代表每个 formal、坐标组合和设备均有穷举对照。简单断点图层以 L2 最小公开使用讲授，完整签名归参考。四个综合案例单独评 L4，不把案例中每个函数的所有参数自动抬到 L4。

参数的静态出现最多证明参考覆盖；人工决定见 `website/teaching-reviews.json`，执行见独立记录。默认表达式来自 `formals()`；NULL 可能在布局阶段解析，候选向量不总是单选，类型与范围须结合合同说明。`same_object_as` 仅记录 R 对象 identical，不能据此把 after_stat/after_scale 当作语义别名。`alias` 只记录源码显式赋值。

API 与参数矩阵、运行日志和图片并非互相替代。source checker 对每个 canonical 入口要求实际课内调用和成功执行页，同时仍不自动决定 L2/L3。

## 验证方式与使用边界

- `PATH=/tmp/ggchord-doc-python/bin:$PATH Rscript tools/documentation-plots.R`：169 PNG 实际导出；使用 view_ggchord、同一登记尺寸和输入，最终无连续点连接颜色条警告。
- `Rscript tools/check-documentation-lessons.R --write-outputs`：47 页分离环境；`--page=slug` 可在独立 Rscript 中复核专题。直接 R 代码也执行；预览/保存只核对对象与参数，不覆盖用户文件，实际导出另有上述证据。两条既有 Seqinfo 警告来自区间适配。
- `python3 tools/check-first-lessons.py`：三首课各自全新进程；`check-documentation-versions.py` 另对真实安装的 CRAN 0.16.6/开发源码作六次渲染和能力探测，保留正式源包 SHA256。
- `Rscript tools/check-layer-compatibility.R`：63 组公开图层渲染，条件拒绝属于预期。不是全部参数组合穷举。
- `Rscript tools/check-teaching-defaults.R`：108 个包内函数的 Rd usage 与实际参数、默认和顺序一致。
- `python3 tools/check-documentation.py`、`generate-parameter-atlas.py --check`、`generate-api-curriculum.py --check`、`reference-translations.py --check` 与 `sync-tutorials.py --check`：源和生成结果一致。
- `python3 tools/build-website.py`：只生成本地 `.site`，不改冻结 docs，不部署；完整站点检查及桌面/手机实际结果见 verification.json。

本轮未运行 testthat、R CMD build/check、vignette knit、PNG/PDF/SVG 发布设备矩阵；仓库明确要求只在获正式 GitHub 发布授权后执行，当前未授权。这些是发布边界，不是本轮重构遗留。一般设备/字体移植、识别结果的实验效力及全部低频参数组合不在本次验证承诺内。未来独立 GC/覆盖度 stat 等未确认候选没有被写成已实现功能。README 已增加开源设计参考并避免商业工具名称，来源追溯记录仍保留。
