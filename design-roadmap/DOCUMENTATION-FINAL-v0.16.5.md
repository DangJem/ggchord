# v0.16.5 第二轮实现记录（新增优化要求尚未完成）

> 用户随后要求继续优化成品质量、真实精简输入与格式说明。本记录仅保留第二轮有限验证；
> 当前全部完成条件以版本文档末尾补充为准，不作为新的完成结论。

日期：2026-11-01。版本分支 `codex/v0.16.5`。
状态：实现、共享例子、站点静态检查和浏览器有限检查已完成；等待最终清洁源码构建确认。
本记录取代首轮的完成结论，不删除其历史验证事实。当前正式版仍为 v0.15.0。

## 实施结果

唯一正式生成基础为 pkgdown 2.2.1 + Bootstrap 5；使用其原生模板和完整 reference，
保留成对 Markdown / 共享 R 代码及图片。首轮自制页面壳和 Rd 转换器退役。
25 对内容入口、93 个 Rd 主题均有英语/汉语导航；主要教程和图注同内容成对，
reference 保留完整权威英语正文，汉语任务说明与原文边界明确。用户原画适量放在
首页、README 和教程入口，RGBA 无损编码记录保留。

18 个主要功能簇按输入、代码、图/结果和解释逐项复核；数量只用于定位，不作为
深度的替代证据。共享源码调用 87 个自有公开函数与 9 个上游入口，其余多为别名、
标度变体或共享解释，不为每种等价形式机械复制教程。

| 功能簇 | 文章 / 共享片段 / 图 | 内容复核 |
| --- | --- | --- |
| 三种坐标模式 | [modes](../website/zh/modes.md)；`chord`、`circular`、`collinear`；chord、circular、collinear | 三种坐标、anchor/自动方向/邻接规则与源坐标边界。 |
| 真实记录比较 | [real-data](../website/zh/real-data.md)；`realimport`、`realprepare`、`realchord`、`reallinear`；realchord、reallinear | 4 条原始序列、266 条注释；8 条筛选连接；弦图明确选择 24 条 CDS / 8 标签，线性图保留全部注释。 |
| 完整圆图注释 | [plasmids](../website/zh/plasmids.md)；`plasmid`；plasmid | pUC19 2686 bp、9 feature、6 universal 精确引物；全量 49845 酶匹配中明确选择 3 个酶位点。 |
| 校验、清理与层次 | [data](../website/zh/data.md)；`validation`、`cleaning`、`hierarchy`；结果/报告 | 错误合成输入 2 行 → 显式修复后 1 行；保留原表及报告，声明 genome/contig 外键。 |
| 文件格式与适配器 | [formats](../website/zh/formats.md)；`adapters`、`formats`、`genbank`、`optionalranges`；结果/报告 | 角色适配、BED/PAF 一次坐标转换；可选 GenBank/S4 代码注明依赖与 SKIP，不计为运行通过。 |
| 连接准备与统计 | [ribbons](../website/zh/ribbons.md)；`ribbons`；density | filter/deduplicate/merge/bundle 及 density/bundle stat；解释 8/8/8/28 输出与 envelope 非重新比对。 |
| 区间与点连接 | [ribbons](../website/zh/ribbons.md)、[modes](../website/zh/modes.md)；`pointlinks`、`collinear`；pointlinks | 点端点、箭头和显式扩展美学 colourbar；合成 support 不解释为生物学测量。 |
| Feature 与 Position | [annotations](../website/zh/annotations.md)；`featurestack`；featurestack | 重叠 feature、未知链、形状映射与共享 stack Position。 |
| 自动与固定标签 | [labels](../website/zh/labels.md)、[plasmids](../website/zh/plasmids.md)；`labellayouts`、`labels`、`plasmid`；labelauto、labelradial、labelarc | 相同输入展示 auto/radial/arc；固定角度与自动布局归属、隐藏/缩小及最终尺寸限制。 |
| 轴、参照条与区域 | [annotations](../website/zh/annotations.md)、[views](../website/zh/views.md)；`decorations`、`focus`；decorations | region、axis、bp scale bar、break 和 source axis 的不同语义。 |
| 定量轨道与数值域 | [tracks](../website/zh/tracks.md)；`tracks`、`multitrack`；multitrack | bar/line 共用数值域；独立分类带；物理宽度、baseline、整数位置和边界诊断。 |
| 引物、酶切及 feature 搜索 | [recognition](../website/zh/recognition.md)、[plasmids](../website/zh/plasmids.md)；`recognition`、`plasmid`；plasmid | 圆点跨原点 primer、酶切计数/筛选、精确参考匹配与标准遗传码/启发式限制。 |
| CDS 与连接投影 | [projection](../website/zh/projection.md)、[homology](../website/zh/homology.md)；`projection`、`linkprojection`；结果/报告 | 分段 CDS feature 投影与成对 protein link 投影；ungapped 必须显式断言。 |
| 显式 cluster 与 loci | [homology](../website/zh/homology.md)；`homologs`；homologs | 显式 membership 与 pairs；报告成员/比较、paralog 默认省略、loci 方向扩展。 |
| 局部、分段与源映射 | [views](../website/zh/views.md)；`focus`、`segment`、`mappingposition`；focus | focus、segment、局部翻转及普通 ggplot 注释的源位置映射。 |
| 标度、图例与主题 | [styling](../website/zh/styling.md)、[adjustments](../website/zh/adjustments.md)；`styling`、`adjustments`；styling | 映射与固定值、用户 scale 优先、角色图例、主题；别名共用说明。 |
| 预览、结果与排序 | [output](../website/zh/output.md)；`save`、`preview`、`provenance`、`layoutresults`；结果/报告 | 标准 ggsave、同尺寸 preview/diagnose、完整 getter/分层 export、count 权重 opt-in 排序。 |
| 迁移与常见问题 | [migration](../website/zh/migration.md)、[faq](../website/zh/faq.md)；`migration`；结果/报告 | 可运行替代、移除版本/守卫、FAQ 和发布/部署边界。 |

## 已执行的有限验证

- `Rscript tools/documentation-plots.R`：执行共享例子并生成 19 张显式尺寸 PNG，
  附带 8 个非图例 build。随后仅将 density 调整为 10 × 8 in，定向重新导出该图。
  图与代码 hash 见 `website/assets/plots.csv`；观察结果见 `example-results.json`。
  已查看新增综合图、点连接、密度、堆叠、多轨道、同源、装饰、三类标签和样式图。
- `Rscript tools/validate-documentation.R`：原始 FASTA/BLAST/GFF3 为 4/32/266 行；
  一次 8 × 5 in PNG preview/diagnose，panel 为 measured。该状态不等于布局完全无碰撞。
- `Rscript tools/freeze-api.R check`：130 导出保持候选签名、别名和上游对象身份。
- `python3 tools/check-documentation.py --source-only`：25 配对、共享代码、资产 hash、
  18 功能簇、130 导出 / 93 Rd 以及 README/vignette/迁移源码同步通过。
- 可选 Biopython 和 GenomicRanges 本环境未安装，相关示例标为 SKIP，未假称运行通过。

- 生成站点静态检查：每语种 123 个 HTML（包含 reference 跳转入口），122 条全文
  搜索记录；完整内部链接、锚点、逐块代码与原文、图片字节及两语种锚点一致性通过。
- 浏览器：1280 × 720 与 390 × 844，专题 #read 和 API #arguments 同页切换通过；
  汉语“标签”搜索 14 页、position_track 搜索 7 页；键盘跳至 main、原生代码折叠、
  Light/Dark/Auto、窄屏菜单/表格和始终可见语言入口通过。修正汉语原生标点锚点后，
  新会话中参考页/首页控制台无 warn/error。证据见 `website/browser-verification.json`。
- 98 条精确历史路径映射与 80 条无当前一一对应的历史提示，未伪造函数别名。
- 构建阶段明确保存代码正文，修正原生下载链接与数字参数重复 ID；本包代码自动
  链接指向本地当前参考；汉语原生标题按英语对应标题归一锚点，避免标点选择器错误。

最终清洁源码构建将在实施提交后复核；以上浏览器检查使用相同实现的本地预览。

## F52 例子暴露的正确性修复

合法空 primer 命中表原先在引物/标签输入转换中作标量赋值，随后空组件映射又可能
求值不存在的 x/y。现在输入转换保留零行表，空组件在构建局部副本取消映射和继承，
源图层与公开导出登记保持。单次公开 build 验证骨架保留、两个空层零行，新增一个
最小严重缺陷回归案例但未运行 testthat。修复提交 `97034f7`，NEWS 已记录；没有
新增、删除、改名函数或形参，非空默认行为不变。

## 来源与限制

- 图片来源、筛选、尺寸、代码和原始输入 hash 可追踪；输入 ZIP 不改写原始记录。
  原画 PNG 保持 1254 × 1254，3048867 字节；无损 WebP 2127126 字节，相对输入
  3195879 字节缩小约 33.44%，逐像素 RGBA 校验不变。favicon 明确为缩略图。
- 源码检查、网页模板构建与相关例子不代表正式发布验收。没有运行 testthat、
  `R CMD build/check`、完整英语/汉语 vignette 渲染或 PNG/PDF/SVG 发布设备矩阵。
- 未 push、tag、GitHub Release、修改远端 Pages 设置或部署。冻结 `docs/` 不变；
  线上切换和旧站删除仍等待独立明确授权，v0.16.6 未启动。
- pkgdown/Pandoc 在构建中发出上游 `--mathml` 弃用提示，不影响本轮 HTML 生成；
  不在本版修改外部软件。浏览器检查为有限的本地检查，不作全浏览器/无障碍认证。
