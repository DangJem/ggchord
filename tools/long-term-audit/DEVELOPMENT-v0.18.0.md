# v0.18.0 本地开发记录

状态：v0.18.0本地开发完成。运行时、定向验证、双语文档、配图、本地网站与新增页面浏览器检查通过；尚未整合main或正式发布。
包版本0.18.0.9000。最近正式发布仍为v0.16.6，本记录不构成正式发布验收。

用户于2026-11-11批准将v0.17.0 f183d67快进整合本地main，并从main建立v0.18.0；
随后要求持续开发至本版完成。全程使用OneDrive原checkout，未创建worktree。
本版工作未获main整合、push、部署、tag或正式发布授权，均不执行。

## 实现结果

- `coord_dotplot()`和`geom_alignment_segment()`直接复用标准序列表与HSP表，
  无需`geom_seq()`或新专用数据结构。保留四端点、逐层数据、稳定HSP ID、元数据及源行。
  两轴独立选择/排序/翻转，采用一基原生bp/aa单位、拼接间隙和区段分隔线。
- `coord_radial_synteny()`复用Collinear源空间Feature、Track、统计及实例路由，
  最后投影到等角辐条；默认统一真实长度比例，equal必须显式选择。正offset为逆时针法向，
  物理Track宽度按最终辐条方向计算，源统计在实例复制之前完成。
- Radial默认环形邻接，仅展示已有首尾HSP，无证据的通道保留为空。N=2每原行一次，
  N=1拒绝。命名重复实例必须明确between；已声明source_hsp_id的重复展示需要显式允许。
- Radial完整前缘采用一基闭区间边界`[min-1,max]`并保留输入方向，单位置区间有非零宽度；
  轴与位置反查仍从1开始。反向前缘、扭转及完整区间不是中点替代。此约定只用于新增视图，
  不重写已有视图区间语义。默认连接留在环形通道内，显式控制点使用最终x/y。
- Ribbon继续默认Identity填充；端点渐变、装饰subject尖头复用共享绘制，仅显式开启。
  数值pident保留0/50/80/90/95/100，用户scale/theme优先。序列名移到辐条外端，
  显式水平标签保持水平；数字Track轴保留原有水平渲染。
- Dotplot导出alignment/axis_layout/connection_report；Radial导出局部切向/法向、角度、
  内孔及完整原始/投影前缘。反查明确选择Dotplot轴；Radial须提供骨架上的x/y，
  支持聚焦源原点及命名实例，不从一个x值猜径向位置。
- 两新视图共享已知/未知方向检查：未知方向警告省略并报告，方向与规范化端点冲突报错。
  修复Track读取tibble时对不存在的可选.track_width字段产生的无关警告。

## 支持边界

| 组合 | Dotplot | Radial |
| --- | --- | --- |
| HSP alignment segment | identity线段 | 明确拒绝 |
| 骨架、gene/feature、region、break | 明确拒绝 | 源局部几何投影 |
| 手动/自动gene与feature标签 | 明确拒绝 | 共享实测队列；按实际尺寸检查隐藏 |
| geom_track与已有Track快捷层 | 明确拒绝 | 共享源Stat、数据Position和物理band |
| Ribbon、点连接 | 明确拒绝 | 默认环形邻接或显式关系/路由 |
| 原生ggplot2注释 | 拼接原生单位 | 最终x/y布局单位 |
| facet、共享branch | 明确拒绝 | 明确拒绝 |
| 圆图专有引物、限制位点、中心标签 | 明确拒绝 | 明确拒绝 |
| 共享直线绝对比例尺 | 明确拒绝 | 明确拒绝；用每条辐条的原生轴 |

任意第三方Geom/Stat、复杂Fusion、逐碱基错配轨迹和多窗口拼为一条辐条不属本版。
Dotplot没有把单序列Track强套到双轴上；Radial负责本版共享Track/Feature扩展。
密集标签可能按图层策略隐藏，来源与可见性报告保留；不声称任意画布都无碰撞。
显式控制点和直线点连接可能穿越中心或其他辐条，不构成默认通道的保证。
旧构造验证器仍可能提示合法反向HSP的降序端点；教程明确不要据此排序并抹去倒位。

## 直接验证证据

使用以下局部脚本，不运行testthat或包级发布流程：

- `Rscript tools/long-term-audit/check-dotplot.R`：真实EMALE PAF 11条、核酸BLAST 7条、
  蛋白BLAST 1条，四端点往返、原始来源不变和实际gtable通过；合成多轴/翻转/多层、
  用户scale/原生theme、选择省略、未知/冲突方向、聚焦、重复授权、空层及旧Chord build通过。
- `Rscript tools/long-term-audit/check-radial.R`：真实两辐条PAF 11条保留一次；
  四序列核酸BLAST按环形邻接保留16/32条；四蛋白BLAST保留8/16条并保持aa。
  完整闭区间宽度及反查通过；合成有/无D-A、N=1/2/3/4/5、单位置、反向、equal、
  A-B-A、源统计成员、聚焦、空关系毫米Track、重复身份授权及一个Ribbon密度Stat通过。
- `Rscript tools/check-layer-compatibility.R dotplot radial`：新增两视图各24个代表图层，
  共19次构建成功与29次预期拒绝，逐行原因在layer-compatibility.csv；此前69行作为
  已完成版本证据保留，不冒充本轮重新执行的默认参数穷举。
- `profile-view-expansion.R`：明确合成复制PAF用于限定规模检查，Dotplot 10,000条及
  Radial 1,000条全部保留。实测机器耗时和built对象大小见view-expansion-profile.json，
  对象大小不是RSS，不设跨机器固定秒数门槛，不代表新增生物学观测。
- `check-view-source-bridge.R`：同一真实PAF的11条HSP在Chord、Collinear、Sankey、
  Dotplot及Radial中逐一实构gtable，源表身份、四原始端点及全部源行均保持不变。
- 单个2.5亿位置长序列的径向端点往返通过；这是数字量级检查，不是压力设备矩阵。
- 已在现有plot API测试文件增加两种新视图各一条最小公开行为，留待获授权正式发布时执行；
  本轮没有执行测试套件。

## 文档、图片与本地网站

- 同步NEWS、两语迁移说明、roxygen/Rd、两语view-expansion课程、支持矩阵、导航、
  API教学映射、任务索引和参数图谱。136个导出、1,545个formals、100个Rd主题；
  114个本包导出函数usage与真实formals匹配，952段汉语翻译保留代码/标识/链接。
- 49个教学页在独立环境准备代码通过；两语代码一致，61对网站页面源检查通过。
  新输入来自原始PAF或随包真实序列/HSP/源链CDS，筛选和计数含义均在教程中说明。
- 从最终运行时重现176张PNG。原有174张与提交基线字节相同，只新增viewdotplot.png
  与viewradial.png；没有重复提交中间PNG。补查并更新离线教程3张已落后于原PNG的JPG
  缩略图及来源指纹。README仅更新当前开发版本，未来方向保持原文。
- 本地网站新增课程及三个接口主题的双语页面：课程在320/768/1440宽度，三个帮助主题
  在320/1440宽度，共18组合无整页横向溢出。高清Radial 1440×1296图加载、原始大小、
  关闭、站内新函数/参数搜索与锚点跳转通过；已恢复临时viewport设置。
- 本地网站仅预览，冻结docs未改。最终静态HTML链接、锚点、搜索与资产复核通过。

## 版本与权限事实

Dotplot首批提交5fd82fd，Radial运行时提交8c3e047。文档/最终配图仅在本版分支集中提交。
未修改已发布历史，未清理Git对象或执行激进gc。当前main仍为f183d67。
未运行testthat、R CMD build/check、vignette编译或PNG/PDF/SVG发布设备矩阵；
未整合本版到main，未push、部署、创建tag/Release或CRAN发布。
