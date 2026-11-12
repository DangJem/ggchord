# v10 长期任务书实施状态

核对日期：2026-11-12。对应任务书：`ggchord_长期开发路线图与统一架构_Codex任务书_v10_独立执行版.md`。
基线审计为33236ef；Foundation、View Expansion、Composition已分别完成本地开发并按授权
整合main；v0.21.0已按授权整合并推送，main/origin/main为69556d8。当前v0.22.0约定范围本地开发完成，尚未整合或推送。
最近正式发布仍为v0.16.6。本文中的“完成”指明确边界内的本地实现与定向验证，
不表示包级发布检测、所有候选功能或任意输入组合完成。

## 按能力核对

| v10条款 | 当前事实 | 证据与剩余工作 |
| --- | --- | --- |
| §1.3、§10.1–8 首轮审计与统一架构方案 | 已完成首轮审计、分类矩阵、迁移设计和隔离原型 | README.md、DESIGN.md、api.csv、parameters.csv、capabilities.csv、prototype-results.json；均为当时基线，不冒充当前API清单 |
| §2、§10.2 来源/实例/路由 | 已实现共享来源、命名实例、Coord前路由与来源导出 | DEVELOPMENT-v0.17.0.md、check-view-source-bridge.R；跨进程任意同值表身份不作推断 |
| §2.5、§10.9 Feature统一 | 已实现feature_data、旧gene_data兼容及互斥、显式分类 | feature-input-results-v0.17.0.json；稳定化阶段仍须决定兼容长期策略 |
| §3.6–3.7、§10.10 Track | 已实现14种首批几何、identity/count/bin/summary、group/after_stat/数据Position、物理band和来源统计成员 | track-grammar-results-v0.17.0.json；不是任意第三方Geom/Stat透传，facets、smooth、boxplot/violin等未支持 |
| §4.2、A.1 纵向Collinear | 已实现源空间到纵向投影、方向/反查与图层组合 | foundation-views-results-v0.17.0.json；本版补缺失角度崩溃修复 |
| §4.3、A.2 Sankey | 已实现query左/subject右镜像、真实区间、默认自身过滤 | 同上；不实现多阶段流量模型，不支持共享branch/重复between路由 |
| §4.4、§10.4 Dotplot | 已实现HSP线段、独立轴、源端点/单位/来源与反查 | check-dotplot.R、DEVELOPMENT-v0.18.0.md；精细错配轨迹及单序列Track不在当前支持范围 |
| §4.5、§10.5/11 Radial | 已实现spokes、首尾邻接、有/无D–A对照及共享Feature/Track/实例 | check-radial.R；N=2每源行一次，N=1拒绝，无虚构首尾连接 |
| §4.7–4.8、§10.12 等长/轴 | 已实现可选equal、真实长度保留、1/L端点优先及圆接缝策略 | foundation-views-results-v0.17.0.json、plasmid-axis-results-v0.17.0.json；显式刻度优先，不提供误导性共享绝对比例尺 |
| §4.11、§10.16 A–B–A | 已实现seq_order/命名副本、between与dplyr筛选，在Collinear/Chord/Radial共用路由 | v0.17/18开发记录；Circular/Dotplot不适配重复实例 |
| §5、A.3 Fusion | F0/F1已实现真实区间、中央Circular/Chord、四方向横纵Collinear、毫米导出与目标尺寸重建 | check-fusion.R、v0.19.0-results.json；F2/F3、外围Sankey/Radial/Dotplot、八方向/自由角度、全局避障未完成 |
| §6.1、A.5 Ribbon | 已实现显式端点渐变、装饰subject尖头、Identity六刻度与用户scale优先 | ribbon-flow-results-v0.17.0.json、identity-guide-results-v0.17.0.json；按后续用户决定默认仍非端点渐变，不能按原提案声称默认渐变 |
| §6.2、A.6 标签 | 已修复30标签/混合层切换、共享实测队列、圆轴字形；真实13质粒/旗舰与代表尺寸有证据 | foundation-label/size-results-v0.17.0.json；不是全局零碰撞保证，未完成与SnapGene同输入导出的完整公平对标，全部引线穿字/内部框指标仍待完善 |
| §6.3.1、§10.15 Guide | 已实现实测预算、Identity保序不等宽分段、最终方向和显式设置优先，已有离散布局原型 | guide-results.json、foundation-guide-results-v0.17.0.json；通用多Guide联合分行/调向求解仍未完成 |
| §6.3、A.7 Adaptive | 部分：共享render context、panel实测、标签/guide按目标尺寸重建；Fusion adaptive分配中央区域 | v0.19开发记录；通用fixed/scale/adaptive样式、有限字号/线宽缩放、全部视图多尺寸质量报告未完成 |
| §6.4 灵活参数 | 已有标量/实体/实例及按链继承，用户映射/scale优先；已同步教程 | v0.17/18记录、website双语flexible-parameters课程；未承诺任意第三方层全兼容 |
| §6.5 Recipes | 已有真实科研综合示例与Fusion配方，继续复用基础API | tools/documentation-examples.R；不是独立模板产品；v0.22.0已提供显式源区间magnify入口 |
| §6.6 适配器 | 本版新增BEDPE保留原字段/方向/缺失策略；FASTA可显式保留真实字符 | check-bedpe.R、check-semantic-detail.R；VCF/SV、BigWig/bedGraph、domain/motif等格式专用导入器未完成，见ADAPTER-AUDIT-v0.20.0.md |
| §4.9、§10.13 短序列 | v0.21.0已实现真实gapped FASTA、源/列双映射、显式源验证和Cartesian局部geom_residue | DEVELOPMENT-v0.21.0.md、v0.21.0-alignment-results.json；不含MSA算法、CIGAR全格式解析、CDS翻译/phase及宏观字符投影 |
| §4.10 语义放大 | v0.22.0已公开magnify_ggchord：真实源区间、独立Collinear细节、毫米高亮与装饰leader、精确导出 | DEVELOPMENT-v0.22.0.md、check-magnify.R；要求显式detail及统计口径。M1–M3完整目标、跨原点选区、任意inset、多级嵌套未完成 |
| §7、A.4、§10.8/14/17 网站/README | 导航/Markdown源/双语教学重构已实施，README未来方向保留，新API同步实际教程 | v0.17/18/19开发记录与浏览器JSON；62对页面为当前源规模。代表页响应式通过，不是全站逐页逐宽度视觉认证；本版不重建旧图或冻结docs |
| §8.3、§9、1.0稳定化 | 待实施 | 未集中冻结全部新API、兼容/弃用策略及完整候选发布验收；版本分支开发完成不等于v10完成 |
| X 非目标 | 持续遵守 | 不开发交互HTML绘图、Shiny/在线编辑器，不运行分析引擎伪造比对，不自动push/部署/Release |

## §10 首轮18项交付的解释

1–8、9–12、14、16–18已有相应审计、设计、原型或后来实现证据，详见上表。
第13项已有真实对齐映射与公开局部字符工作流，不能视为完整MSA生态或宏观联动；第15项有Guide规格/原型与部分运行时，
不能视为通用自适应Guide求解器。第18项由本表持续区分完成、部分、研究和未开发。
首轮交付完成与最终产品全部能力完成是不同状态；不按函数数、版本数或行数计算总体百分比。

## 后续依赖

真实逐列对齐映射与局部字符层已在v0.21.0落地；v0.22.0已提供显式源区间语义细节组合；后续评审宏观字符投影、CDS细节与自动语义重建；新格式逐项确定
来源与缺失语义。Fusion F2/F3、全局标签质量及通用Adaptive各有独立门槛。
实际纳入范围确定后进入1.0前集中API治理。后续具体版本号仍须按路线图与整合授权确定，
本表不把所有候选强制塞入v0.20.0，也不把延后写成已完成。

v0.21.0仅完成相关新图/示例/源文档的定向验证；整站旧lesson engine与178张旧图依赖指纹过期，完整整站检查未通过。后续须重现，不将旧图冒充当前输出。

本轮main推送成功但远端CI失败；限制位点直线引线回归已在v0.22.0本地修复，尚未推送验证。旧文档资产待办仍保留，见CI-69556d8.json。
