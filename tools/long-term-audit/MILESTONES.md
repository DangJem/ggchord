# v10 后续里程碑与依赖

2026-11-10 用户选择以 v10 重排。系列号仅组织工作，不是日期或全部功能交付承诺。
v0.17.0 首轮审计已完成；用户随后已授权开始正式开发，首批统一 Feature 输入与来源桥接已落地。下版创建与版本
分支整合继续按 AGENTS.md，不能未经批准把本分支合并后自动串行开发所有版本。

| 阶段 | 工作包与优先级 | 进入条件 | 退出证据 |
| --- | --- | --- | --- |
| 首轮 v0.17.0 | 当前源码清单、兼容设计、来源/实例桥接原型、网站计划、README | v10 及冲突选择已确认 | 本目录报告、CSV、定向结果；公开运行时不改 |
| Foundation / v0.17.x | P0 来源/实例注册、feature_data 兼容；Track col/identity 与 count；P0 标签根因控制变量；P1 纵向 Collinear、Sankey MVP、轴/equal、pident/Guide | v0.17.0 开发已授权；按已审阅迁移方案和输入契约逐批推进 | 旧调用可解释迁移；公开最小成功路径；真实来源/数值与局部输出；未实现部分明确 |
| View Expansion / v0.18.x | Dotplot、spokes Radial；重复实例完成 Collinear/Chord/Radial；Track 扩展组合和导出 | 基础注册/路由已验证，轨道局部统计与值域稳定 | HSP 来源与端点、DA 缺失对照、跨图层支持矩阵、无静默重路由 |
| Composition / v0.19.x | Fusion F0→F1，panel 注册/变换；尺寸/标签/Guide 联合排版 | 单图坐标映射和测量可信，真实跨图数据可用 | 实际区间锚点，多尺寸无漂移，失败诊断；F2/F3 单独门槛 |
| Semantic Enrichment / v0.20.x 及后续 | 适配器、SV/domain/motif；短序列双层方案、区域放大 M0→M3、Recipes | 有真实字符/对齐样本、统计范围和局部坐标契约 | 不虚构 mismatch；source↔alignment 映射；局部重绘与 leader 语义准确 |
| API Consolidation / 1.0 前 | 集中命名/责任/兼容/弃用与支持范围；文档迁移 | 实际实现范围已确定，不把所有候选强制塞进 1.0 | 正式 API 契约、迁移与完整候选验收计划，未完候选明确延后 |
| v1.0.0 | 稳定 API 正式发布 | 上述已纳入稳定范围的门槛完成；用户明确发布授权 | 按仓库要求完整包/设备/双语/平台验收与发布事实 |

依赖链：source/provenance → instances → routed endpoints → Coord projection → panel registry
→ global composition。Track source Stat、标签/字体测量、Ribbon 路径可在各自隔离边界推进；
联合 Guide/Panel 重排贯穿多个阶段，不单独维护 viewer 算法。逻辑可并行不代表委派代理授权。

优先稳定落地的是入口兼容桥接、col 迁移、首批局部 count 与纵向主轴策略；高风险先做
原型的是共享实例来源键、曲线统计组合、反向 Ribbon 与渐变、全局标签回溯、跨 panel
物理变换和自适应。暂不做无真实用途的新 geom、任意第三方 Geom/Stat 透传、完整组装图
Coord、Logo 字形重造、任何交互图/在线编辑器。

网站独立工作包：先导航/安装兼容事实与开发记录归位，再 Markdown/共享生成源，接着
双语任务课程与已实现 API 示例，最后构建 HTML 与浏览器 320/375/650/768/1024/1280/1440
回归。绘图运行时目标未实现前不写进正式可运行教程；Docs 历史归档不重建。

包级检测门槛仍是用户明确 GitHub 正式发布授权。日常各阶段只做直接相关的局部验证；
长期完整质量矩阵登记为验收计划，不将本地设计交付冒充发布验收。

当前实现、验证和未完成项见 [v0.17.0 开发记录](DEVELOPMENT-v0.17.0.md)。
