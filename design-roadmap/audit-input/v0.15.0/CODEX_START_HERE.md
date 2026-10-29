# Codex：从这里开始

请先通读 `GGCHORD_API_ARCHITECTURE_AUDIT_AND_REFACTOR_PLAN.md`，这是一份**审计规范和阶段性重构方案**，不是删除/重命名全部 API 的自动授权。

## 明确的最新用户决定

1. v0.15.0 已合并进 `main`；以实时最新 `main` 和 `AGENTS.md` 为准。
2. ggchord 已具备 chord/circular/collinear 三种主要大类，暂不加入第四类。
3. **从 v0.16.0 开始进行真正的 API 重构**，分版本实施；可适度增加必需的新 API。
4. v1.0.0 是 API、源码、迁移、文档与网站全部稳定的**正式发布**，不是首次启动重构的版本。
5. 本轮先做**源码级全量审计 + 目标架构提案 + v0.16.0 范围计划**，不要未经批准批量删除接口。

## 推荐工作顺序

- 先检查 `AGENTS.md`、`DESIGN-ROADMAP.md`、v0.99.0/v1.0.0 既有规范，与最新用户决策的冲突显式记录。
- 使用 `Rscript audit_ggchord_api.R . design-roadmap/audit-generated`（必须在可用的 R 环境执行）；生成导出、内部函数、全部形参、S3 清单。脚本输出是候选**静态**台账，需补充别名、ggproto、`quote`、非标准求值、真实默认值/单位及调用路径。
- 阅读 `api_export_baseline.csv`，逐行按 KEEP/FIX/RENAME/MERGE/SPLIT/INTERNALIZE/REMOVE/ADD/DEFER 给出**有证据**建议。
- 逐参数生成 type/default/units/NULL/NA/semantic/validation/dependency/coord-support/migration 列；一项也不能跳过。
- 对每一项破坏性变更先提供：旧调用、新调用、迁移例子、受影响用例、错误提示、删除版本，以及用户需要作出的明确决策。
- 首批优先修复可复现正确性问题；避免一开始重写整个图形引擎。
- 不执行自动推送或发版，也不未经明确授权运行完整发布验收。

## 应先向用户提交什么

1. v0.16.0 版本计划（精确列出范围与非目标）。
2. 全量 API/参数清单与结构性问题总表。
3. 第一批低风险实施方案与回归样例。
4. 待用户决策的破坏性 API 命名、删除和合并问题清单。
