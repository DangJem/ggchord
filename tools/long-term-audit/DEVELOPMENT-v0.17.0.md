# v0.17.0 开发记录

## 统一 Feature 输入与来源桥接

2026-11-10 用户要求恢复 README 原详细未来方向，并开始 v0.17.0 正式开发。
本批从 `v0.17.0` 已有审计提交继续，不新建版本分支或 worktree。包版本已设为
`0.17.0.9000`；首轮基线 CSV/JSON 保持冻结，不作为新实现的最新能力表。

已实现：

- `ggchord(feature_data=...)` 和 `validate_ggchord_data(feature_data=...)`；两个非 NULL
  注释入口互斥。gene_data 第三位置及旧无类型、已知链语义保留，只在构造时 message
  一次，后续 build 不重复，未安排删除版本。
- 缺少 strand 列按未知 `.` 处理，通用注释接受 `+/-`，可输入 GRanges 和 aa 特征；
  原始表不被改写，额外列与 list-column 保留。
- 默认基因/基因标签从 feature_class/type/role（依此次序）选择 gene/CDS；tRNA/rRNA
  留在通用特征。缺少类别信息时不猜测，提示显式提供 gene 图层数据。
- geom_feature 与 feature labels 默认使用全表；独立图层数据/数据函数优先。显式
  互斥过滤可以组合 gene 与其他 feature，不根据另一层的存在隐式删行。
- 默认输入贯通几何、标签、避障、导出。目录注释颜色继承与显式输入一致，用户标度
  仍优先。source_row 元数据覆盖计算来源行号的缺陷已复现并修复。
- `source_registry` 注册构造来源；`layer_sources` 保留 resolved input 行到源表/源行的
  旁表；`annotation_duplicates` 报告同一构造记录出现于多个几何层，标签不计重复。
  来源 ID 只在图对象内有效；显式表/函数为 layer scope，不按相同内容猜源身份。

主要源码：`R/annotation-input.R`、`R/ggchord.R`、`R/validate.R`、`R/zzz.R`、
`R/ggchord-build.R`、`R/link-obstacles.R`、`R/layout-export.R`、`R/diagnostics.R`。
相关 roxygen/Rd、NEWS、双语迁移同步。复用/扩充现有 feature 测试案例，只增加一个
统一输入的代表性公开测试；本轮没有执行 testthat。

## 可运行的混合注释语法

下列是明确标为合成的最小 API 输入，不冒充实际科研观测：

```r
library(ggchord)
library(dplyr)
seqs <- data.frame(accver = "A", length = 1000)
annotations <- data.frame(
  accver = "A",
  start = c(100, 400), end = c(300, 450),
  strand = c("+", "."), type = c("CDS", "promoter")
)
p <- ggchord(seqs, feature_data = annotations) +
  geom_seq() +
  geom_gene() +
  geom_feature(
    data = annotations |> filter(type != "CDS"),
    shape = "block"
  ) +
  coord_collinear()
```

统一输入不要求转换为复杂对象。类型缺失/NA 作为未分类通用记录保留；不能据已知链
就把它认作基因。几何 source_row 仍索引解析后输入；用户同名元数据保留在 original_data。

## 定向验证

执行 `Rscript tools/validate-feature-input.R`，结果保存在
[feature-input-results-v0.17.0.json](feature-input-results-v0.17.0.json)。覆盖混合表互斥
图层、全表/标签继承、默认基因筛选和重复显示诊断、单独基因标签、旧位置参数与消息、
错误输入/结构化校验、无类型/无链、aa、相同源记录不去重、空基因子集、用户 source_row
元数据隔离、真实 pUC19 配色等价和 GRanges（本环境已安装）。这些是直接相关的局部
验证，不代表全部参数或发布矩阵。

真实 pUC19 识别出的9条注释通过统一输入构建圆图、导出来源，并在一个8×6英寸 PNG
设备下绘制和人工查看；继承/显式输入的实际 feature_fill 相同。此次查看只核对新输入
路径，既有密集标签排版问题没有在本批解决，不宣称视觉质量专项通过。

README 同源一致性与 Git diff whitespace 检查通过。未执行 testthat、R CMD build/check、
vignette 重建、网站构建或 PNG/PDF/SVG 发布设备矩阵。

## 当前边界与后续依赖

整版仍在开发。本批不声称 A–B–A/实例路由、通用 Track/count/col、纵向 Collinear、
Sankey、equal/端点轴、Ribbon 渐变/箭头、通用自适应或标签算法已经完成。
下一批沿 Foundation 依赖继续来源/显示实例及 Track 统计空间；工作流 helper 的统一
输入贯通和网站示例迁移也需后续推进，现有 gene_data 调用继续兼容。

首轮审计报告与设计保留其基线事实；README 详细版的恢复优先于任务书中简写建议，
未来方向仍保留“规划/探索、未必已可用”的说明。网站存量版本标识、教学输出与配图
指纹需在后续文档工作包同步，不把旧生成物称为新版本已验证输出。

发布状态：仅本地版本分支开发；未整合 main、推送、部署、打 tag/Release 或提交 CRAN。

## Foundation 运行时推进（2026-11-10，持续开发中）

已落实通用 `geom_track()` 的14种几何、源空间 identity/count/bin/summary 与原生数据
Position；Stat 在实例复制前执行。数值 bar 已迁移 col，旧物理 Position 写法兼容。
带出原始成员与统计结果。Feature 的 prepare/clean/focus/segment 统一入口已贯通。

已落实命名/重复显示实例、实体继承与实例覆盖、between 名称/序号、歧义拒绝、声明 HSP
身份的重复显示检查、路由/省略报告；不按相同值猜观测身份。同一 A–B–A 对象可切换
Chord/横纵 Collinear。新增镜像 Sankey，自身过滤与互反输入区分，坐标反查已贯通；
首批不支持共享分支，显式请求会报错。equal 为主动显示投影，导出源长与投影系数，
禁止误导性共享绝对比例尺。默认1/末位、中间物理空间选择和闭圆接缝合并已接入。

Ribbon 默认沿曲线参数由 query 色变为 subject 色，subject 端装饰收尖；固定 fill、映射
与用户 scale 优先，tip=none 保留全端面。显式 pident 使用固定域和五段数值比例 guide。
自动离散 Guide 测量文本选择列数，并沿用联合 panel/gtable 测量；显式主题/guide 优先。

直接验证：`validate-track-grammar.R`、`validate-foundation-views.R` 与首批 Feature
脚本通过。不运行 testthat。新增最小公开行为案例保留待正式发布授权后的包级执行。

标签专项仍在本轮持续修正：已移除30标签/混合primer的隐式切换，短特征局部失败交给
外围队列，长同位酶名按真实字体换行，按实际空间平衡接缝附近负载。首轮发现小画布
物理迭代未收敛，旧“估计框零碰撞”不足以证明最终文字不重叠；正在用最终设备度量
复核。139条真实酶切来源合并为同位标签，应报告 shared-label 成员，不误作隐藏。
双语教学、参考翻译和网站仍在同步，当前不宣称整版已完成或已发布。
