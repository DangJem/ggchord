# 三模块教程重构：本地审计记录

首批基线：`647cd3b4655b299346237c84b804a531bb7ec45a`；当前矩阵基线：`1dc62fee9320ae57b022364b14485584cd9471ce`。版本 `0.16.7.9000`，2026-11-05。
矩阵读取当前工作区，包含本轮未提交文档；SHA 标记其父基线，不表示这些文档已经提交。
正式版为 0.16.6，核查来源为 https://stat.ethz.ch/CRAN/web/packages/ggchord/refman/ggchord.html。
不修改 R/、默认值或运行依赖；未整合、推送、部署和发布。

## 如何使用矩阵

- `teaching-api-v0.16.7.csv`：当前命名空间所有 130 个导出，区分 owner、源码显式别名、相同对象候选、Rd、页面与共享代码候选、教学等级、模式和人工判断依据。
- `teaching-parameters-v0.16.7.csv`：1474 个公开 formal 条目（含重导出、同名参数及各函数的 `...`），默认值来自实际 `formals()`，不是手抄旧签名。
- `teaching-summary-v0.16.7.json`：版本、基线与等级统计。
- `../../website/teaching-reviews.json`：人工确认达到 L2/L3 的少量主讲 API 与参数证据。
- `../../website/teaching-tasks.json`：任务、关键参数、模式、主讲位置、学习层级；生成双语索引。

这里的 L0/L1 是证据下限，不能解读为已经逐项人工判定“完全没教”。L2/L3 仅用于有明确人工审阅记录的条目；不把函数出现、共享代码被调用、图像存在或一次执行自动提升为充分教学。没有宣称任何 L4 或全部核心 API 已达 L3。模式尚未逐项确认的行明确标为未审阅，不从函数名猜支持情况。参数表保留完整默认表达式，允许值与类型须回查其 Rd，尚未自动抽取完整参数语义。

## 优先级与处理结果

| 优先级 | 问题 | 本轮结果 / 未完成部分 |
| --- | --- | --- |
| P0 | 安装页声称使用 0.16.6，但示例含开发版行为 | 修正双语安装页及同步离线源；明确本地源码与 CRAN 对应文档，不假定未推送的远端已可安装 |
| P0 | 首课共享对象可能造成跨页依赖 | 三个全新 Rscript 进程独立执行并 build 成功；38 个页面独立环境执行成功 |
| P0 | 兼容性容易按函数名误猜 | 新增基于实际 Rd/build 校验契约的兼容表，明确单圆图拒绝连接、直线标尺与断点仅限共线性 |
| P1 | 对齐、锚点、自动方向没有完整主讲课 | 新增双语 collinear-layout，六参数表、四种方向、冲突、练习、自备数据与五幅同尺寸实际图 |
| P1 | tile 未在共享课程直接绘制 | 使用原 pUC19 GC 窗口新增真正 tile 示例；与 bar/line 区分输入、数值域和编码方式 |
| P1 | 缺少任务到参数的检索入口 | 结构化 32 项任务生成双语索引并接入导航、搜索和生成一致性检查 |
| P1 | 旧审计依赖 0.16.4 清单 | 新工具读取当前 exports/formals，并展开当前课程代码；保留旧历史工具和记录 |
| P1 | 分支、箭头样式、固定标签、repel、引物、酶切主要参数缺少逐项对照 | ribbon 边界参数表已补；系统视觉对照仍未全部完成，不能记为 L3 |
| P1 | 所有公开入口的语义审计与类型、允许值、适用模式尚不完整 | 已列全入口和实际默认值，10 项人工判断；剩余逐项审阅仍开放 |
| P2 | 全部 scale/guide/identity 家族辅助入口缺少独立解释与例子 | 沿用 styling 主讲，任务索引提供入口；未穷举补图 |
| P2 | 全部完整案例达到研究问题—来源许可—变换—限制规范 | 保留真实来源案例和 13 图库，尚未逐例重写或宣称 L4 |
| P2 | 桌面/移动导航和稳定版对照运行 | 链接/锚点/搜索由工具检查；浏览器人工和安装后的 0.16.6 跨版本运行另记未验证 |

## 课程架构与旧页映射

不批量制造约 50 个页面。三种首课保留原 URL；共线性在行距与跨行课之间新增布局课；其余旧页不移动 URL。

| 知识单元 | 层级 | 现有或新增主讲位置 |
| --- | --- | --- |
| A01–A04 | 基础 | chord-start / chord-data / chord-build |
| A05–A07 | 进阶 | annotations / chord-tuning / chord-geometry；箭头对照仍开放 |
| A08–A11 | 进阶 | chord-routing / chord-labels / labels |
| A12–A14 | 进阶 | adjustments / styling / ribbons；分支与统计深度待补 |
| A15–A16 | 实战 | chord-import / chord / real-data |
| B01–B05 | 基础→进阶 | collinear-start / collinear-data / collinear-build / annotations |
| B06–B07 | 进阶 | 新增 collinear-layout；独立输入、五幅对照 |
| B08–B10 | 进阶 | collinear-routing / collinear-tuning / labels |
| B11–B12 | 进阶 | views / homology |
| B13–B15 | 专题→实战 | tracks / formats / collinear / real-data |
| C01–C05 | 基础→进阶 | plasmid-start / plasmid-data / plasmid-tuning；中心文字细节待补 |
| C06–C10 | 进阶 | plasmid-build / annotations / labels；系统样式对照待补 |
| C11–C14 | 专题 | recognition / plasmids；引物与酶切显示参数对照待补 |
| C15–C16 | 专题 | tracks（新增 tile）；轨道声明与编码边界 |
| C17–C19 | 实战 | formats / plasmids / plasmid-catalog；图库为可选参考 |
| D01–D06 | 基础→进阶 | modes / adjustments / styling；共用概念一处主讲 |
| D07–D09 | 进阶 | labels / annotations / tracks / views |
| D10–D12 | 专题 | output / faq / compatibility |
| 全部任务入口 | 检索 | 新增 task-index，从单一任务表生成 |

## 新增对照目录

| 图片 ID | 输入 | 目标改动 | 设备 |
| --- | --- | --- | --- |
| alignleft | 随包四条完整真实序列与 HSP | 默认 left/manual，无 anchor | PNG，9 × 6 in，144 dpi |
| aligncenter | 同上 | 仅 align=center | 同上 |
| alignright | 同上 | 仅 align=right | 同上 |
| alignanchor | 同上 | 从 left 基线仅增加按 ID 命名的中点 anchor；不是同源声明 | 同上 |
| alignauto | 同上 | 从 left 基线仅启用 auto；0.1 显式重述默认阈值 | 同上 |
| tracktile | pUC19 同一组 100 bp GC 窗口，末窗用实际碱基数 | bar/line/tile 编码比较，非单参数对照 | PNG，9 × 8 in，144 dpi |

实际路径、对象、生成器、输入说明及 MD5 由 `website/assets/plots.csv` 记录；图注双语同源。没有下载额外数据或修改原始来源。

## 验证命令与状态

已通过：

- `python3 tools/check-first-lessons.py`：三类首课分别在全新 R 进程运行并构建。
- `Rscript tools/check-documentation-lessons.R --write-outputs`：38 页独立环境，实际表格与控制台输出刷新。两条已有 Seqinfo 警告来自可选区间转换，非新课错误。
- `python3 tools/audit-teaching.py`：当前导出与函数签名盘点；不运行包级测试。

图像生成、站点构建、全站检查与视觉核对的最终结果见本文件后续记录。包级 testthat、R CMD build/check、vignette 重建和发布设备矩阵未执行，遵守仓库对普通文档任务的边界。此次仅同步离线 Rmd 源码，不 knit。

### 发现的契约差异

`R/coord-collinear.R` 称中心标签、引物与酶切为 circular-only；但定向运行单条 pUC19 + `coord_chord()` + 这三个图层分别 build 成功。兼容页保留这项差异，标注条件与未做视觉验收，不修改核心 API，也不把它自动解释为官方保证。自动引物标签没有在该定向脚本单独运行。后续应由维护者确认公开承诺后统一实现与帮助。

可选依赖限制：本机 `/usr/bin/python3` 未安装 Biopython，GenBank 可选片段跳过，实际控制台快照为 NULL；不得宣称该导入成功。其余已有 Seqinfo 警告如上。这个差异反映当前环境，不是新版修改了输入文件。

新增图像对象独立于既有图片；配图依赖工具显式排除这些独立 chunk 对旧图的影响，
保留所有旧输入、R 引擎和数据依赖。默认新课使用基础外观，已登记到 basic_ids，
不要求它继承完整案例的背景主题。首次全量 PNG 导出在进入新图时因遗漏此登记
停止；已修正并改为定向导出六幅新图，旧图文件由 MD5 与依赖指纹校验。

130 个导出中，108 个为包自身函数，22 个为上游重导出或对象。same_object_as 只说明 R 对象 identical，不代表语义别名（例如 after_stat/after_scale 受 staged mapping 语义区分）；alias 只记录源码中的显式名称赋值。

轨道页虽新增参数表与实际编码图，尚缺单独 width 等几何参数的渲染对照，因此 geom_track_tile 与 position_track 仍保守记为 L2，不计为 L3。

### 最终本地验证

- 定向导出 alignleft/aligncenter/alignright/alignanchor/alignauto/tracktile 成功，6 幅新增图与 98 幅资产登记一致；源码检查通过图片依赖、MD5 与坐标模式核对。
- `python3 tools/build-website.py` 构建 48 对页面、95 个参考主题，无部署。
- `python3 tools/check-documentation.py` 通过双语代码等价、875 个参考译文块、图片、旧链接、锚点与搜索完整性。
- 浏览器检查共线性新课、前后章节与任务索引；390 × 844 手机视口下标题、正文与页内目录可用，目录展开后可见六个锚点；视口已恢复。
- 查看 aligncenter、alignright、alignauto 与 tracktile PNG，布局及图例可读；没有对全部 98 图重新做人工视觉验收。
- 最终复核修正任务索引中 Position 的参数归属与原生蛋白比较模式；索引列的是任务参数族，不应误认为全部属于同一个函数。

完整任务仍未全部完成：剩余主要参数的系统对照、全部 API 的语义逐项审阅、辅助入口 L2 覆盖与完整案例 L4 审核均保留在上方缺口清单。本记录区分本轮已交付改造和原任务的全面验收。

`Rscript tools/check-teaching-defaults.R` 定向比对全部 108 个包自身导出函数的 Rd usage 与当前 formals，参数顺序、必需参数和默认表达式全部匹配；此检查不运行 testthat 或包级发布检测。

## 第二批：几何对照与设计参考（基线 1dc62fe）

用户要求继续审计重构，并在 README 介绍参考的软件、避免提及商业工具。本轮沿用 v0.16.7，无运行时/API/依赖变更。

- 双语 README 增加设计参考：ggplot2、BLAST、Circos、circlize、clinker、gggenes、genoPlotR、DNA Features Viewer、gggenomes、ggview；每项说明借鉴方面并链接原项目。由 website/readme-inspirations-{en,zh}.md 同步生成，检查不含商业产品名称。没有删除必要的数据来源、许可或历史记录。
- 新增 feature-geometry：真实 OQ646790.1 integrase 43,765–46,140 bp，保留原始正链，局部展示 43,000–47,000 bp 并显示源坐标。三张 9 × 4 in 图固定共同头宽 1.6，仅切换 shouldered/flush/triangle。
- 新增 chord-branches：随包 HSP 没有完全相同的 query 端点；明确使用三条 1000 bp 简化序列、两条共享 query 区间的示意记录。对照 none/query、0.15/0.35，随后显式计算 qpos/spos 比较 curve/straight，不篡改实际比对。
- tracks 增加同一真实 GC 数据的 10/4 mm tile 带宽对照，固定 9 × 8 in、窗口、数值域和颜色。position_track 的 width 与数值高度区别可直接观察。
- 通用箭头、标签、灵活参数、映射和样式集中到“通用绘图技巧”。既有页面 URL、章节锚点和三类独立首课保留；任务索引增至 33 项。
- website/parameter-comparisons.json 记录新增 9 图的输入、独立准备片段、改动参数、当前默认、比较基线、尺寸、主讲页与运行状态。配图指纹只依赖该组实际代码、包数据与引擎；定向渲染先验证未请求图片的字节和依赖，不能跳过旧图校验。
- 修复审计 CSV 的行尾为 LF；参数静态证据同时检查展开代码，避免仅因参数写在共享片段而误记为完全未出现。人工教学等级仍单独维护，不用代码出现自动升级。

本轮运行：九张图由 view_ggchord 实际生成；40 个教学页独立环境执行成功；源码检查通过 50 对双语页面、107 张图、130 个导出与 95 个 Rd 主题。基因肩部、平齐头、共享主干、tile 带宽已查看图像。站点最终构建及链接检查结果另见下方。

已审阅的 geom_gene、geom_link_ribbon、geom_link_line、geom_track_tile、position_track 与 coord_collinear 为 L3；等级只覆盖记录中的主讲关键行为，不宣称全部参数或全部模式组合已完成验收。未解决重点仍有短特征回退、固定/自动标签、引物/酶切的系统单参数对照、其他公开入口的人工语义审阅及完整案例 L4。

第二批最终验证：`python3 tools/build-website.py` 成功构建 50 对页面和 95 个帮助主题；`python3 tools/check-documentation.py` 通过双语代码等价、图片指纹、链接/锚点/搜索检查。浏览器确认箭头课程三图、共用参数说明和“通用绘图技巧”导航；390 × 844 视口正文可读，章节抽屉显示新分组并可收起，随后恢复默认视口。未运行包级发布检查、合并、推送或部署。
