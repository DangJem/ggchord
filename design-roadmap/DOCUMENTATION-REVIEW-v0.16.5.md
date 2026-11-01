# v0.16.5 网站设计与内容复核

日期：2026-11-01。复核基线：`85d7765`。用户认为首轮截图风格不符合预期、内容
不足，要求参考 gggenomes 以及其他 ggplot2 扩展包网站的设计和写作方式。

**结论：用户指出的问题成立。首轮是网站框架与入门内容实现，不能认定为完整文档
重构完成。v0.16.5 重新打开；本文件不表示新版网站已经重做、部署或通过内容验收。**

## 参考网站实际如何构建

以下是官网、源配置和文章页的核查结果；不依据视觉相似性猜生成器。

| 网站 | 已核实技术 | 内容与写作方式 | ggchord 可借鉴之处 |
| --- | --- | --- | --- |
| [gggenomes](https://thackl.github.io/gggenomes/) | [源配置](https://github.com/thackl/gggenomes/blob/main/_pkgdown.yml)使用 pkgdown、Bootstrap 5、Bootswatch Minty | 首页给出真实六条病毒基因组的综合案例；教程从序列逐步增加注释、连接与数值轨道；reference 按任务分组 | 用有来源的综合图说明能力，正文保持宽阔、图和代码占主导，导航有 Get started / Reference / Articles / Changelog |
| [ggplot2](https://ggplot2.tidyverse.org/) | [源配置](https://github.com/tidyverse/ggplot2/blob/main/_pkgdown.yml)使用 pkgdown、Bootstrap 5、tidytemplate | reference 按 grammar 职责分类，专题、开发文档与 FAQ 分层 | 数据映射、图层、Position、scale、guide、coord、theme 各自有说明入口与互链 |
| [patchwork](https://patchwork.data-imaginist.com/) | 官网页脚明确为 pkgdown；没有凭首页推定具体主题 | 首页简短例子说明核心能力，文章按组图、布局、注释、跨页对齐组织；reference 有分组解释 | 面向“我要完成什么”的任务教程，再链接到精确接口 |
| [cowplot](https://wilkelab.org/cowplot/) | 官网页脚明确为 pkgdown | 首页概述与安装，intro 展示主题、排列、注释；详细 articles/reference 分开 | 每项主要能力有解释、代码和对应效果，入口页不承载全部帮助正文 |

gggenomes 的 Minty 是已核实的参考实现，不是用户要求必须复制的主题。ggplot2 的
家族模板不是 ggchord 身份，不应直接复制其品牌标识或额外服务脚本。
本轮只研究这些资料，没有复用其具体代码、数据、图像或字体文件。

## pkgdown 的写作链路

根据 [pkgdown 官方说明](https://pkgdown.r-lib.org/)及
[定制说明](https://pkgdown.r-lib.org/articles/customise.html)：

| 来源 | 生成内容 | 人需要完成的工作 |
| --- | --- | --- |
| README / 专用首页源 | 首页 | 定位、安装、旗舰案例、学习入口 |
| `R/` roxygen → `man/*.Rd` | Reference | 准确签名、参数、返回、约束、例子和交叉引用 |
| `.Rmd` / `.qmd` 文章源 | Get started / Articles | 完整任务叙述、输入、逐步代码、输出图和解释 |
| `NEWS.md` | Changelog | 实际已实现变化与迁移事实 |
| `_pkgdown.yml` | 导航、reference 分类、模板、图像设置 | 编辑信息架构与生成配置 |

这些站点的内容是维护者编写的 R 文档，pkgdown 负责把来源组织成网站；换主题不会
自动补齐缺失教程。gggenomes 的[真实数据教程](https://thackl.github.io/gggenomes/articles/emales.html)
提供预处理数据并说明可选的外部分析步骤，各阶段都有代码、图片和解释。

[语言设置](https://pkgdown.r-lib.org/articles/translations.html)可翻译 pkgdown 生成的界面词，
**不会自动把文章、函数说明或示例解释翻译成另一种语言**。双语正文仍须成对维护。

## ggchord 首轮实际覆盖了什么

`python3 tools/audit-documentation-coverage.py` 只解析源码，不运行示例，生成
[130 项覆盖台账](audit-generated/documentation-coverage-v0.16.5.csv)。实际结果：

- 130 个导出均能对应到 Rd，包含 108 个自有入口与 22 个上游重导出；reference
  入口确实齐全，这项工作保留。
- 共用文档脚本直接调用 38 个自有入口与 6 个上游入口。**这不是“仅实现 38 个
  功能”，也不是要求为每个别名各写一遍教程**；它说明首轮的小型共用示例没有覆盖
  许多主要任务，必须做下面的任务级复核。
- 17 对页面包含首页、安装、FAQ、术语、reference 索引、变更等入口页，并不等于
  17 篇深入教程。部分高级模块只是一段概念加一个基础调用。
- 六张共用图以 A/B/C 合成输入为主。现有原始 FASTA/GFF3/BLAST 被读取验证过，
  但网站没有把它们串成有图、有解释的真实数据主教程。
- 汉语路径的 reference 正文仍是完整英语原文。首轮已对此标注，但不能把“双路径”
  统计成“参考正文已完整双语翻译”。

## 按主要任务复核的缺口

“需补”表示当前网站缺少充分演示，**不表示包没有实现该功能**。源码和已有 tools/
examples/ 可以作为写作输入；它们存在不等于读者已能从网站完成任务。

| 任务/功能簇 | 首轮网站证据 | 复核结论及应补内容 |
| --- | --- | --- |
| 三种模式与基本组成 | 三个最小组合和三张合成图 | 基础有；需说明各场景选择、支持边界及图形解读 |
| 真实多序列比较 | FASTA/GFF3/BLAST 读取片段 | 缺主线；应串起读取、检查、清理、映射、出图与保存，不改原始记录 |
| 单序列生物学圆图 | 三个短合成 feature，另有短序列搜索 | 缺完整场景；需要来源明确的序列、注释、primer、restriction 与标签共存说明 |
| 共线性锚定与方向 | 一个 anchor/auto 例子 | 缺逐步比较与诊断解释；不能把显示方向直接解释为生物学倒位 |
| 清理、准备、层次校验 | validate/clean/prepare 基础调用 | 缺问题输入、报告解读、hierarchy、显式保留/删除策略和来源追踪 |
| 输入格式与适配器 | 运行 FASTA/GFF3/BLAST；其他主要为文字 | 缺 BED/PAF/GenBank、S4 适配、分段/原点及未知链的数据实例 |
| ribbon 准备与统计 | 筛选/去重/合并/聚束/密度主要被提及 | 缺有输入、结果解释的工作流及聚合来源/估计标记说明 |
| 区间连接与点连接 | 共享脚本只有 ribbon | 缺 line、选对、邻接策略、anchor 和两种几何用途比较 |
| feature / gene / plasmid / Position | 基础 gene/feature 和逐链位移 | 缺 feature stack、plasmid preset、映射 shape、同类多层和位置继承实例 |
| 标签策略 | gene auto 图及 feature 默认 | 缺 gene 三种布局及 feature 策略的代表对照、隐藏/适配诊断、固定层替代 |
| 轴、参照条、高亮、中心标签 | 基础源坐标轴；其余没有专门演示 | 缺 region、seq scale、center label、命名参数与坐标支持边界 |
| 定量轨道 | 合成正弦 bar | 缺 line/tile、多个轨道、共享 ID/domain、baseline、focus track、数据单位解读 |
| primer / restriction / common features | 极短序列的 primer/EcoRI | 缺 catalog、显示筛选、通用 feature 识别、已知参考与未知样本的边界 |
| CDS、同源关系与局部工作流 | 一个 feature 投影，focus/segment 基础 | 缺 link 投影、显式 cluster 归属、loci 查找、局部翻转及源位置映射 |
| scale、guide、theme | 一条 ribbon 梯度与基础 theme | 缺命名手动标度、图例组织、预设对照、固定值/映射与用户 scale 优先的实用说明 |
| 输出、完整布局与优化 | preview/diagnose/export 示例；getter/optimize 提及 | 缺字段消费实例与结果解释，需说明 measured 不等于视觉认证及设备差异 |
| Gallery 与教学数据 | 六张最小图与尺寸清单 | 缺场景导向的完整案例、数据出处、可下载代码和参数效果对照 |
| 迁移与 FAQ | 当前迁移表和基础 FAQ | 有可复用基础；需从具体教程互链，并按新任务缺口补充 |

## 修订方向与新验收标准

用户已明确的方向：接近成熟 ggplot2 扩展包的文档站观感，并充分呈现 ggchord 现有
能力。gggenomes 是示例而非唯一模板；不能只更换配色后再次宣称完成。

1. **页面形态**：顶栏放入门、函数参考、专题、图集、变更及搜索/语种。白色或浅色
   正文以文章、图、代码为主；右栏辅助目录/项目资料。首页 logo 适量展示，让位于
   真实综合图；避免把科研包官网设计成密集卡片仪表盘。
2. **内容主线**：先有来源可追踪的真实综合案例，再拆出三模式入门和上表专题。
   每个主要任务给出输入、前置步骤、可复制代码、实际图/结果、解释、边界及接口链接。
3. **覆盖标准**：维护“主要功能簇 → 教程章节 → 可运行来源 → 图片/结果 → reference”
   台账。别名和重导出允许共享解释；低频格式不机械穷举。页面数、链接全绿、函数名
   被提及均不能代替内容验收。
4. **双语**：专题与图注按相同任务成对，代码和资产共用；明确 reference 英语原文与
   汉语任务讲解的关系。若进一步逐条翻译 reference 正文，需明确为额外内容选择，
   不把导航语言设置当作全文翻译。
5. **资产**：保留用户原画和已验证的无损编码。图来自实际 ggchord 代码与有来源输入，
   不复制参考网站图片，不伪造同源关系或生物学数值。
6. **验证边界**：普通开发继续做相关例子的有限验证和静态预览。新增网站方案须保证
   不写冻结 `docs/`，不默认执行整套 vignette/发布检测；远端切换与发布仍另行授权。

**技术推荐，尚未执行切换**：采用 pkgdown + Bootstrap 5 作为正式文档构建基础，
优先使用通用模板并针对 ggchord 定制。其优势是 API 自动链接、reference/文章/NEWS
一体化和生态一致性；本机已安装 pkgdown 2.2.1。具体主题、双语构建路径、仅构建
网站文章而不自动重跑完整 vignette 的方式仍须在实施时确定，不能把本轮建议记为
已经实现。也不把用户的外观示例误写成“用户指定必须使用 Minty”。

## 本轮实际完成与未完成

已完成：四个参考站的首页/内容结构调查、gggenomes/ggplot2 源配置核对、gggenomes
实际页面视觉检查、130 项本地入口/示例来源盘点、主要功能簇内容复核，及版本状态
纠正。首轮参考索引、代码/资产共享和无损图标工作保留。

尚未完成：按上述方向重做前端、补齐主要任务的双语文章与实图、确定并实现构建方案、
完成新内容验收。没有修改 R 的公开 API、冻结 `docs/`、网站部署或 GitHub 发布状态。
