# 开发工具入口

从仓库根目录运行命令。先按要回答的问题选脚本，再查看脚本开头的输入、依赖和输出位置。

| 任务 | 入口 | 使用时机 |
| --- | --- | --- |
| 检查一处数据或图层修改 | 直接在 R 中加载当前源码，调用受影响的公开函数。按主题查找 `validate-data-*.R`、`validate-clean-intervals.R`、`validate-collinear*.R` 等已有脚本 | 默认使用最小输入。部分旧脚本有独立的样本、版本或存储要求，先读其开头 |
| 检查文档来源和链接 | `python3 tools/check-documentation.py --source-only`，单查仓库链接用 `check-repository-links.py` | 更新网站源、函数帮助或迁移说明后。输入与生成顺序见[网站维护说明](../website/README.md) |
| 更新双语说明 | `sync-website-languages.py`、`sync-tutorials.py`、`reference-translations.py`、`generate-parameter-atlas.py`；其他 `generate-*.py` 生成相应教学索引 | 对照[双语维护说明](BILINGUAL-MAINTENANCE.md)更新汉语源和保存的译文。Rd 来自 `R/` 中的 roxygen |
| 复核教学示例或图片 | `check-documentation-lessons.R` 执行示例，`documentation-plots.R` 绘制常规网站 PNG。共享示例、课程与综合案例分别见 `documentation-examples.R`、`site-course-examples.R`、`site-flagships.R` | 教学执行回执在 `teaching-audit/`，图片输入与哈希在 `website/assets/plots.csv`。示例输出与真实图片分别核对 |
| 本地构建网站 | `build-website.py`。它使用 `site_*.py`、`render-website.R`、`finalize-website.R` 和专题渲染器。CI 的本地构建依赖见 `website-dependencies.R` | 网页预览在 `.site/`。下载源码须受 Git 跟踪。运行时及缓存输入明细在 `.cache/website-diagnostics/` |
| 查看设备、几何或耗时证据 | [profile-showcase-render.R](profile-showcase-render.R) 分开测量完整出图的构造、build、gtable 和写图时间。`benchmark-layout.R` 和 `profile-foundation-ribbons.R` 用于连接布局 | 选定图、设备与尺寸后再比较。测量结果见[图层接口与性能复审](long-term-audit/v0.29.0-layer-grammar.md)和[二次审计](long-term-audit/v0.29.0-second-review.md)。新过程图和布局先写临时目录，留存规则见[本地存储说明](LOCAL-STORAGE.md) |
| 正式发布 | `validate-release*.R`、`validate-rc*.R`、`render-rc.R`、`benchmark-rc.R`、`check-vignette-resources.py`，以及 `release-audit/`、`cran-audit/` | 按该版本的发布授权与验收清单调用。[默认测试范围](TESTING.md) |
| 继续 v0.50.0 架构重构 | 从[设计总表](long-term-audit/v0.50.0-design.md)进入导出去向、六个使用场景和基线引擎探针 | 先读已确认的直接切换政策及待评审接口，再按阶段记录实现证据 |
| 调查历史接口或实验 | `audit-*.R`、`freeze-api.R`、`catalog-api-migrations.R` 和相应时期的 `validate-api*.R`。版本结论、脚本、测量与图像可从[历史资产索引](long-term-audit/archive-index.tsv)按用途检索 | 历史结论适用于记录的提交与输入。当前开发结果见 [v0.29.0 审计整改记录](long-term-audit/v0.29.0-review.md) |

## 历史证据与依赖

历史索引记录 `main@ab415637` 的 918 个原始路径、文件用途、字节数、SHA-256 和恢复路径。`storage` 标示当前 checkout 或历史 Git 提交。路径以 `tools/long-term-audit/` 为基准。检索某次实验时，先找其报告或测量，再沿文件名核对图片。对字节相同的独立实验保留原名映射。需要完整还原旧工作树时，可以从该基线提交提取原文件。从旧提交恢复一份文件到临时目录可运行 `git show 8e805060:tools/long-term-audit/<文件名> > /tmp/<文件名>`。本轮逐个核对 420 个移出文件的 Git 对象，保留了 496 个有引用或承担来源、许可及索引用途的文件。重新运行实验产生的过程输出请使用临时目录。

部分时期的工具依赖维护者本地的 `design-roadmap/`，该目录未随普通 clone 分发。例如 `audit-documentation-coverage.py` 需要 `design-roadmap/audit-generated/candidate-exports-v0.16.4.csv`。需要复现当期结果时，使用与该脚本版本匹配的本地输入。要查看当前导出和参数，可直接运行 `Rscript tools/teaching-inventory.R`。随包数据的再生成输入见[生成器索引](../data-raw/README.md)。

网站命令保留在 `tools/`。搜索、导航、页面组件和链接转换的共享 Python 实现集中在 `tools/website/`。网站下载的教学证据表由双语兼容页引用，其余教学清单供开发工具读取。

包声明 R ≥ 4.1.0、ggplot2 ≥ 4.0.0。现有常规 CI 在 macOS、Windows 的 R release 以及 Linux 的 release 和 devel 上运行。本次使用 R 4.6.0，在独立库核对 ggplot2 4.0.0、4.0.3 和 3.5.2，具体行为见[依赖下限与绘图核验](long-term-audit/v0.29.0-second-review.md#dependency-floor)。R 4.1.x 与最低可安装依赖组合仍缺少实机验证记录。
