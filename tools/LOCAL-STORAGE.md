# 本地生成物与空间维护

包运行使用源码与随包数据。旧安装副本、检查目录、网站预览和临时工具依赖不应
长期堆在 checkout；`.Rbuildignore` 已将这些开发目录排除出发布包。

```sh
# 只盘点，不删除
python3 tools/clean-local-artifacts.py

# 明确清理已核实的忽略生成物
python3 tools/clean-local-artifacts.py --clean
```

工具仅处理本地 `.site`、带检查日志的根级 `*.Rcheck`、旧 OCR 临时依赖，以及
`.codex-tmp/release-v*` / `cran-v*` 中已知的构建、安装、下载和设备输出子目录。
目标必须被 Git 忽略、没有跟踪文件且不是符号链接。清理前后核对所有原跟踪文件
内容与 Git refs；保留精简检查日志和网站来源快照到 `.codex-tmp/cleanup-records`。
恢复本地网站预览需显式重建 `.site`，清理不会部署或删除远端网站。

以下内容保留：`R/`、`data/`、原始 `examples/`、`man/` 及正式引用图像、冻结
`docs/`、设计文档、审计脚本与既有验证材料、个人备份和 Git 分支/标签/历史。
工具不清理整个 `.codex-tmp`，不执行 `git clean`、历史重写或自动垃圾回收。

日常诊断图、失败候选和大布局缓存应输出到系统临时目录或 `.codex-tmp`；需要
保留证据时，优先记录精简指标、来源和生成命令，只提交必要最终图。相同图像
按校验和复用已有文件。文档图仍在原位更新，避免为每次重跑另存整套图片。
旧发布日志、输入和来源记录不应为节省少量空间而一并清空。

本次实际清理与校验见
[v0.27.0-storage-cleanup-results.json](long-term-audit/v0.27.0-storage-cleanup-results.json)。
该记录是清理事实，不是包级发布认证。

## 2026-11-18 的目录与忽略规则审计

盘点基线为本地 `v0.50.0` 的 `c08cf2af`。用 `git ls-files` 分开统计索引、已忽略和其余未跟踪文件，再按工作区文件逻辑字节数汇总。索引有 **2,315 个文件、131.2 MiB**。项目与本地规则已忽略 **3,533 个文件、202.1 MiB**。当时普通未跟踪文件为零，已跟踪且命中现行忽略规则的文件也为零。以下容量以 1 MiB = 1,048,576 字节计算。

| 路径 | 当前状态 | 建议 |
| --- | --- | --- |
| `.site/`、`.codex-tmp/`、`ggchord.Rcheck/` | 已忽略，分别为 115.0、34.0、16.0 MiB | 临时绘图、构建及检查结果继续写到既有临时位置。只读清理工具此次实际列出 `.site/` 和 `ggchord.Rcheck/`，两者文件系统分配合计 134.84 MiB。 |
| `AGENTS.md`、`design-roadmap/` | 由本 checkout 的 `.git/info/exclude` 排除。路线图目录 213 个文件、22.2 MiB | 作为当前工作区的开发规则和设计资料保留。共享版的 v0.50.0 实施记录在 `tools/long-term-audit/`。 |
| `tools/long-term-audit/`、`tools/teaching-audit/` | 已跟踪，分别为 535 个文件、75.7 MiB 和 62 个文件、8.3 MiB | 保留归档基线、配方、指标及被引用的图。新一轮试画和中间快照进入临时目录，选定证据后再提交。 |
| `man/`、`data/`、`data-raw/`、`examples/` | 帮助、文档图、随包数据、重建配方和原始示例均有跟踪文件。特定本地目录或第三方比较文件已有单独规则 | 按路径保留可追溯的输入和输出。按 PNG、CSV、JSON、RDS 或日志后缀统一忽略，会掩盖这些目录里应当提交的新文件。 |
| `docs/`、`website/` | 已跟踪，分别为 278 个文件、6.8 MiB 和 487 个文件、6.0 MiB。前者是冻结归档，后者含正文源、翻译及配套资源 | 按当前网站暂停决定保留，后续网站重构时统一审定生成物的归属。 |

**建议下次调整的规则，按优先级排列。** 已用 `git check-ignore -v --no-index` 校验现状。本次审计保留了现行规则与索引。

1. **补齐 R 会话文件。** 将限定根目录和个别子目录的规则统一为按文件名匹配的 `.Rhistory`、`Rplots.pdf`，再加 `.RData`、`.Ruserdata` 和根目录本机配置 `/.Renviron`。当前 `examples/.Rhistory`、`examples/Rplots.pdf` 及后三种根级文件都未命中。这些是规则探针，此次未发现这五个路径有待提交文件。
2. **收窄通配范围。** 根规则 `.*/` 会连 `.devcontainer/devcontainer.json` 一起隐藏。可明确列出 `/.codex-tmp/`、`/.cache/`、`/.site/`、`.Rproj.user/` 和实际使用的编辑器缓存目录。`vignettes/.gitignore` 的 `*.R` 也会隐藏新的教学辅助脚本。可只列两个教程生成的根级 `/ggchord_vignette.R`、`/ggchord_guide_hans.R` 及对应 HTML，让 `tutorial-assets/` 下的新源码正常可见。[根规则](../.gitignore) · [教程规则](../vignettes/.gitignore)
3. **归档体积优先按内容核对。** 以 `HEAD` 的 Git blob ID 和文件大小比较，`tools/long-term-audit/` 内有 25 组内容完全相同的文件，多出的路径合计占工作区逻辑大小 11.55 MiB。Git 已复用这些相同 blob。若以后归并，应先核对历史报告、脚本和文档对各路径的引用，再决定用哪一处作为来源。[现有渲染报告及证据引用](long-term-audit/v0.50.0-render.md)

**Git 对象库另计。** `du -sh .git` 约为 319 MiB。`git count-objects -vH` 报告三项配套文件不完整，合计 19.63 MiB。一组 `.pack`、`.rev` 缺对应 `.idx`，另一个 `.rev` 缺对应的 `.idx` 和 `.pack`。`git fsck --connectivity-only --no-dangling` 通过。这些文件属于本地对象库，忽略规则只作用于工作区路径。本次只读盘点保留了文件和历史。
