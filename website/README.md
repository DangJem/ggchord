# ggchord 文档网站源 / Documentation source

本地开发版本 0.16.5，API 基线为 0.16.4；**尚未部署**。正式发布版仍是 0.15.0。
冻结的 `docs/` 保持原样。网站生成物不作为当前 R 接口的事实来源。

Local development documentation for 0.16.5, based on the 0.16.4 candidate API.
**Not deployed.** The formal release remains 0.15.0. Frozen `docs/` is untouched.

## 生成基础 / Generator

唯一正式生成基础为 **pkgdown + Bootstrap 5**：顶栏、文章布局、主题切换及完整
函数参考使用 pkgdown 的模板和公开生成接口。英语与现代汉语分别位于 `en/` 和
`zh/`；所有正文、专题、图集和函数页均有同页语言切换，保留对应锚点。手机端直接
显示切换入口。根入口记忆上一次切换的语言，无法使用存储时仍可正常跳转。

The sole generator is **pkgdown + Bootstrap 5**. Native templates provide the
navigation, article layout, theme control and complete function reference.
English and modern Han pages live in `en/` and `zh/`; each page switches to its
counterpart, preserving matching anchors, including on mobile. The root entry
remembers the chosen language when browser storage is available.

成对 Markdown 经 Pandoc 转换，展示已有代码、结果与图，不执行其中的 R。完整
reference 从当前 Rd 生成，不运行帮助示例。帮助正文保留完整英语原文；汉语页面
提供本地化导航、分类和任务说明，不声称全文帮助已翻译。首次原型的自制页面壳和
Rd 转换器已退役，不同时维护另一套网站实现。

Paired Markdown is converted by Pandoc using existing code, outputs and figures;
R chunks are not executed. Current Rd generates the complete reference without
running examples. Canonical help remains English with localized Han navigation
and workflow framing. The first prototype's custom page shell/Rd renderer is retired.

## 本地生成 / Local build

需要 R、Python ≥ 3.9、Pandoc（rmarkdown 可发现）、包依赖及 `pkgdown`、`pkgload`、
`rmarkdown`、`jsonlite`。本轮使用 pkgdown 2.2.1。无 Python 第三方构建依赖。
系统字体；依赖资产由 pkgdown 复制到产物，界面不依赖 CDN。构建禁用可选的远端
favicon 服务，直接使用本地图标，禁用仓库外 Sass 缓存。

Requires R, Python ≥ 3.9, Pandoc (discoverable by rmarkdown), package dependencies,
`pkgdown`, `pkgload`, `rmarkdown`, `jsonlite`. Developed with pkgdown 2.2.1. No
third-party Python modules. System fonts and locally copied dependencies avoid
CDN UI dependencies. The optional remote favicon service and external Sass cache
are disabled; the existing local icon is used.

```sh
# 修改共享例子后才更新 PNG / Only when changing shared examples
Rscript tools/documentation-plots.R
# 仅同步源码，不渲染 vignette / Sync sources, without knitting
python3 tools/sync-tutorials.py
python3 tools/build-website.py
python3 tools/check-documentation.py
python3 -m http.server 8765 --directory .site --bind 127.0.0.1
```

打开 / Open `http://127.0.0.1:8765/`。产物位于忽略的 `.site/`，不会写入 `docs/`。
这不是 `pkgdown::build_site()` 的完整 vignette 构建，也不是 R 包发布检查。

## 来源与维护 / Sources and maintenance

| 来源 / Source | 职责 / Role |
| --- | --- |
| `_pkgdown.yml` | Bootstrap、导航与完整 Rd 分类 / Template, navigation and complete reference categories |
| `website/{en,zh}/*.md`、`pages.json` | 配对文章、显式锚点和标题 / Paired articles, matching anchors and titles |
| `tools/documentation-examples.R` | 唯一共享可运行代码 / Shared executable source |
| `tools/documentation-plots.R` | 图形尺寸、结果和输入声明 / Dimensions, outputs and input declarations |
| `man/figures/*.png`、`assets/plots.csv` | 两种语言及 README/vignette 共用的图片字节和 hash / Shared figure bytes and hashes |
| `coverage.json` | 18 个功能簇的章节、代码及输出证据 / Article/code/output evidence for 18 capability families |
| `tools/sync-tutorials.py` | 从同源同步 README、迁移和双语 Rmd / Synchronize README, migration guides and paired Rmd sources |
| `tools/render-website.R` | pkgdown 原生页面与 Rd 参考生成 / Native pages and Rd reference |
| `tools/finalize-website.R`、`assets/site.*` | 同页语言切换、代码折叠和本地全文搜索 / Language switching, code disclosure and local full-text search |
| `tools/check-documentation.py` | 配对、代码、图 hash、功能簇、链接、锚点和搜索覆盖 / Pairing, code, figures, capability evidence, links, anchors and search |

真实案例使用仓库原始 FASTA/GFF3/BLAST；选择显示记录的规则写入代码和正文，原表
不变。质粒来自随包 pUC19 原始序列。其余小表、轨道数值与 cluster 归属明确为合成
教学输入。构建产物提供原始输入 ZIP、SHA256 清单、共享 R 脚本和构建来源清单。
可选 Biopython/GenomicRanges 例子在环境已有依赖时才执行，结果注明 PASS 或 SKIP。

Real workflows use original repository FASTA/GFF3/BLAST records with explicit
display filters; source tables remain unchanged. The plasmid uses source-backed
pUC19. Other small tables, track values and membership are declared synthetic.
Builds include an original-input ZIP, SHA256 list, shared script and provenance
manifest. Optional Biopython/GenomicRanges examples report PASS or SKIP.

原画无损优化、第三方模板与资产来源见[来源记录](ASSET-SOURCES.md)。用户原画用于
首页、README 与教程入口，不替换为重新生成的图；`assets/icon.json` 保留逐像素
RGBA 校验和体积记录。Pillow/cwebp 仅用于资产维护，不是构建依赖。

See [asset sources](ASSET-SOURCES.md) for the lossless user artwork and third-party
licenses. `assets/icon.json` records exact RGBA identity, sizes and provenance.
Pillow/cwebp are optional asset-maintenance tools, not build dependencies.

## CI 与上线边界 / CI and deployment

普通 push/PR 只检查已跟踪源码。手动 workflow 构建预览 artifact，默认不部署；只有
明确获得网站切换授权后才启用部署输入。`--publication` 只改变候选文档的预览标记，
构建器本身不部署、不推送、不发 Release。旧 URL 有精确对应则重定向，否则明确
显示历史文档提示，不凭空映射到新函数。发布与旧站删除仍需独立授权。

Ordinary push/PR checks tracked sources. Manual workflows produce preview
artifacts and default to no deployment. Explicit cutover authorization is needed
before enabling deployment. `--publication` only changes preview labeling;
the builder never deploys, pushes or releases. Exact historical counterparts
redirect; unmatched URLs receive a historical notice rather than invented aliases.
