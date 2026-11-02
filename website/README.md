# ggchord 文档网站源 / Documentation source

当前文档源为0.16.6发布候选；API候选基线为0.16.4。实际发布与部署状态以 GitHub Release/Actions 为准。
旧 `docs/` 保留为历史归档，不再作为 Pages 部署源；历史地址由版本化清单生成跳转。网站生成物不作为当前 R 接口的事实来源。

Current documentation sources track the 0.16.6 release candidate, based on the 0.16.4 candidate API.
Consult GitHub Releases/Actions for publication and deployment status. Legacy `docs/` is retained as an archive; Pages uses Actions and versioned historical redirects.

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
reference 从当前 Rd 生成，不运行帮助示例。帮助正文从英语 Rd 生成，汉语函数参考逐项翻译并绑定源摘要，缺译或过期会阻止构建。首次原型的自制页面壳和
Rd 转换器已退役，不同时维护另一套网站实现。

Paired Markdown is converted by Pandoc using existing code, outputs and figures;
R chunks are not executed. Current Rd generates the complete reference without
running examples. Canonical R help remains English; the Han website reference translates every prose section and checks source hashes. The first prototype's custom page shell/Rd renderer is retired.

## 本地生成 / Local build

需要 R、Python ≥ 3.9、Pandoc（rmarkdown 可发现）、包依赖及 `pkgdown`、`pkgload`、
`rmarkdown`、`jsonlite`、示例所需的 `dplyr`。本轮使用 pkgdown 2.2.1。无 Python 第三方构建依赖。
系统字体；依赖资产由 pkgdown 复制到产物，界面不依赖 CDN。构建禁用可选的远端
favicon 服务，直接使用本地图标，禁用仓库外 Sass 缓存。

Requires R, Python ≥ 3.9, Pandoc (discoverable by rmarkdown), package dependencies,
`pkgdown`, `pkgload`, `rmarkdown`, `jsonlite`, and `dplyr` for examples. Developed with pkgdown 2.2.1. No
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
| `website/reference-zh/*.json`、`tools/reference-translations.py` | 完整汉语 Rd 说明与源摘要校验 / Complete Han prose with source-hash validation |
| `tools/finalize-website.R`、`assets/site.*` | 同页语言切换、代码折叠、全文搜索及图片放大 / Language switching, code disclosure, search and image viewing |
| `tools/check-documentation.py` | 配对、代码、图 hash、功能簇、链接、锚点和搜索覆盖 / Pairing, code, figures, capability evidence, links, anchors and search |

真实案例使用仓库原始 FASTA/GFF3/BLAST；选择显示记录的规则写入代码和正文，原表
不变。质粒来自随包 pUC19 原始序列。小型教程也使用随包记录，GC 由实际 DNA 计算，同源分组核对真实编码序列。构建产物提供原始输入 ZIP、SHA256 清单、共享 R 脚本和构建来源清单。
可选 Biopython/GenomicRanges 例子在环境已有依赖时才执行，结果注明 PASS 或 SKIP。

Real workflows use original repository FASTA/GFF3/BLAST records with explicit
display filters; source tables remain unchanged. The plasmid uses source-backed
pUC19. Other examples also use bundled records; GC is calculated from DNA and membership checks actual coding sequences.
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

汉语 Rd 接入使用 pkgdown 2.2.1 的 classed Rd 解析适配；该单点通过显式类检查和生成页
正文/代码检查保护。升级 pkgdown 时需复核解析适配，其余生成使用公开接口。

The Han Rd adapter uses pkgdown 2.2.1 classed parsing with explicit structural
guards. Recheck this one adapter when upgrading pkgdown; page generation uses
public builders. Generated reference prose and code are checked.

## 阅读架构与局部验证 / Reading structure and scoped verification

`navigation.json` 将入门、完整工作流、任务指南及参考分组。桌面保留章节导航，宽屏显示
右侧页内目录；中等及窄屏使用可展开目录。文章的 `prepare` 标记将该页实际依赖的共享
代码就地展开，避免跳回首页寻找变量。首页按科研问题选择三类真实案例，图集每张图
链接其准备、代码及解释。调研依据与取舍见
[跨产品调研](../design-roadmap/WEBSITE-RESEARCH-v0.16.5.md)。

`navigation.json` groups learning, complete workflows, task guides and reference.
A persistent chapter sidebar and responsive page outlines support long articles.
Page-local `prepare` markers expand their actual shared-code prerequisites.
The home page selects real cases by scientific question; every gallery figure
links to preparation, code and interpretation. Research and decisions are recorded
in the linked design document.

```sh
# 只检查文档任务的独立准备，不运行包检查 / Independent lesson objects only
Rscript tools/check-documentation-lessons.R
```

检查 20 页的共享准备与步骤，使用独立环境，不保存用户输出、不打开预览，也不执行
完整 testthat、vignette 或包发布检查。准备代码默认折叠，展开后可按页顺序复制运行。

The scoped script checks 20 pages in separate environments, without user exports,
interactive previews, test suites, vignette rendering or package release checks.
Preparation is initially folded and can be expanded to copy in page order.

折叠页内目录在生成时直接读取实际标题，完整参考页无需等待脚本填充目录。
图片使用资产清单的固有宽高预留比例，懒加载不挤动章节。文档顶栏保持可见，
窄屏移除模板行负边距；数据图形仍按页面宽度缩放，并支持原始尺寸画布。

Compact outlines are generated from actual headings, including reference topics.
Intrinsic figure dimensions reserve their aspect ratios before lazy loading.
The documentation header stays visible; narrow layouts remove inherited negative
row margins. Figures remain responsive with an actual-size canvas available.

代码复制按钮只复制 `pre code` 正文，避免把折叠标题一并写入剪贴板；反馈随页面语种显示。
Copy controls use only `pre code` text, excluding disclosure labels, with localized feedback.

当前示例优先使用 dplyr 动词和 R ≥4.1 的原生 `|>` 管道；保留清晰的向量计算与文件 I/O。
安装示例依赖时不需要加载整个 tidyverse。共享例子改动后应实际重建相应 PNG，
不得仅更新代码摘要来复用未验证的图片。

Current examples prefer dplyr verbs and the native R ≥4.1 `|>` pipe. Keep direct
vector calculations and file I/O where clearer; loading the entire tidyverse is
unnecessary. Regenerate PNGs after shared example changes rather than changing
code hashes without verified exports.

## 当前导航与展示 / Current navigation and showcases

首页与README共用三张旗舰PNG，完整示例源仍只有documentation-examples.R。
首页提供三张带图的模块入口和完整大图；README将长代码折叠保留，入门首图继续使用
最小数据。桌面左侧是本页目录，右侧是全站目录；小屏分别折叠，不混成一个菜单。

Home and README share three flagship PNGs from the same executable source.
Home offers three illustrated task entries and full-size figures; README keeps
long runnable code in disclosures. First-plot lessons remain minimal. The left
outline describes the current page; the right guide covers the whole site.
They become separate disclosures on compact screens.

website/navigation.json是全站任务分类的单一来源。顶部“可视化”菜单读取其中的
弦图/共线性/质粒模块，其他顶栏链接为首页、入门、教程总览和函数参考。全站目录
列出各模块下的具体流程和通用技巧，不在顶栏重复全部页面。

website/navigation.json owns the task taxonomy. The top Visualizations menu
uses its chord/collinear/plasmid modules; Home, Get started, All guides and
Reference are primary entry points. The right guide expands these paths into
workflows and shared techniques rather than repeating every page in the header.

综合示例使用dplyr、purrr、tidyr整理原始记录、展开GC窗口和逐窗口计算；这些包是
Suggests中的可选文档依赖。密集pETDuet-1示例使用显式视野与label_fit="shrink"，
按16.5×8英寸通过标准ggsave保存；自动view预览可以选择不同尺寸，不能混用其布局。

Complete examples use dplyr, purrr and tidyr to prepare source records, expand
GC windows and calculate per-window values. These are optional Suggests for
documentation. The dense pETDuet-1 example uses an explicit viewport and
label_fit="shrink", saved with standard ggsave at 16.5×8 inches; automatic view
previews can choose different dimensions and should not be mixed with that layout.

三个模块从数据表、导入、骨架与图层逐步进入完整示例，顺序为弦图、共线性、质粒。13个独立质粒参考均由当前代码生成；共同数据/图层的coord替换示例与圆形专属参数限制分开说明。

文档图像manifest保留每张PNG及实际生成时共享脚本的MD5，同时记录当前依赖摘要。R实现或共同准备变动会要求重新生成所有受影响图；已审阅的独立共线性/三种标签示例块可局部重绘。只有完整当前输出和实际PNG哈希全部匹配时，才允许为旧manifest补录依赖摘要，不能用补写哈希代替渲染。

## 离线教程图片

`python3 tools/prepare-vignette-figures.py`（需要开发工具Pillow）从当前PNG生成
共享阅读资源到vignettes/tutorial-assets；原始PNG与网站输出不覆盖。先运行sync-tutorials.py，
再准备图片。两份html_vignette内嵌脚本与样式，图像用data-external共用本地资源，放大链接也
使用同一路径。发布前须核对源码包内inst/doc/tutorial-assets及两个HTML的实际资源闭包。

源码包生成后运行 `python3 tools/check-vignette-resources.py ggchord_版本.tar.gz`，
检查两个安装HTML中的全部本地src/href目标，而不只检查源Rmd路径。
