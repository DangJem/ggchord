# ggchord 文档网站源 / Documentation source

本分支文档对应的包版本见 [DESCRIPTION](../DESCRIPTION)，稳定版和开发通道见 [version.json](version.json)。正式发布与部署记录分别见 GitHub Releases 和 Actions。
旧 `docs/` 保留为历史归档，不再作为 Pages 部署源；历史地址由版本化清单生成跳转。网站生成物不作为当前 R 接口的事实来源。

The tutorials use the public API in this checkout. DESCRIPTION identifies the package version and version.json records the documentation channel and stable release.
Consult GitHub Releases/Actions for publication and deployment status. Legacy `docs/` is retained as an archive; Pages uses Actions and versioned historical redirects.

## 单一正文源与翻译

网站正文唯一手工源为 `website/zh/*.md`，英语文件 `website/en/*.md` 是生成物。
页面标题、导航、首页和图注等结构化文案编辑 `website/content-zh/`；其余同名配对
JSON由生成器输出。英语译文保存在 `website/translations/en/`，以汉语摘要绑定并
标记校对状态。[生成清单](content-policy.json)按 `source`、`translations`、`output` 列出每组路径。例如更新首页应编辑 `website/content-zh/home.json` 和 `website/translations/en/content/home.json`，随后由同步命令写入 `website/home.json`。

修改汉语正文后，先运行 `--prepare` 标出待翻译段，校对对应英语和保护占位符，再把状态改为 `reviewed`。摘要由准备步骤更新，未经校对的旧译文不会自动通过：

```sh
python3 tools/sync-website-languages.py --prepare
# 校对 website/translations/en/ 中 needs_review 的译文后再生成
python3 tools/sync-website-languages.py
python3 tools/sync-website-languages.py --check
```

构建前强制检查：缺译、汉语摘要过期、保护的代码/链接/锚点被改动，或生成文件
被独立修改，都会失败。翻译在 Markdown/结构化内容阶段完成，不翻译HTML。共享
段落骨架、代码块、指令、资产和导航标识。已有英语译文迁移保存，内容差异先校对。
参数图谱由当前API/契约生成汉语源，再走同一英语生成流程；R包英语Rd及其汉语
函数参考仍遵循既有帮助源和摘要译文机制，不能当作可独立编辑的网站正文。

正文链接使用仓库中真实存在的文件：同语言页面用 `.md`，函数帮助用
`../../man/*.Rd`，科学图像用 `../../man/figures/`；章节保留明确的 HTML 锚点。
`tools/website/source_links.py` 只在展开正文时将它们转换为网站的 `.html` 和资产路径，
因此 GitHub 源码阅读与生成网站共用一份正文。不要在源页引用只存在于构建目录的文件。
源码检查同时核对本地目标、跟踪状态与锚点；Release 的版本链接按完整 tag 引用检查，
不能用同名分支代替标签。

Prose links point to real repository files: `.md` pages, `../../man/*.Rd` help,
and `../../man/figures/` images, with explicit HTML section anchors. Expansion
converts these to site URLs through `tools/website/source_links.py`; build-only
files must not be linked from source pages. Source checks cover tracked targets
and anchors. Versioned Release links are checked against explicit tags rather
than same-named branches.

`python3 tools/check-repository-links.py` 同时检查当前与归档文档的真实文件目标、
锚点及 GitHub 提交引用；共享 README 模板按两种语言的生成位置验证。
`--report /tmp/ggchord-links.json` 可列出待联网复核的外链。联网检查应区分 404、
登录跳转、访问拦截及临时网络错误，不改写历史观测记录或测试用假地址。

The repository link check covers current and archived documentation; the shared
README template is validated through both generated outputs. Its optional JSON
report lists external URLs for separate online review. Distinguish missing
pages from login redirects, access restrictions and network errors, and preserve
historical observations and deliberately invalid test URLs.

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

汉语正文与保存译文生成的英语 Markdown 经 Pandoc 转换，展示已有代码、结果与图，不执行其中的 R。完整
reference 从当前 Rd 生成，同时运行帮助示例。帮助正文从英语 Rd 生成，汉语函数参考逐项翻译并绑定源摘要，缺译或过期会阻止构建。首次原型的自制页面壳和
Rd 转换器已退役，不同时维护另一套网站实现。

Han-source Markdown and generated English are converted by Pandoc using existing code, outputs and figures;
R chunks are not executed. Current Rd generates the complete reference with evaluated examples. Canonical R help remains English; the Han website reference translates every prose section and checks source hashes. The first prototype's custom page shell/Rd renderer is retired.

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
# 检查双语网站源 / Check the paired website sources
python3 tools/sync-website-languages.py --check
python3 tools/build-website.py
python3 tools/check-documentation.py
python3 tools/website/copy_examples.py
python3 -m http.server 8765 --directory .site --bind 127.0.0.1
```

打开 / Open `http://127.0.0.1:8765/`。产物位于忽略的 `.site/`，不会写入 `docs/`。
这不是 `pkgdown::build_site()` 的完整 vignette 构建，也不是 R 包发布检查。

## 来源与维护 / Sources and maintenance

| 来源 / Source | 职责 / Role |
| --- | --- |
| `_pkgdown.yml` | Bootstrap、导航与完整 Rd 分类 / Template, navigation and complete reference categories |
| `website/zh/*.md`、`content-zh/` | 唯一手工正文与结构化源 / Canonical Han prose and shared structure |
| `website/translations/en/`、`tools/sync-website-languages.py` | 校对译文、生成英语和配对JSON / Reviewed translations and generated English/paired JSON |
| `tools/documentation-examples.R` | 唯一共享可运行代码 / Shared executable source |
| `tools/documentation-plots.R` | 图形尺寸、结果和输入声明 / Dimensions, outputs and input declarations |
| `man/figures/*.png`、`assets/plots.csv` | 两种语言及 README/vignette 共用的图片字节和 hash / Shared figure bytes and hashes |
| `coverage.json` | 18 个功能簇的章节、代码及输出证据 / Article/code/output evidence for 18 capability families |
| `tools/sync-tutorials.py` | 从同源同步 README 和双语 Rmd / Synchronize README and paired Rmd sources |
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

## 工具书内容与验证 / Handbook content and verification

汉语正文和已校对英语译文共用导航、代码、图片与锚点。顶栏的五个入口面向开始使用、绘图、数据与操作、案例和函数查询。首页保留四组旗舰预览，图库按生成图像的坐标类型分组。默认显示代表图，可展开整组并使用历史图片锚点。函数索引和参数表取自本分支的帮助与公开签名。

Both editions share navigation, code, figures and anchors. Five entries lead to getting started, plotting guides, data and operations, the gallery and the API reference. The dashboard retains four families of flagship previews. The gallery groups figures by coordinate system and preserves individual figure anchors.

普通绘图示例使用安装包内的数据和公开接口。单个复制按钮、可下载的脚本及网页上的代码来自同一段源。定向验收从生成的 HTML 读取每个独立示例，分别在空工作目录与新 R 进程运行。复杂案例在对应条目保留输入、尺寸和完整图层。文件导入条目在使用处指出依赖和输入文件。

Each standalone example uses the same source for its on-page code, copy button and download. The clipboard check reads the built HTML and executes the examples in fresh R processes and empty working directories. Cases record their input, layers and output size. File import entries identify their exact files and optional packages.

```sh
python3 tools/build-website.py --render-jobs 2
python3 tools/check-documentation.py
python3 tools/website/copy_examples.py
```

页面内容从汉语源和译文生成。修改网站显示图后，按对应生成器重绘并更新图像清单。网页构建输出在本地 `.site/`；离线教程图片有独立的同步命令。

## 离线教程图片

`python3 tools/prepare-vignette-figures.py`（需要开发工具Pillow）从当前PNG生成
共享阅读资源到vignettes/tutorial-assets；原始PNG与网站输出不覆盖。先运行sync-tutorials.py，
再准备图片。两份 html_vignette 将图片、脚本和样式一起内嵌，使 CRAN 单独提供 HTML 时
也能显示图片。网站保留指向独立图文件的放大链接。发布前检查两个 HTML 的资源闭包。

源码包生成后运行 `python3 tools/check-vignette-resources.py ggchord_版本.tar.gz`，
检查两个安装HTML中的全部本地src/href目标，而不只检查源Rmd路径。

汉语代码注释由 r-comments-zh.json 集中维护。汉语下载脚本使用 -zh.R 文件名；检查器比较可执行语句，并逐语言核对实际渲染代码。
