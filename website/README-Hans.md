# 维护 ggchord 文档

[English](README.md) · [现代汉语](README-Hans.md)

网站说明当前检出版本的公开 API。比较本地构建与已部署网站时，先核对包版本和文档通道。[包元数据](../DESCRIPTION) · [文档通道](version.json)

## 更新文章或导航标签

文章正文维护在 `website/en/`。英语标题、导航标签、图注及其他结构化文案维护在 `website/content-en/`。可翻译字段采用 `{"en": "Get started"}` 这样的格式。内容清单为每项列出英语源、保存在 `website/translations/zh/` 的汉语译文及生成文件。渲染器从 `website/zh/` 读取汉语 Markdown，并从双语 JSON 读取结构化文案。[内容清单](content-policy.json)

修改英语原文后，在仓库根目录准备翻译草稿。

```sh
python3 tools/sync-website-languages.py --prepare
```

对照当前英语逐项校对标有 `needs_review` 的条目。草稿记录了英语原文及其 SHA-256 摘要。更新汉语译文，将状态改为 `reviewed`，然后生成并检查输出。

```sh
python3 tools/sync-website-languages.py
python3 tools/sync-website-languages.py --check
```

检查会指出缺译、校对所依据的原文已变，以及直接修改生成文件等问题。两种语言共用代码、指令、图像标识、导航标识、链接目标和显式章节锚点。翻译时应保持标题层级和表格列一致。

## 导入完整的校对稿

更新已有文章时，使用清单中登记的英语源路径，并提供已校对汉语 Markdown 的路径。批量新增页面时，将汉语校对稿放在同一目录，文件名与各自的英语源一致。目录导入命令会先检查所有提交的双语稿，全部通过后再登记新页并写入生成文件。

```sh
python3 tools/sync-website-languages.py --reviewed-page en/installation.md /tmp/installation-reviewed.md
python3 tools/sync-website-languages.py --reviewed-directory /tmp/reviewed-han
```

修改结构化文案时，先编辑英语 JSON。校对稿的结构、英语字符串和共享字段均应与原文一致，每个可翻译字段同时包含 `en` 和 `zh`。新增结构化文件时，先在内容清单中登记，再导入译稿。下列命令用于登记已校对的首页文案。

```sh
python3 tools/sync-website-languages.py --reviewed-content home.json /tmp/home-reviewed.json
```

文章链接指向仓库里实际存在的文件，例如同目录的 `.md` 页面、`../../man/` 下的帮助主题或 `../../man/figures/` 下的图片。渲染时会将这些链接转成网站网址。撤下文章时，一并移除清单条目、英语源、译文记录和生成文件。更新指向该文章的链接，并记录原网址和章节锚点的去向。

## 更新 API 参考

API 参考说明本地源码的公开接口及其分类。先编辑 `website/reference-manual/en.json` 中的英语正文，再校对 `website/reference-manual/zh.json` 中的对应译文，并在汉语文件中登记已校对英语源的 SHA-256。共用的 `structure.json` 保存调用形式、参数与角色名称、示例和关联链接。两种语言共用代码和锚点。原生分面页面链接到 ggplot2 官方文档。

运行下列命令可独立核对和构建 API 参考章节。

```sh
python3 tools/check-documentation.py --reference-only
python3 tools/build-website.py --reference-only
```

本地英语入口为 `.site/en/reference/index.html`，对应汉语入口位于 `.site/zh/reference/`。API 参考沿用工具书的 pkgdown 模板、导航和外观设置。构建器运行页面提供下载的 R 脚本，在双语页面显示相同的图形和控制台输出，并用同一检出版本的包源码刷新 `.site/assets/ggchord-source.zip`。示例运行器使用 downlit、ragg、png 和读取器的可选依赖 rtracklayer、VariantAnnotation、Biopython 及 MUMmer 的 `show-coords`。如果解码器安装在默认搜索路径之外，可以设置 `GGCHORD_PYTHON` 或 `GGCHORD_SHOW_COORDS` 指定位置。

## 更新函数帮助、包 README 或离线教程

- **原生 R 帮助。** 修改英语 roxygen 源并重新生成 Rd。公开接口行为有变化时，核对相应的 API 参考主题。
- **包 README。** 两种语言共用 `website/readme-template.md` 的版式和代码。修改 `website/content-en/readme-text.json`，校对汉语译文并重新生成双语 JSON。运行 `python3 tools/sync-tutorials.py` 更新两份 README，然后加上 `--check` 验证结果。
- **离线教程。** 修改英语源 `vignettes/ggchord_vignette.Rmd`，同步更新汉语版 `vignettes/ggchord_guide_hans.Rmd`。校对后，将英语文件的 SHA-256 写入汉语文件头。README 同步命令也会检查教程摘要、R 代码块、章节锚点及本地资源引用。
- **离线图片。** 运行 `python3 tools/prepare-vignette-figures.py`，为 CRAN 单独分发的教程 HTML 准备嵌入图片。这一步需要 Pillow。发布验收时，用 `tools/check-vignette-resources.py` 核对源码包中教程 HTML 引用的本地资源。

## 构建和检查本地网站

站点面向电脑设计。文章和参考页的全站章节与本页目录在左侧并排。表格和代码置于正文内，过宽的内容可在各自区块内横向滚动。网站验收使用源码、构建和静态 HTML 检查。从构建结果核对导航、样式表、章节链接和示例输入。

基本构建需要 Python 3.9 或更高版本、R、Pandoc、包依赖及示例依赖。网站使用 `pkgdown`、Bootstrap 5、`pkgload`、`rmarkdown` 和 `jsonlite`。完整构建会生成双语页面，回放工具书和 API 示例，并在核对产物后更新本地预览。在仓库根目录执行下列命令。

```sh
python3 tools/sync-website-languages.py --check
python3 tools/check-documentation.py --source-only
python3 tools/build-website.py --render-jobs 2
python3 tools/check-documentation.py
python3 -m http.server 8765 --directory .site --bind 127.0.0.1
```

打开 `http://127.0.0.1:8765/` 可查看 `.site/` 中的文件。两种语言的首页共用四组预览图。点击预览可查看完整图，由 `website/examples/` 下的案例脚本绘制。首页预览沿用同组案例，通过 `Rscript tools/render-site-flagships.R` 重新生成。修改网站、包 README 或离线教程的图像时，保留输入数据、过滤规则、绘图代码及图像来源记录。[原画、第三方资产及许可](ASSET-SOURCES.md)

## 持续集成与部署

推送或拉取请求涉及文档构建输入时，会触发源码检查。手动运行 GitHub Actions 可生成预览产物。部署需要明确授权。`--publication` 控制构建产物的标记。历史网站归档在 `docs/`。包的发布记录见 GitHub Releases，网站的部署记录见 GitHub Actions。
