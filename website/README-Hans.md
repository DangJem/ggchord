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

## 更新函数帮助、包 README 或离线教程

- **函数参考。** 修改英语 roxygen 源并重新生成 Rd。汉语参考使用 `website/reference-zh/` 下的译文，按主题核对英语原文的摘要。`_pkgdown.yml` 中的英语 `title` 或 `desc` 变动时，校对 `title_zh` 和 `desc_zh`。分组摘要为英语标题、一个换行符及英语说明依次拼接后的 SHA-256。
- **包 README。** 两种语言共用 `website/readme-template.md` 的版式和代码。修改 `website/content-en/readme-text.json`，校对汉语译文并重新生成双语 JSON。运行 `python3 tools/sync-tutorials.py` 更新两份 README，然后加上 `--check` 验证结果。
- **离线教程。** 修改英语源 `vignettes/ggchord_vignette.Rmd`，同步更新汉语版 `vignettes/ggchord_guide_hans.Rmd`。校对后，将英语文件的 SHA-256 写入汉语文件头。README 同步命令也会检查教程摘要、R 代码块、章节锚点及本地资源引用。
- **离线图片。** 运行 `python3 tools/prepare-vignette-figures.py`，为 CRAN 单独分发的教程 HTML 准备嵌入图片。这一步需要 Pillow。发布验收时，用 `tools/check-vignette-resources.py` 核对源码包中教程 HTML 引用的本地资源。

## 构建和检查本地网站

网站默认面向电脑设计。核对可用空间时采用 1280、1366、1440 和 1920 CSS 像素。文章和函数参考页的章节列表与本页目录在左侧并排。页面没有本页目录时，相应空间交给正文。表格和代码置于正文内，过宽的内容可在各自区块中横向滚动。逐幅写明示例应展示的图形、标签或对照，核对输入记录，再按桌面实际显示宽度查看结果。检查表格时比对末列和外框，遇到确实较宽的表再验证局部滚动。布局推算与浏览器实测分别记录。手机和平板布局不列入默认设计和验收范围。

基本构建需要 Python 3.9 或更高版本、R、Pandoc、包依赖及示例依赖。网站使用 `pkgdown`、Bootstrap 5、`pkgload`、`rmarkdown` 和 `jsonlite`。构建时会渲染文章，并运行原生 Rd 示例。在仓库根目录执行下列命令。

```sh
python3 tools/sync-website-languages.py --check
python3 tools/check-documentation.py --source-only
python3 tools/build-website.py --render-jobs 2
python3 tools/check-documentation.py
python3 -m http.server 8765 --directory .site --bind 127.0.0.1
```

打开 `http://127.0.0.1:8765/` 可查看 `.site/` 中的文件。两种语言的首页共用四组预览图。修改网站、包 README 或离线教程的图像时，保留输入数据、过滤规则、绘图代码及图像来源记录。[原画、第三方资产及许可](ASSET-SOURCES.md)

## 持续集成与部署

推送或拉取请求涉及文档构建输入时，会触发源码检查。手动运行 GitHub Actions 可生成预览产物。部署需要明确授权。`--publication` 控制构建产物的标记。历史网站归档在 `docs/`。包的发布记录见 GitHub Releases，网站的部署记录见 GitHub Actions。
