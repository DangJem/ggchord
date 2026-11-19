# ggchord documentation / 文档网站

The website follows the package API in this checkout. [DESCRIPTION](../DESCRIPTION) identifies the version and [version.json](version.json) records the documentation channel. The deployed website and the current development checkout can have different versions.

## English sources and Han translations / 英语源与汉语译文

Edit English Markdown in `website/en/` and English structured content in `website/content-en/`. Each localizable structured string has the form `{"en": "Get started"}`. Page identifiers, destinations, figure keys, example code and related metadata occur once. The [content policy](content-policy.json) lists the English source, saved translation and generated output for each item. The English R help remains the source of the function reference.

网站以英语为主源。页面正文维护在 `website/en/`，标题、导航、首页文案、图注及 README 文案维护在 `website/content-en/`。汉语译文保存在 `website/translations/zh/`。构建器读取生成的 `website/zh/` Markdown 和配对 JSON。两种语言共用代码、图片、导航标识、链接目标及章节锚点。

1. Edit the English source and run `--prepare`. Changed strings become `needs_review`. Each draft carries the current English text and a SHA-256 fingerprint.
2. Compare the Han translation with the current English. After reviewing an entry, update its translation and change its status to `reviewed`.
3. Generate the outputs, then run `--check`. The site build uses this check to reject missing translations, outdated source fingerprints, changed code, links or anchors, and edits to generated files.

```sh
python3 tools/sync-website-languages.py --prepare
# Review entries in website/translations/zh/
python3 tools/sync-website-languages.py
python3 tools/sync-website-languages.py --check
```

已校对的完整 Markdown 稿可以直接与英语源配对登记。共有代码块与指令应逐字一致，标题层级、表格列及显式锚点应对应。整目录导入会检查每份 `.md` 的同名英语源，全部通过后再更新相应记录，并将新页加入清单。退役页面时，同时移除清单项、英语源、已保存译文和生成页，并处理原网址与锚点。

```sh
python3 tools/sync-website-languages.py --reviewed-page en/installation.md /tmp/installation-reviewed.md
python3 tools/sync-website-languages.py --reviewed-directory /tmp/reviewed-han
```

For structured changes, edit the corresponding English JSON first. A complete reviewed draft must match its structure, English text and shared metadata. Localized fields in that draft hold both `en` and `zh`. Register any new structured file in the policy before importing it. Ordinary cell-level translation can use the same `--prepare` workflow as Markdown.

```sh
python3 tools/sync-website-languages.py --reviewed-content home.json /tmp/home-reviewed.json
```

最初切换以分支提交 `85b1b065` 的 94 个页面、2 个 README 附页和 14 个结构化文件为对照，先验证旧生成器能复现原英语，再以其为新源逆向比对已有汉语。旧稿中未能逐项匹配英语内联标识的段落保留待校对状态。后续版本使用当前清单和英语摘要。

## Source links, reference and README / 链接、函数参考与 README

- Article links use real repository files such as a sibling `.md`, a `../../man/*.Rd` help topic or a figure in `../../man/figures/`. The source-link resolver converts these to site destinations during rendering. Keep explicit section anchors when moving content.
- The reference renders English Rd from the package and checks the saved Han translation in `website/reference-zh/` against each source digest. In `_pkgdown.yml`, each reference group has its canonical English `title` and `desc`, corresponding `title_zh` and `desc_zh`, and the SHA-256 of the English title, one newline and the description. Review both localized fields when English changes.
- The two READMEs share layout and code in `website/readme-template.md`. Edit `website/content-en/readme-text.json`, review the matching Han strings, regenerate the paired JSON and run `python3 tools/sync-tutorials.py`. Run it with `--check` to verify the result.
- The installed offline tutorial has its own English source in `vignettes/ggchord_vignette.Rmd` and reviewed Han edition in `vignettes/ggchord_guide_hans.Rmd`. The Han source header records the reviewed English file's SHA-256. After reviewing an English edit, update the Han source and this digest. The sync command checks the digest, R chunks, section anchors and local asset references.

## Local website build / 本地构建

The site uses pkgdown, Bootstrap 5 and Pandoc. Install the package dependencies, `pkgdown`, `pkgload`, `rmarkdown`, `jsonlite` and the example-specific dependencies. Python ≥ 3.9 supports the core build scripts. Optional image maintenance uses Pillow. The renderer reads the Markdown sources and renders native Rd reference pages with evaluated examples.

```sh
python3 tools/sync-website-languages.py --check
python3 tools/check-documentation.py --source-only
python3 tools/build-website.py --render-jobs 2
python3 tools/check-documentation.py
python3 -m http.server 8765 --directory .site --bind 127.0.0.1
```

网站生成至 `.site/`。首页使用四组旗舰预览，两种语言引用同一份图像。原画、第三方资产与授权见[来源记录](ASSET-SOURCES.md)。README、网站及离线教程涉及的原始数据、过滤规则、绘图代码和图像来源须可追溯。

离线图像维护由 `python3 tools/prepare-vignette-figures.py` 单独执行。离线教程在 CRAN 单独提供 HTML 时仍可使用已嵌入的图片。发布验收时，`tools/check-vignette-resources.py` 核对源码包内 HTML 的资源闭包。

## CI and deployment / CI 与部署

Ordinary pushes and pull requests run source checks. Manual Actions runs can build a preview artifact. Deployment requires explicit authorization. The `--publication` build flag controls the candidate label. The old `docs/` tree is a historical archive. Release and deployment status are recorded separately in GitHub Releases and Actions.
