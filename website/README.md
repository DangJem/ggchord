# ggchord 文档网站源 / Documentation source

本地开发版本 0.16.5，API 基线为 0.16.4；**尚未部署**。正式发布版仍是 0.15.0。
冻结的 `docs/` 保留原样。这里不使用旧站生成物作为当前接口来源。

Local development documentation for 0.16.5, based on the 0.16.4 API candidate.
**Not deployed.** The formal release remains 0.15.0. Frozen `docs/` is untouched.

## 2026-11-01 复核状态 / Review status

当前网站是首轮文档原型。用户要求参考成熟 ggplot2 扩展包网站重做风格并补足内容，
复核确认 reference 入口齐全不等于完整任务教程已完成。v0.16.5 已重新打开，当前
缺口与技术建议见[复核记录](../design-roadmap/DOCUMENTATION-REVIEW-v0.16.5.md)。
推荐的 pkgdown 方案尚未替换下方自制生成器；不要把两者误认为同时完成的正式方案。

This is the first documentation prototype. Review confirmed that complete
reference coverage does not establish complete task tutorials. Version 0.16.5
has been reopened for package-documentation styling and deeper content; see the
[review](../design-roadmap/DOCUMENTATION-REVIEW-v0.16.5.md). The recommended pkgdown
migration has not replaced the generator described below.

## 本地生成 / Local build

在仓库根目录准备已安装的 R、Python ≥ 3.9、Pandoc（rmarkdown 可发现）及包开发依赖
`pkgload`、`rmarkdown`、`jsonlite`。无需安装前端框架或 Python 第三方包。

From the repository root, use R, Python ≥ 3.9, Pandoc (discoverable by rmarkdown),
package dependencies, and `pkgload`, `rmarkdown`, `jsonlite`. No frontend framework
or third-party Python modules are needed.

```sh
# Only when changing shared plot code / 修改共享绘图代码后才运行
Rscript tools/documentation-plots.R
# Update paired source documents, without knitting / 同步双语源码，不渲染 vignette
python3 tools/sync-tutorials.py
python3 tools/build-website.py
python3 tools/check-documentation.py
python3 -m http.server 8765 --directory .site --bind 127.0.0.1
```

打开 / Open `http://127.0.0.1:8765/`。英文与汉语首页分别为 `en/index.html`、
`zh/index.html`。生成目录 `.site/` 被 Git 与包构建忽略；不会写入 `docs/`。

## 来源与维护 / Sources and maintenance

- `pages.json`：配对导航和标题 / paired navigation and titles.
- `en/`、`zh/`：相同章节锚点、代码标记的专题 / matched topic anchors and code markers.
- `tools/documentation-examples.R`：共用可运行代码 / shared runnable code.
- `tools/documentation-plots.R`、`assets/plots.csv`：显式尺寸、设备、输入、hash / dimensions, device, inputs, hashes.
- `assets/captions.json`、`terms.json`、`migrations.json`：配对图注、术语、当前迁移 / paired captions, terminology, migration.
- `R/` → `man/` → `tools/render-reference.R`：完整参考，API 正文为英语原文 / complete reference with canonical English API text.
- `tools/sync-tutorials.py`：从网站专题与代码生成 README 和两份 vignette 源码 / source-only synchronization.
- `.site/build-manifest.json`：版本、提交、dirty 状态及源 SHA256 / version, commit, dirty state, source SHA256.

修改专题应同时更新对应语种；代码只改共享来源；帮助只改 roxygen 源码后生成 Rd。
日常不运行完整 vignette、包级检查或发布设备矩阵。网站生成只转换 Markdown/Rd，不执行
reference examples。`_pkgdown.yml` 保留完整分类供未来需要，但本版新站不依赖 pkgdown。

Update both topic languages together; edit code only in its shared source, and help
only through roxygen/Rd. Ordinary development does not knit full vignettes, run
package checks or a release device matrix. The site converts Markdown/Rd without
executing reference examples. `_pkgdown.yml` retains complete fallback categories;
the new site does not depend on pkgdown.

## 部署切换方案 / Deployment cutover

使用 GitHub Pages Actions；`.github/workflows/documentation.yaml` 的普通 CI 只检查源码，
手动运行生成预览 artifact。只有明确授权后才勾选部署，并在仓库 Pages 设置中选择
GitHub Actions、确认 `github-pages` 环境保护和目标源提交。默认 `deploy=false`。
该工作流不会在 push 时部署；无需向 main 写生成 HTML。

Use GitHub Pages Actions. Ordinary CI checks sources; a manual run builds a preview
artifact. After explicit authorization, select deployment, configure the repository
Pages source as GitHub Actions, confirm environment protections and the exact source
commit. `deploy=false` is the default. Push events never deploy generated HTML.

上线前核对预览 artifact、双语路径、旧 URL 映射、资源与源 hash，记录部署 URL、run、
commit。旧站 `docs/` 只在新站上线并完成链接核验后，于独立授权切换中移除。
本轮不更改远端设置、不运行工作流、不 push、不部署、不删除旧站。

Before cutover, inspect the preview artifact, language pairs, legacy URL mappings,
assets and source hashes; record the deployment URL, run and commit. Remove frozen
`docs/` only after successful deployment and link verification in an authorized
cutover. This local task changes no remote settings, starts no workflow, pushes or
deploys nothing, and retains the old site.

方案依据 / Design sources:
[GitHub Pages custom workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages),
[roxygen Markdown](https://roxygen2.r-lib.org/articles/markdown.html),
[ggplot2 grammar](https://ggplot2.tidyverse.org/).

## 图标 / Artwork

用户于 2026-11-01 提供的图标用于首页、README 和教程入口，完整图像为 1254 × 1254。
`man/figures/logo.png` 是保留透明度的无损 PNG，网页优先加载无损 WebP；两者与输入
RGBA 逐像素一致。64 × 64 导航图像是另行注明的缩略图。输入/输出 hash 与体积见
`assets/icon.json`。`tools/optimize-documentation-icon.py` 可从原始 PNG 重建，额外需要
Pillow 和 cwebp；这些仅为资产维护工具，常规网站构建不依赖它们。

The supplied 1254 × 1254 artwork appears on the home page, README and tutorial.
Full PNG and WebP preserve identical decoded RGBA pixels; the 64 × 64 navigation
image is an explicitly resized thumbnail. See `assets/icon.json` for provenance,
hashes and sizes. Re-encoding with `tools/optimize-documentation-icon.py` requires
Pillow and cwebp only for asset maintenance, not for normal site builds.

手动授权部署使用 `--publication` 标识候选文档，仍不宣称 R 包正式发布；构建器自身
不部署。授权构建清单注明应以 Pages run 核对实际部署结果。
Authorized deployment builds use `--publication` to label unreleased candidate
documentation. The builder does not deploy; verify actual deployment through the Pages run.
