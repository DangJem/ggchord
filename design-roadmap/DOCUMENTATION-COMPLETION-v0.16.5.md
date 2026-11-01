# v0.16.5 本地开发完成记录

日期：2026-11-01；开始于 2026-10-31。基线 `24930b0`，实现提交 `c78723e`。
本版全部本地开发内容及相关轻量验证完成，版本分支为 `codex/v0.16.5`，按既定规则
整合本地 main。未 push、tag、GitHub Release 或部署；未自动启动 v0.16.6。
正式发布版仍是 v0.15.0，本地完成不等于正式发布验收。

## 交付

- F23 已解决：保留包名，DESCRIPTION Title、README、constructor 帮助和首页均涵盖
  chord、single circular 和 collinear。未新增、删除或改名公开函数/参数。
- [网站源与维护说明](../website/README.md)：17 组一一对应的英语/汉语页面，按
  安装、数据、三模式、常用调整、轨道、标签、识别、投影、局部视图、输出、FAQ、
  迁移、图集、术语、reference、开发变更组织。当前 93 个 Rd 主题在两种导航下均
  完整提供，包含 130 个导出及其别名、上游重导出、参数、返回和例子。
- API 参考正文保留完整英语原文并显式标明；不是声称逐条翻译了全部 R 帮助。
  汉语专题、教程、迁移、术语及关键工作流与英语对应，章节锚点与可执行代码同步。
- 两份简洁 README、两份完整 vignette 源、[English 迁移](../MIGRATION.md)和
  [汉语迁移](../MIGRATION-Hans.md)由共享来源生成。迁移文件可直接在仓库阅读，
  历史设计记录不覆盖，旧调用只用于迁移解释。
- `Roxygen: list(markdown = TRUE)` 使现有代码名、链接正确转成 Rd；帮助从源码生成。
  `R/ggchord.R` 只有注释变更，NAMESPACE、公开签名、默认和返回字段保持候选基线。
  `_pkgdown.yml` 补全全部 Rd 分类，但不借此重建冻结网站。
- [共享 R 例子](../tools/documentation-examples.R)、单一 PNG 导出脚本与
  [资产清单](../website/assets/plots.csv)覆盖六张图片；README、网站、vignette、
  两种语言共用同一资产。教学输入明确为合成值，不充当生物学观测。
  圆图例子显式隐藏轴以避免与 feature 标签争用空间，包默认样式不变。
- 清理 10 张当前源文档不再引用的旧 PNG（1,628,163 字节），来源、hash 和冻结
  站点副本记录在[退役清单](../website/assets/retired-plots.csv)。原始业务示例不改写，
  冻结 `docs/` 的页面与图片逐文件保持不变。
- 原生 CSS/JavaScript 实现本地全文搜索、页内目录、语种切换、原生代码折叠、
  系统/浅色/深色主题、键盘跳至正文和响应式导航。无 CDN、下载字体或前端框架。

## 用户图标

用户于 2026-11-01 提供的原画作为首页品牌展示、README 与教程入口图片，导航使用
独立缩略图。完整图像仍为 1254 × 1254、RGBA；未裁剪、重绘、量化或缩小分辨率。

| 编码 | 字节 | 相比输入 |
| --- | ---: | ---: |
| 用户输入 PNG | 3,195,879 | 基线 |
| 无损优化 PNG | 3,048,867 | 减少 4.60% |
| 网页优先无损 WebP | 2,127,126 | 减少 33.44% |

两份完整图片均已与输入逐像素核对 RGBA，透明度保持。原始 RGBA SHA256 为
`85d36ae6f6e618f64e4c5d909314bd64072b0087bf093d24f2d41ac4dab133b3`。
[图标清单](../website/assets/icon.json)记录来源及输出 hash；64 × 64 缩略图明确标为
缩略资产，不把它称为全尺寸无损图。此前的原创六边形研究标识保留在设计来源页作
补充，不替代用户原画。来源和使用授权记录与源码的 MIT 条款分开表述。

## 实际验证

| 检查 | 实际证据与边界 |
| --- | --- |
| `Rscript tools/documentation-plots.R` | 六张显式尺寸 PNG 导出；当前几何/标签/轨道/局部视图/迁移组合执行，另五种组合单次 build；不是发布设备矩阵 |
| `Rscript tools/validate-documentation.R` | 仓库原始 FASTA/BLAST/GFF3 读取：4 序列、32 连接、266 feature；一个 8 × 5 in PNG 预览/诊断通过，panel 为 measured，不解释为全面视觉认证 |
| `Rscript tools/freeze-api.R check` | 130 个候选导出签名、默认、别名及上游对象身份保持；没有刷新 v0.16.4 快照 |
| `python3 tools/sync-tutorials.py --check` | README、迁移、vignette 源均同步；只写源码，不 knit |
| `python3 tools/build-website.py` | 34 个专题页面 + 186 个完整参考页面；98 条对应旧地址和 80 条无一一对应入口的历史提示页，共 398 HTML；生成到被忽略的 `.site/` |
| `python3 tools/check-documentation.py` | 双语配对、相同代码、内部链接/锚点、图片 hash、全部导出/Rd 分类、完整搜索索引与语种元数据通过；不声称已扫描全部外部网址 |
| 源码静态检查 | 86 个 R、93 个 Rd 无警告解析；JavaScript 语法、工作流 YAML 和 diff 空白检查通过 |
| 本地浏览器 | 1280 × 720 与 390 × 844；搜索、同页语种锚点、三态主题、Enter 跳转/折叠/菜单通过；窄屏正文无横向溢出，宽表自身滚动；预览 tab 无 error/warn 日志 |

环境为 R 4.6.0、ggplot2 4.0.3。机器可读记录见
[例子验证](../website/verification.json)和[浏览器记录](../website/browser-verification.json)。
图标首页截图随本次会话交付，不成为网站运行依赖。基础 R 手册链接使用官方
`R/refmans` 路径，第三方包使用 `CRAN/refmans`；已核对 grid/unit 与 utils/read.table
对应官方页面。

未运行 testthat、包 build/check、完整 vignette 渲染或 PNG/PDF/SVG 发布矩阵。
未运行远端 Actions，也不把本机浏览器检查写成多浏览器、WCAG 或平台认证。

## 网站发布与后续边界

本版完成的是当前授权的本地文档/网站源开发；原版本设计明确保留的网站冻结和
“迁移上线另行确认”边界继续有效，不因本地完成而擅自改动线上站点。

`.github/workflows/documentation.yaml` 的普通 CI 只做源码检查。手动工作流默认
`deploy=false`，构建可审阅的 preview artifact；明确网站切换授权后才使用部署输入，
并在远端启用 Pages Actions、确认环境保护和源提交。授权构建使用 `--publication`
标记尚未正式发布的候选文档，不把网站上线称为 R 包正式发布。

构建清单记录源 commit、dirty 状态和 SHA256；实际部署 URL/run/commit 必须在
授权切换时另行核对。新站上线并验证旧链接之后才移除 main 中冻结 `docs/`；本轮
不执行这一步。网站设计、源文件和本地验证已完成，不把尚未授权的线上动作混作
“已发布”。下一阶段为 v0.16.6，等待用户启动；其正式发布审计按另行授权执行。
