# 网站说明书重构与 ribbon 样式修复实况

## 当前授权与边界

开发分支为 `v0.27.0`，基线为 `0f5dd96`，仍使用原 OneDrive checkout。
仅本地分阶段提交，不合并 `main`、不推送、不部署、不发布。冻结 `docs/` 未重建。
用户于 2026-11-14 明确取消以后所有电脑端交互验收；移动端也不测试。
本记录仅报告源文件、公开 R build、正常尺寸 PNG、静态 HTML、链接与下载核验，
不声称浏览器布局、菜单、搜索键盘操作或视口验收通过。

三份审计材料仅作为需求来源：原文断言与本轮实际验证分开。材料哈希、旧入口、
参考条目、问题账本和 582 个 URL 见 `baseline.json`。仓库长期规则及本地设计文档
保持适用；最新要求已记入被本地忽略的 `design-roadmap/versions/v0.27.0.md`。

## 已提交阶段

- `bdce3bc`：建立材料、URL、接口与审计基线。
- `e081ca7`：投影、GFF3、bar/col、位置类型、Fusion 范围和公开别名路由。
- `f9c67e6`、`62bbb05`、`a6cfe48`、`fe82652`：初始模块注册表、可执行说明、拆页、查询与图库。
- `95f01b8`、`8ff87bc`、`1eb0891`、`7e6a4c0`、`393ed6d`：旧锚点、原生分批渲染、参数契约及前阶段核验记录。
- `0fee412`：跨行 ribbon 自动样式限定在 Collinear；Sankey 角度元数据的 NA 安全处理。
- `e9b88aa`：按最新复核意见重构四模块说明书、导航、首页和真实旗舰图。
- `7cba801`：完成说明书措辞补查，修正旧模块编号和汉语兼容表多列问题。
- `4b059e3`：历史汉语入口改连通用安装说明，去掉生成模板中的学习总览文案。
- `c6067f4`：开始使用菜单的稳定双语标识、按注册表复制图像及严格受限的失败恢复/静态输出修正。

前阶段的 94 对页面、53 页准备、178 张旧图字节未变和等待桌面解锁，只是当时事实；
以下状态以最新要求和本轮执行结果为准。

## 最新十三项复核意见的落地

| 用户意见 | 当前实现与事实来源 |
| --- | --- |
| 四大模块 | `pages.json`、`navigation.json`、基本概念和首页统一为多序列关系、单序列圆图、HSP 点图、多视图组合；四种多序列 Coord 与 Fusion 组合层职责明确。 |
| 紧凑首页与旗舰图 | 首页移走长篇说明，放大 logo，四模块能力卡；七张首页图均用 `#FFFCF5`。原 Chord/Collinear 的数据选择、图层与构图保留，Radial/Sankey 从它们派生。Chord 仅因本次样式 bug 实际重绘。 |
| 工具说明书口吻 | 双语正文、导航、图注、搜索与生成器改用功能/输入/参数/结果措辞；内部 recipe/lesson 元数据键和历史文件名保留。 |
| 顶栏 | 新增首页入口、统一控件高度、去掉重复版本文字、加宽菜单。仅源规则检查，未做浏览器高度/交互验收。 |
| 开始使用 | 删除 learn 两语正文与索引入口，旧 URL 只跳转；安装改为通用未固定版本的 GitHub 方法；基本概念覆盖四模块、六 Coord 和组合层。 |
| 表格 | 表头、间距、斑马纹、分隔线及数字对齐统一；宽表局部滚动。 |
| 模块标题与编号 | 四模块采用“模块一 · …”风格；多序列先共享内容，再分 Chord/Collinear/Radial/Sankey；左右导航和正文使用同一章节编号源。 |
| 13 个质粒 | 采用保留完整来源特征、独立参考引物及声明酶 profile 全部位点的 13 张已有正常尺寸 PNG；停止三切点抽样。`manual-figures.json` 与 `reference-catalog.json` 记录来源和计数。 |
| 圆图案例分工 | pUC19 专讲完整图层构成；pETDuet-1 独立说明完整输入下的密集标签与尺寸，去掉不存在的上文引用。 |
| HSP 点图 | 补四序列、32 条真实 HSP 双轴矩阵，以及方向、Identity/长度映射、拼接坐标、局部查看和输入来源。 |
| 功能命名 | 简化导航/页名，局部细节改为“局部放大与字符细节”等直白名称。 |
| 食谱措辞 | 改为“嵌套组合与跨原点放大”；根因是内部 recipe 术语直译进读者文案，现统一修正展示层和生成器。 |
| 两侧字号 | 两侧导航提高至 .96rem，并提高分组/当前章辨识度。 |

首页按普通桌面高度尽量一屏，矮窗口保留自然滚动，不通过隐藏溢出裁掉内容。
没有浏览器验收，因此不保证某一分辨率完全无滚动条。

## 最新 ribbon 缺陷

Chord 原布局会按非相邻关系设虚线；Radial/Sankey 又复用线性预布局，连同跨行
自动线宽一起继承了 Collinear 约定。现在默认跨行虚线和加粗只用于 Collinear，
Chord/Radial/Sankey 默认实线与普通线宽；显式固定及映射线型继续优先。
定向脚本通过四种公开 Coord build，比较修复前后 x/y/group，几何逐项哈希一致。
Collinear 仍为 solid/dashed/dashed，其他三种均 solid/solid/solid，普通线宽 .12；
固定 dotted、映射 dotdash/longdash/dotted 均保留。见 `ribbon-styles-before/after.json`。
复用现有最小回归案例扩充四 Coord 断言，未运行 testthat 或包级检查。

Sankey 实际 PNG 绘制还暴露混合文字/引线元数据中的 NA 角度索引错误，现用 `which()`
限定有效角度，保留 NA。新 Sankey 旗舰只选每条序列一个来源标签做固定放置，
所有特征几何仍保留；这不代表全部 Sankey 引线布局已做质量认证。

## 已执行的轻量验证

- `python3 tools/check-documentation.py --source-only`：93 对页面、157 导出、121 Rd、
  双语可执行代码、源锚点、API/参数生成同步、1149 段保护译文与图像依赖检查通过。
- `Rscript --vanilla tools/check-documentation-lessons.R --write-outputs`：当前 52 页在
  新 R 进程的独立环境准备成功；包命名空间可共享。最终 Sankey 标签调整后又以
  `--page=multi-sequence-views --update-outputs` 单独新进程复核。
  预览/保存参数只检查参数，PNG 另行真实绘制。
- `Rscript --vanilla tools/site-audit/check-ribbon-coordinate-styles.R`：四种 Coord 的
  默认/固定/映射线型和线宽检查；与修复前几何比较通过。
- `Rscript --vanilla tools/site-audit/refresh-ribbon-figures.R`：原 178 图中 61 张
  Chord/Radial PNG 按登记正常尺寸实际重绘，117 张按精确运行时差异及源/像素哈希保留。
  `verify-existing-assets.py` 对所有依赖逐项核对，未将全部旧图冒称本轮重绘。
- 22 张前阶段新增 PNG 保留各自真实渲染记录；本次受影响的 coursechord、courseradial、
  courseradialangle 三图重绘。未受影响图不刷新虚假的生成 stamp。
- `Rscript --vanilla tools/render-site-flagships.R flagradial flagsankey flagdotplot flagfusion`：
  四张新旗舰 PNG 实绘，每图保存绘制前后的源哈希与 PNG SHA。未重做原两张多序列旗舰。
- `python3 tools/site-audit/check-manual-assets.py`：20 张说明书图、13 完整质粒图和七张
  首页图的统一背景检查通过。已有完整质粒图复用此前测量/布局证据，不冒称重新绘制。
- `python3 tools/check-site-design.py --write`：76 个受控颜色对与六条静态规则通过；
  不作为 WCAG 或浏览器认证。渲染模板已通过 R parse，`git diff --check` 通过。

四模块最小示例曾各在新 R 进程执行 public build/export，记录在 `first-lessons.json`，
该记录只对应当时脚本/引擎，不作为本次引擎的重新验收。安装指令未联网安装。
可选 GenBank 的本机 NULL 属跳过；旧 Biopython 成功证据只代表历史记录。
反向 HSP 和 Seqinfo 原警告保留，不能因页面准备成功推断所有可选解码器均已运行。

## 构建中补查与修正

本次有两轮构建在源复核时主动中止，均未计为成功：先发现少量“本课 / Exercise”残留、
共线性旧模块编号以及兼容表的一行额外列。修正纯文案和表格后重新完整构建；
随后发现历史汉语入口模板还指向旧总览，改为通用安装说明。源检查增加正文措辞、
表格列数和历史入口检查，重新通过。未修改原始输入或本次绘图引擎。

## 最终构建状态

双语原生构建与收尾已完成，`native-build.json` 记录实际过程：每种语言各 18 个初始
worker，汉语全部成功；英语前 17 个成功，最后一个在 fusion-details 异常退出且无
明确 R 错误消息。保留此失败，没有推断为内存问题或记为通过。随后以同一源快照
将该批四页分别在新 R 进程重跑，全部返回 0，才执行最终处理。共 40 次实际 worker
调用（35 次初始成功、1 次初始失败、4 次隔离重跑成功）。

完整静态校验随后发现开始使用菜单的原生双语 ID 不同，以及输出目录缺少一张
showcasecollinear PNG。修正生成器，并通过 `finish-manual-static.py` 证明只有两个
精确源差异：开始使用菜单 ID/ARIA 规范化、按注册表逐项复制图像。现有输出只更改
388 个 HTML 的两处菜单属性并补一张原图，正文、可执行代码和输入包未重绘/改动。
源快照保留此次受限派生修正链，未把它冒称全站重新渲染。

之后全部重新通过：`check-documentation.py` 的 93 对页面/121 Rd/157 导出、双语
代码与生成链接/锚点/搜索；`check-built.py` 的 582 地址、1146 旧章节、14 公开别名、
178 原图 ID、22 张新增例图与下载哈希；`check-manual-assets.py --built` 的 20 张
说明书图片；ZIP 76 输入与展开代码 30 路径均与当前源快照一致。

旧预览地址/系统代理最初返回 503，不计为通过。另启仅绑定 127.0.0.1 的本机服务器，
绕过代理后，582 个真实 HTTP 请求均为 200，未知地址为 404，双语恢复页可取。
没有访问或部署远程网站。尝试用未安装的 bs4 读取 HTML 后立即改用已有 xml2，
未安装依赖；静态读取确认模块标题及左右编号来自同一章号，未代替浏览器验收。

最终只更新核验报告并本地提交。未合并 main、未推送、未发布，未重建冻结 docs/；
未进行浏览器/移动端验收、包级检查、vignette 重建或发布设备矩阵。
