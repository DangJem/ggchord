# 网站审计后重构实况

本轮由用户于 2026-11-14 授权，基线为 0f5dd96，开发分支 v0.27.0。
仅本地提交，不合并、不推送、不部署或发布。仅验证电脑端，不测试移动端。
审计材料来源哈希、62 旧入口、122 参考条目、15 条问题记录（V01–V14 加 MIG-14）
及 582 个 URL 见 baseline.json；材料中的断言不自动等于本轮验证结果。

## 已提交阶段

- bdce3bc：审计基线账本；本地被忽略的设计文档同步当前授权与范围。
- e081ca7：安装、投影、GFF3、bar/col、位置类型、Fusion 范围及公开别名路由。
- f9c67e6：页面注册表、显式栏目首页、四模块和 learn、版本元数据及 API 回链。
- 62bbb05、a6cfe48：共同首课、四视图比较、Sankey/Radial 正式课程及独立准备修正。

## 已实现的课程与查询重构

Circular 首课直接消费已核对的 pUC19 长度和 ori 注释；识别属于后续独立工作流。
Dotplot 有独立核心课；Fusion 按首图、来源路由、物理布局、局部细节和进阶食谱拆分。
formats、tracks、output、foundation 等长页按任务归属。共同表格、几何/参数/标签
与统计主课保留原始解释及对照；ribbon 处理与统计合并。旧页/章节映射在
website/page-migrations.json，渲染与检索使用同一注册表；不将旧链接全部送往首页。
76 个原字面 R 块提取到 tools/site-course-examples.R，与原共享示例共同生成下载脚本。

图库保留 178 个原 figure_id、原图与来源，增加精选/模块/类型/任务/难度筛选及 hash
显露；搜索按任务意图、精确函数及参数选择入口，同页结果合并。辅助页与 API 类型分开。
README 与 vignette 源文件同步生成；没有编译 vignette。

## 已执行的轻量验证

- `python3 tools/check-documentation.py --source-only`：94 对入口、157 导出、121 Rd、
  双语代码及锚点、参数/API 生成同步、图像哈希及源链接检查通过。
- `Rscript --vanilla tools/check-documentation-lessons.R --write-outputs`：53 页在一个
  新 R 进程的独立环境执行；并对之后的 data/ribbons/multi-sequence-links 变更分别
  用 `--page=… --update-outputs` 复核。预览/保存调用只检查参数；PNG 另行实际绘制。
- `python3 tools/check-first-lessons.py`：四模块首课各自新 R 进程，完成公开 build/export。
- `python3 tools/render-site-courses.py`：22 张新增教学 PNG 按正常物理尺寸实际输出；
  逐图 SHA 和代码哈希见 course-rendering.json，不作为发布设备矩阵或整体几何认证。
- `python3 tools/site-audit/verify-existing-assets.py`：178 原 PNG 逐项 MD5 相同；
  R/api-contract.R 只改 roxygen，经 R parse 表达式相同；原共享代码仅将一个真实
  GFF3 读取显式改为 keep，已证明 24 条 CDS 输入逐项相同。源数据、运行时代码
  和原依赖/图像哈希逐项比较，证明等价。未刷新旧生成 stamp 或冒称本轮重绘。
- 1149 段 Rd 译文及保护标识检查通过；76 个颜色 token 对比度与 6 项静态规则通过。

安装指令未联网重新安装。可选解析器的缺失/NULL 不计导入成功；本机本轮 GenBank
块返回 NULL，原 Biopython 1.85 成功记录仍仅为历史证据。其他可选分支也不因页面
执行成功推断全部解码器均已运行。原反向 HSP 与 Seqinfo 警告保留。

## 已发现并修正的本轮问题

四视图等显示长度小节缺 foundation_data，圆图布局缺 firstcircle/firstfeatures 准备，
以及表格迁移的正则过度匹配，均由独立执行或源锚点检查发现并修正后复核。
首课验证工具最初让 R 顶层自动打印图形，出现 write failed；改为显式 eval 构造并
仅执行指定 build/export 后，四进程完整通过。前述失败没有计为通过。

## 本轮中断恢复后的提交与核验

- fe82652：课程拆分、旧章节迁移及搜索/图库。
- 95f01b8、8ff87bc：旧锚点、分批原生渲染、独立 Fusion 子图及 Markdown 分隔修正。
- 1eb0891：停止按参数名/默认值猜类型，引用同一 API 的 Rd 契约；高风险参数显式元数据；
  所有当前 GFF3 教学读取显式 keep；蛋白课将 dplyr 正文提醒改成可执行代码。
- 7e6a4c0：reference 迁移收尾顺序、单并发四页 worker、桌面顶栏和真实图像替代文字。

修改 GFF3 后，52 个当前可执行课页重新在新 R 进程中的独立环境运行；蛋白课和
aligned-sequences 又各自新 R 进程复核。前文 53 页是较早课程归属阶段的事实，
当前登记执行页为 52。表格/控制台 JSON 值逐项相同，仅键顺序改变。

最终运行 `python3 tools/build-website.py --render-jobs 1 --render-batch-size 4` 成功，
英语和汉语各 19 个 worker 全部退出 0；源输入快照贯穿原生渲染和 finalization。
随后用严格受限的 `refresh-derived.py` 仅重新生成搜索索引，修正空 formal 被误称
必需参数的摘要；HTML/示例没有再渲染，刷新范围及前后源哈希独立记录。

- 完整 `python3 tools/check-documentation.py` 通过：94 对入口、157 导出、121 Rd，
  双语代码/锚点、全部生成链接、搜索目标和资产。类型/单位未凭默认值推断；未逐项
  声称已独立执行所有 1671 条参数记录，精确契约仍以对应 API 的 Rd 为准。
- `python3 tools/site-audit/check-built.py` 通过：582 个审计地址、1146 个声明旧章节、
  14 个误退役公开别名、178 原图库 ID、22 新图代码/PNG 哈希和实际下载文件。
- `python3 tools/long-term-audit/check-native-archive.py` 通过：ZIP 中 76 个原输入、
  30 个由下载脚本字面路径展开的文件与当前源字节/快照一致；本轮结果另存为
  native-archive-check.json，未覆盖历史记录。
- `python3 tools/site-audit/check-http.py` 对本机预览通过：582 个 URL 返回 200，
  未知地址实际返回 HTTP 404，独立 404.html 提供双语找回入口。本地 Python 服务器
  的未知地址使用其默认 404 响应体；托管平台自定义错误页/重写尚未部署或实测。
- `Rscript --vanilla tools/site-audit/check-semantics.R` 通过：真实 rps16 分段 CDS 上，
  一个 feature 返回一行/两段，明确的同残基成对控制返回两行；windows、bar 拒绝
  非零基线、两个数值 position 入口直接验证通过。选定残基控制不冒称生物学比对。

构建失败原因和恢复策略见 build-recovery.json。两并发十二页批次曾在 real-data
退出 1，日志没有 R 错误条件；单页渲染通过，不据此断言 OOM。最终缩小 worker
后全量成功。原生 finalizer 的重复调用不能作为可靠恢复路径，最终使用干净构建。

## 仍待实际电脑端验收

原生 Chrome 曾在 1440px 发现顶栏版本重复/文字换行，相关 CSS 已修正。最终构建
完成后准备复核时，Mac 返回锁屏状态；已请用户手动解锁。浏览器扩展此前不可用，
未绕过系统锁屏。1280/1440/1920px 最终版布局、两语切换、搜索排序/快捷键/焦点、
图库筛选/旧 hash 显露、图片弹层尚未通过本轮最终浏览器验收，见 desktop-check.json。
截图或静态检查不能替代这项待办；不将本轮描述为全部验收完成。
不运行 testthat、R CMD build/check、vignette 编译或 PNG/PDF/SVG 发布矩阵。
