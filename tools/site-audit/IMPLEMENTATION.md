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

## 本阶段已实现，最终 HTML 与电脑端复核待完成

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
- `python3 tools/render-site-courses.py`：20 张新增教学 PNG 按正常物理尺寸实际输出；
  逐图 SHA 和代码哈希见 course-rendering.json，不作为发布设备矩阵或整体几何认证。
- `python3 tools/site-audit/verify-existing-assets.py`：178 原 PNG 逐项 MD5 相同；
  R/api-contract.R 只改 roxygen，经 R parse 表达式相同、原共享代码/数据/运行时代码
  及原依赖哈希逐项比较，证明图像依赖等价。未刷新旧生成 stamp 或冒称本轮重绘。
- 1149 段 Rd 译文及保护标识检查通过；76 个颜色 token 对比度与 6 项静态规则通过。

安装指令未联网重新安装。可选解析器的缺失/NULL 不计导入成功；本机本轮 GenBank
块返回 NULL，原 Biopython 1.85 成功记录仍仅为历史证据。其他可选分支也不因页面
执行成功推断全部解码器均已运行。原反向 HSP 与 Seqinfo 警告保留。

## 已发现并修正的本轮问题

四视图等显示长度小节缺 foundation_data，圆图布局缺 firstcircle/firstfeatures 准备，
以及表格迁移的正则过度匹配，均由独立执行或源锚点检查发现并修正后复核。
首课验证工具最初让 R 顶层自动打印图形，出现 write failed；改为显式 eval 构造并
仅执行指定 build/export 后，四进程完整通过。前述失败没有计为通过。

最终 HTML、旧 URL/章节、搜索、下载包、实际桌面宽度和键盘操作尚待下阶段实核。
不运行 testthat、R CMD build/check、vignette 编译或 PNG/PDF/SVG 发布矩阵。
