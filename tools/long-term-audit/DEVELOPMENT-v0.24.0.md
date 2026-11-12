# v0.24.0 — Guide 联合布局与可选尺寸策略

状态：约定范围本地开发与定向验证完成，尚未整合、推送或正式发布。
2026-11-12按用户明确批准，将v0.23.0（4a9d755，含4项审计修复）快进整合本地main，
再从main创建v0.24.0。main保持4a9d755，origin/main仍为69556d8；本轮没有远端操作。
实现提交4b3db69；DESCRIPTION为0.24.0.9000，最近正式版仍为0.16.6。
本版增加可选公开能力，使用minor开发版本；修复与文档整理没有另建patch/minor版本。

## 实现与决策

依据v10 §6.3、§6.3.1、A.7，复用既有渲染上下文、原生Guide和gtable，不新增数据模型。

- `adapt_ggchord()`附加图对象策略，默认不启用。fixed保持默认主题/图例部件的参考
  物理尺寸；scale有界调整默认主题文字、线宽、间距及图例部件；adaptive进一步联合
  选择原生Guide行列、方向与侧边，复用已有标签布局并重测完整图形。
- 参考尺寸默认8×6英寸；文字倍率为目标/参考面积比的四次方根，限制在0.85–1.35；
  线宽和默认间距倍率为文字倍率平方根。字号下限默认2.5mm，已有更小默认保持原值。
  用户显式geom样式、映射、主题字段（包括继承及列表添加）、Guide行列/方向/位置/
  键尺寸、顺序、scale和show.legend优先；不改原始数据、统计或源区间。
- 先绘制原生Guide候选并测量毫米尺寸，以共享侧边预算和保留主面板空间为约束，
  优先保留现布局。有界确定性搜索保留有限联合候选，再构建实际gtable核对边界和
  面板可测性。搜索宽度/次数是内部策略，不是稳定API或全局最优证明。
- 连续及分段色条不拆分，保留原生刻度、顺序和数值比例；不通过删类别或无界缩字
  容纳图例。未知自定义Guide只测量。无法满足预算时返回insufficient_space并警告，
  报告约束与实际溢出。feasible只针对图例画布/预留预算，不认证全部标签。
- 普通print、ggsave、view、diagnose共享构建/渲染路径。Fusion策略在各子面板实际
  分配尺寸内执行，再计算跨面板端点；不自动重排外层面板或其独立关系图例。
- `diagnose_ggchord()`增加adaptive_layout与guide_layout，记录有效主题样式、倍率、
  显式锁定、最终行列/侧边、候选尺寸和最终图例框毫米边界。纯布局导出只记录策略及
  target_device_required，不把估算当作目标设备结果。分面在策略添加或构建时明确拒绝。

没有删除公开接口、增加运行时依赖或改变未启用策略的默认设计。

## 定向证据

`Rscript tools/long-term-audit/check-adaptive.R`通过，结果为
[v0.24.0-adaptive-results.json](v0.24.0-adaptive-results.json)。

- 明确标注的合成契约数据：10组Seq、10组Feature类别、Identity六刻度，6×4英寸。
  最终图例框均在画布内，类别数量及原始序列/HSP表保留。同尺寸再次构建选择相同，
  原图与未添加策略的父图不变。
- 真实seq_data_example/ribbon_data_example：fixed的倍率为1，scale在10×8英寸
  增大默认样式但仍受界限约束；纯布局来源导出保持一致。
- 显式4mm geom标签、主题18pt序列标签及12pt图例文字保持对应优先级；后续仅改
  颜色不丢弃字号锁定。父级title字段、列表添加主题及显式底部/水平/两列Guide均保留。
- Chord之外，Collinear、Sankey、Radial、Dotplot、Circular各一个代表性共享构建成功。
  无图例的图返回no_guides，不把空表当作图例测量证据。
- pUC19源图加局部Collinear的Fusion案例：两个策略各使用分配尺寸，真实101–1900
  选区及两条装饰引线保留；未改变源区间或统计口径。
- 显式锁定100mm宽图例键的小画布返回insufficient_space且报告溢出。

两张代表PNG已人工查看：[合成多Guide](v0.24.0-adaptive-guides.png)、
[真实四序列与scale策略](v0.24.0-adaptive-real.png)。仅用于所记录场景，非设备矩阵。
另以相同真实输入运行一次标准ggsave到临时PNG，验证实际保存路径；后加facet的拒绝检查通过。
新增一个最小公开API案例留待发布授权后运行testthat，本轮未运行测试套件。

## 文档同步

NEWS、roxygen/Rd、英语/汉语样式课、支持矩阵、安装版本、API课程、参数清单与参考已同步。
同时补齐此前审计修复涉及的Track和layout contract汉语参考段落，未改写历史事实。

- `Rscript tools/check-documentation-lessons.R --page=styling`通过；独立R进程执行
  真实数据策略及诊断，证据[v0.24.0-styling-lesson.json](v0.24.0-styling-lesson.json)。
- `python3 tools/reference-translations.py --check`：1045译段及受保护代码/链接一致。
- `generate-parameter-atlas.py --check`、`generate-api-curriculum.py --check`、
  `sync-tutorials.py --check`通过；147导出入口清单含新策略。
- 62对页面原始R块和锚点一致，git diff --check通过。Rmd仅同步源码，没有编译。

## 剩余范围与发布边界

本版完成已约定的原生Guide联合布局和可选策略，不等于v10整体完成。
新全局标签求解器、任意facet/grob、自定义Guide重排、所有几何常量缩放、字体嵌入认证、
Fusion外层联合Guide仍不支持。M2/M3、F2/F3、新格式适配器、1.0前API治理继续保留。

没有运行testthat、R CMD build/check、vignette重建或PNG/PDF/SVG发布矩阵，没有重建网站。
此前5个旧lesson engine与178旧图指纹过期待重现，旧CI记录仍未由新推送验证；不能将这些
旧资产声称为本版重新生成。源码、定向证据、完整发布验收和正式发布是不同状态。
