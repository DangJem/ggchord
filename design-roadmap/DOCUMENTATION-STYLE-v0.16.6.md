# v0.16.6 当前文档示例风格审阅

2026-11-02，在 `codex/v0.16.6` 实施；比较基线为 `bfd07d6`。用户要求当前网站、
说明文档的代码优先使用 tidyverse 风格，并审阅代码是否清晰、简洁、易于维护。
长期偏好已写入 AGENTS.md，本版范围与检查边界已写入版本文档。

## 实际修改

共享源 `tools/documentation-examples.R` 使用 dplyr 的 select、filter、transmute、
mutate、rename、group_by、slice_max、distinct、join 和 bind_rows；采用 R >=4.1
支持的原生 `|>`，不使用更新 R 才支持的管道占位符。保留有意义的输入、整理结果、
绘图对象及数据来源说明，不把单个管道拉成整个工作流。DNA 反向互补和 GC 比例
保留短小的向量函数，文件读写和绘图对象列表迭代保留直接的 base 写法。

同源 CDS 由一次提取加明确筛选得到，显式恢复输入序列的先后顺序；按最长区间选取
保留单一并列结果。DNA 仍按真实链方向计算，不用同名基因替代序列证据。定量轨道
的移动平均显式使用 `stats::filter()`，避免与 dplyr 的同名函数混淆。

网站两个语种的首图、安装和独立准备源已同步；README、迁移指南和双语 Rmd
从共享源重新生成，未编织 vignette。README 保留顶部 320px 图标和三种视图的
独立代码/图片，没有恢复已移除的 reallinear 展示。

36 个 R 文件中的帮助示例已整理，49 个 Rd 的修改仅涉及例子；汉语参考的 865 个
正文块和内联 token 保持原值，仅在核对正文后更新绑定摘要。57 个帮助示例块可解析。
删除了定量轨道、同源聚类等专门帮助中无关的完整绘图准备；focus、BED、CDS
投影及 primer 帮助也只保留其使用的输入与步骤。

当前独立教学脚本中的表格使用 tibble，步骤和图层分行。可选 Bioconductor 的 S4
构造保持其原有接口。原始 FASTA/GFF3/BLAST、随包数据、历史验收和冻结 docs
未改写。dplyr 已在 Suggests，无新增包运行时依赖；独立 GFF3 脚本明确加载
其使用的 dplyr/readr/stringr，不再加载整个 tidyverse。

## 审阅中发现并处理的问题

| 编号 | 问题与处理 | 验证 |
| --- | --- | --- |
| D01 | 真实数据共线性示例重复加入相同 doc_axis 标度；保留最终有效的一处 | 重新输出的对应 PNG 与原图字节相同 |
| D02 | 独立 GFF3 转换器的文件正则未锚定扩展名，product 提取要求末尾分号；改为明确文件匹配、清晰分步和接受最后一个 product 属性 | 当前 266 条 CDS 的列值及顺序相同；有/无末尾分号的最小输入都保留 AmpR；未写出或覆盖原输入/既有表 |
| D03 | 可选列通过 `$` 访问，tibble 对缺失列发出多余警告 | validation 的 length/pident、focus 来源列、轨道截断标志改用精确 `[[ ]]` 读取；保留 NULL 回退和既有必需列检查 |

D03 涉及 `R/validate.R`、`R/focus.R`、`R/focus-track.R` 的五处读取，没有改变公开
签名、数据结果或引入依赖。最小 tibble/普通表格对照保留缺失值、截断侧别和重复
focus 的来源行；不以压制所有警告处理问题。真实反向连接的既有校验提示仍保留。

## 实际验证及边界

机器可读摘要见 [documentation-style-v0.16.6.json](audit-generated/documentation-style-v0.16.6.json)。

- 以基线共享源和当前共享源分别准备数据，对照 57 组输入/整理/投影结果的业务列、
  数值和行序，全部一致。比较将表格统一为普通列列表；不把 tibble 类、自动行名、
  自带属性或转换报告中的自动 source_names 当成业务标识。显式来源行字段仍比较。
- 帮助改写先通过 R 表达式解析对照确认实现未变；之后另作上述五处 D03 修复。
  当前运行实现差异限于这三文件的五处可选读取。
- `Rscript tools/documentation-plots.R` 实际运行全部共享示例和 21 张既有 PNG；
  图片尺寸、全部 PNG 字节及输出计数与基线一致。manifest 的 code_md5 对应本次
  实际执行的共享源，没有仅改摘要复用旧图。不是发布 PNG/PDF/SVG 设备矩阵。
- `Rscript tools/check-documentation-lessons.R` 的 20 页准备通过；隔离是同一进程
  中的独立环境，不代表干净安装或完整跨平台运行。
- 五个当前教学脚本完成对象准备；57 个帮助块解析通过，未逐个执行全部帮助。
  GFF3 转换器仅在内存对照实际数据和最小属性例，不覆盖 gene_track.tsv。
- `Rscript tools/freeze-api.R check` 确认 130 个导出保持候选签名、别名和上游身份。
- `python3 tools/check-documentation.py --source-only` 通过：27 对语种、同源代码、
  图片摘要、130 exports/93 Rd，以及 865 个汉语正文翻译块一致。
- Biopython、GenomicRanges 在当前环境不可用，对应可选示例明确为 SKIP。

没有运行 testthat、包 build/check、完整 vignette、发布设备矩阵或网站部署。
没有修改 tests；此前测试组织调研的五类归并建议仍为待确认方案。
main 保持 `d697345`，本轮不整合、不 push、不发布。

风格参考：[tidyverse 管道规范](https://style.tidyverse.org/pipes.html)、
[语法规范](https://style.tidyverse.org/syntax.html)与
[dplyr filter 的缺失值/行序语义](https://dplyr.tidyverse.org/reference/filter.html)。
采用这些习惯时以示例语义和 R 下限为约束，不机械替换所有 base 表达式。
