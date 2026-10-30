# 逐参数人工复核输入

本目录保存以具体源码为依据的人工语义记录，不是静态脚本根据参数名称推测的合同。
当前包括构造/数据四入口、三个坐标入口、四种自有 Position、seq/gene、两个
feature 入口、link/ribbon/bundle/density、主题、scale/guide和输出诊断，共 1,050 行、64 个公开入口；这些入口的
当前显式形参均已有源码语义记录。其余公开参数视图仍待复核，不能把本批入口覆盖
当作全部 API 审计完成。来源提交、复核证据和实际验证范围均保存在各 CSV 行中。

本轮新增的三表分别为 `seq-gene-geoms.csv`（27 行）、`feature-geoms.csv`（33 行，
data 两行沿用既有表）与 `links-stats.csv`（54 行）。F33 相关行在修复后重新阅读
共用输入列附加路径并单独更新来源；其他行保留原复核提交。构造器 ribbon_data 的
消费者说明亦已更正：point link 需要自身输入，不自动回退到构造器连接表。

持续开发批次又补入 `themes.csv` 的705行（141个共同形参的五份预设视图）、
`scales-guides.csv` 的140行和 `output-diagnostics.csv` 的22行。共享实现采用相同
契约记录，不将重复签名行理解为独立功能。F35/F36位置标度和F37 alpha默认palette
相关行分别在修复后重审并更新来源；其余行保留原复核SHA。

`SOURCE_REVIEWED` 表示已经阅读相关源码路径，不能解释为所有组合均运行通过。
`test_case` 标记 NOT_RUN 的场景仍待验证；KEEP 只表示本轮建议保留对应现状，不是
稳定版兼容承诺。文档与实际行为的差异单独写明，未获选择的冲突不自动覆盖。

运行 `Rscript tools/audit-ggchord-api.R . design-roadmap/audit-generated` 时，本目录
CSV 按 function/parameter 唯一键合并到公开参数行；其余参数仍保持 REVIEW_REQUIRED。
签名、参数顺序、默认表达式和本次静态源码版本不允许由人工表覆盖。原始静态引用
证据与人工说明分别保留；缺少必需复核信息、重复键、错误列或无对应公开参数会使生成失败。

CSV 用 UTF-8 保存；R 读取时使用 `na.strings = character()`，避免把代码中的
字面表达式 `NA` 当成表格缺失值。改变源码后必须评估对应人工记录是否过期，不可
只更新 source_revision 字符串而不复核实际内容。
