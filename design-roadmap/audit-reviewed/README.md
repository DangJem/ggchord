# 逐参数人工复核输入

本目录保存以具体源码为依据的人工语义记录，不是静态脚本根据参数名称推测的合同。
当前首组覆盖 ggchord、clean_ggchord_data、prepare_ggchord_data 和
validate_ggchord_data。来源提交、复核证据、实际验证范围均保存在各 CSV 行中。

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
