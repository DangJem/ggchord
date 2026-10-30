# ggchord API 静态审计基线

此目录由 `tools/audit-ggchord-api.R` 生成。完整覆盖状态见 manifest.csv；源码复核、局部运行验证与发布验收分别记录。

起始提交：`69b907a6ddffea60a73ed1d4b1392c8cc5f7da22`；生成时 HEAD：`4724092f5fdc03e7bb0367537faf6fdb52f8c794`。
输入指纹（MD5，仅作可重复性标识）：`4cc5474f5dd8ea45b1347b7d5857ef06`。
源码、已保存 Rd 与人工复核输入的逐文件散列见 source_hashes.csv；工作树状态单独记录，不能把带改动源码冒充原始提交。
解析前固定输入文件列表、散列与 HEAD，写表前再次核对；并发改动会使生成失败并要求重跑，不写入不一致的新产物。

## 重建

```sh
Rscript tools/audit-ggchord-api.R . design-roadmap/audit-generated 69b907a6ddffea60a73ed1d4b1392c8cc5f7da22
```

读取 CSV 时使用 `read.csv(path, na.strings = character(), check.names = FALSE)`，以免把代码默认表达式 NA 误读为 CSV 缺失值。

## 覆盖

- 83 个 R 文件全部成功解析；897 个 FUNCTION token 与访问数逐文件相等。
- 130 个导出名称，22 个重导出，9 个 S3 注册；导出不等于独立函数。
- 897 个函数语法定义（含匿名、嵌套及 ggproto 方法）；3302 行参数含导出别名/预设的签名视图。
- api_objects.csv 补充顶层绑定、别名、ggproto 构造和工厂调用；api_dynamic_registration.csv 标记 quote/eval、注册、函数体/环境修改等。
- 颜色拼写别名通过源绑定追踪；theme 预设的签名来自工厂复制源码，运行时形参另由108入口签名脚本核对。
- api_arguments.csv 保存全部计划字段、默认表达式及 Rd 参数说明；未知语义显式 REVIEW_REQUIRED，文档说明不等于行为证据。
- 默认表达式占位：<required> 表示无语法默认值，真实是否必填仍须核查 missing() 等分支；<ellipsis> 表示可为空的 ...，不标为必填。
- api_aesthetics.csv 记录静态美学声明/映射或待查表达式；api_return_contracts.csv 保留文档返回说明及语法尾表达式，不推定稳定字段。
- api_compatibility_matrix.csv 合并人工三坐标源码合同；runtime_validation 单独说明有限例子或 NOT_RUN，没有自动执行绘图。

## 人工语义复核

- 本次纳入 37 份 audit-reviewed/*.csv，3302 行参数、108 个公开入口标记为已做来源或运行复核；其余语义仍待复核。
- 参数表按 function+parameter 匹配公开参数；有 function_id 时必须精确匹配唯一内部函数及其参数。重复键、未知字段、签名字段或未定位键使生成失败。
- 默认表达式、形参顺序、签名证据、当前 source_revision 不接受覆盖；原始 static_confidence/static_actual_usage 与语义复核分列。
- review_source_revision、review_input_path、review_status、review_evidence 保留人工复核来源；复核提交与当前提交不同须核查相关源文件差异，不能自动把旧证据升级为当前行为。
- audit-reviewed-details/ 下的 export-decisions 与 coordinate-compatibility 按不同键合并返回、去留和兼容记录；要求完整且无多余键，保留独立 review SHA。
- internal-functions 按私有 function_id 合并角色、返回、状态与动态依赖证据；未覆盖的函数仍保留 REVIEW_REQUIRED，不将整文件解析冒充语义复核。
- layout-stages 记录七个共享帧阶段的读写、顺序、状态和风险；精确挂接对应的 top-level quote 绑定及 eval 调用，其他动态记录仍待查。
- aesthetic-contracts 按文件/owner/声明种类/美学或占位/表达式精确对应当前声明；input_role、draw aesthetic与动态default集合分开，源码复核不等于所有运行集合已穷举。
- 对象表区分已复核函数绑定、公开别名/工厂和延迟阶段内的局部赋值语法；依赖/动态表只挂接caller合同，保持SYNTAX_ONLY，不将其变成删除或运行可达性结论。
- object-contracts逐项记录非函数原型/常量/环境及local临时绑定；动态记录分别挂接函数、延迟阶段或原型构造作用域。状态表示源码拥有者已审，不表示动态结果已穷举。
- 没有相应人工覆盖文件时照常生成，语义保持 REVIEW_REQUIRED；两个人工目录均纳入输入快照校验。

## 局限与后续工作

- 全部文件被解析不等于全部动态对象被穷举。运行生成的函数、动态注册、字符串回调、NSE、ggproto 继承及语言对象注入必须人工/定向运行复核。
- 函数行号为 FUNCTION token 起点；顶层绑定和调用证据行号为所属表达式/函数的起点，location_precision 或证据类型说明精度。
- 普通调用语法、被引用调用语法和符号出现次数均不是真实运行调用图。零引用既不证明无用，也不构成删除授权。
- 声明型 default/required aesthetic 可提取；经运行修改/继承得来的最终集合仍待复核，重复声明不机械合并。
- 不加载包，不解析或执行原始业务数据，不运行测试、包检查、vignette 或设备矩阵。语义、验证路径、NULL/NA、单位、迁移及去留决定仍需逐项证据。
- 部分公开参数通过外部 ggplot2/grid 的 ... 透传，其完整外部签名不由本仓库静态函数列表重建。
