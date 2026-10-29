# ggchord API 静态审计基线

此目录由 `tools/audit-ggchord-api.R` 生成。完整记录见 manifest.csv；不代表逐项语义审计已完成。

起始提交：`69b907a6ddffea60a73ed1d4b1392c8cc5f7da22`；生成时 HEAD：`2c90371ce2a27ee4083c07a396150e78c46040b0`。
输入指纹（MD5，仅作可重复性标识）：`3e2d3affc91c2bc96df7854606d59d2c`。
源码、已保存 Rd 与人工复核输入的逐文件散列见 source_hashes.csv；工作树状态单独记录，不能把带改动源码冒充原始提交。
解析前固定输入文件列表、散列与 HEAD，写表前再次核对；并发改动会使生成失败并要求重跑，不写入不一致的新产物。

## 重建

```sh
Rscript tools/audit-ggchord-api.R . design-roadmap/audit-generated 69b907a6ddffea60a73ed1d4b1392c8cc5f7da22
```

读取 CSV 时使用 `read.csv(path, na.strings = character(), check.names = FALSE)`，以免把代码默认表达式 NA 误读为 CSV 缺失值。

## 覆盖

- 83 个 R 文件全部成功解析；896 个 FUNCTION token 与访问数逐文件相等。
- 130 个导出名称，22 个重导出，9 个 S3 注册；导出不等于独立函数。
- 896 个函数语法定义（含匿名、嵌套及 ggproto 方法）；3301 行参数含导出别名/预设的签名视图。
- api_objects.csv 补充顶层绑定、别名、ggproto 构造和工厂调用；api_dynamic_registration.csv 标记 quote/eval、注册、函数体/环境修改等。
- 颜色拼写别名通过源绑定追踪；theme 预设的签名来自已核对工厂复制源码，仍标为需要运行时复核。
- api_arguments.csv 保存全部计划字段、默认表达式及 Rd 参数说明；未知语义显式 REVIEW_REQUIRED，文档说明不等于行为证据。
- 默认表达式占位：<required> 表示无语法默认值，真实是否必填仍须核查 missing() 等分支；<ellipsis> 表示可为空的 ...，不标为必填。
- api_aesthetics.csv 记录静态美学声明/映射或待查表达式；api_return_contracts.csv 保留文档返回说明及语法尾表达式，不推定稳定字段。
- api_compatibility_matrix.csv 仅建立三坐标与 geom/stat 的待验证矩阵，没有执行绘图。

## 人工语义复核

- 本次纳入 2 份 audit-reviewed/*.csv，29 行参数、4 个公开入口标记为已做来源或运行复核；其余语义仍待复核。
- 覆盖文件按 function+parameter 唯一键，只能匹配唯一 exported 参数；重复键、未知字段、签名字段或未定位键会使生成失败。
- 默认表达式、形参顺序、签名证据、当前 source_revision 不接受覆盖；原始 static_confidence/static_actual_usage 与语义复核分列。
- review_source_revision、review_input_path、review_status、review_evidence 保留人工复核来源；复核提交与当前提交不同须核查相关源文件差异，不能自动把旧证据升级为当前行为。
- 没有人工覆盖文件时照常生成，语义保持 REVIEW_REQUIRED；人工文件也纳入前后快照校验。

## 局限与后续工作

- 全部文件被解析不等于全部动态对象被穷举。运行生成的函数、动态注册、字符串回调、NSE、ggproto 继承及语言对象注入必须人工/定向运行复核。
- 函数行号为 FUNCTION token 起点；顶层绑定和调用证据行号为所属表达式/函数的起点，location_precision 或证据类型说明精度。
- 普通调用语法、被引用调用语法和符号出现次数均不是真实运行调用图。零引用既不证明无用，也不构成删除授权。
- 声明型 default/required aesthetic 可提取；经运行修改/继承得来的最终集合仍待复核，重复声明不机械合并。
- 不加载包，不解析或执行原始业务数据，不运行测试、包检查、vignette 或设备矩阵。语义、验证路径、NULL/NA、单位、迁移及去留决定仍需逐项证据。
- 部分公开参数通过外部 ggplot2/grid 的 ... 透传，其完整外部签名不由本仓库静态函数列表重建。
