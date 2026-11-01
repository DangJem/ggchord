## 建立显式编码映射 {#map}
这个合成例子的一个密码子跨越两段编码区间：

{{code:projection}}

`genomic$segments` 保留实际编码片段，不使用包含内含子的外包区间。
`project_ggchord_links()` 同样投影显式蛋白连接，不计算新的同源关系或翻译 DNA。

## phase、不完整末端与例外 {#phase}
同一 CDS 的片段共享 `feature_id`，位于同一序列且链已知。GFF3 `phase` 从各段
生物学 5' 端计数到下一个密码子边界。只有第一段 phase 被跳过；后续 phase 校验
累计编码长度，保留跨片段密码子的碱基。显式 `segment_order` 可声明跨原点连接。

默认拒绝不完整末端。显式 `partial = "trim"` 排除不能组成完整密码子的碱基，
残基一表示第一个完整观测密码子。投影中遗传密码 ID 是元数据，不是翻译选择。
已声明的移码、滑移、翻译例外和不支持的位置运算符会被拒绝，不能推断未观测蛋白的
长度或相似度。

## 逆向映射 {#inverse}
使用 `project_ggchord_features(genomic, cds, from = "genome")` 反向投影。
非密码子边界需要明确 `codon_boundary` 策略。结果保留源坐标、父对象与投影来源。
处理不完整或多段生物学模型前，应阅读完整[函数参考](reference/project_ggchord_features.html)。
