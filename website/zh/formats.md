## 从标准角色开始 {#adapters}
先运行[共享准备代码](index.html#quick-start)。列重命名是明确的角色映射，不要求改写
accession version，也不能借此推断完整序列长度。

{{code:adapters}}

检查返回的普通表格和 `ggchord_conversion` 属性。`as_ggchord_seq()` 要求已知完整长度，
`as_ggchord_feature()` 保留未知链，`as_ggchord_link()` 保留锚点及已有元数据。
图层接受的直接 S4 输入使用同一组适配器。

## BED 与 PAF：只转换一次 {#zero-based}
这里的文件均为显式合成输入，用来演示格式转换，不代表真实序列上的实测比对或注释。

{{code:formats}}

BED `[0,100)` 转为一基闭区间 `1..100`。PAF 相对负方向产生递减的 subject 端点，
mapq 255 表示不可用，原值与转换后的缺失值分别保留。PAF 相似度来自工具报告计数，
不重新按碱基计算。保留导入属性和物理来源行，绘图前不要再加一或减一。

## GenBank：记录、分段与来源 {#genbank}
{{code:genbank}}

只有所选 Python 已安装 Biopython 时才执行这个可选例子，函数本身不安装环境。
结果含 `seq_data`、命名的 `features$annotation` 和 `report`。joined location
保留真实片段及记录/feature 来源索引。未知长度不能由 feature 末端推断。
`unsupported = "report"` 把不支持位置保留在报告中，而不是猜测并压平为本地区间。

## 可选区间对象 {#ranges}
{{code:optionalranges}}
{{results:optional}}

未知链 `*` 保持未知。完整长度来自 Seqinfo，而非最大观测区间。
InteractionSet 的 `GInteractions` 还可适配为连接，须保留统一的区域/长度背景及
显式生物学归属。缺少可选依赖的例子在记录中标为跳过，不声称已经验证。
详见完整[数据契约](reference/ggchord_data_contract.html)。
