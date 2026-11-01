<div class="task-intro"><p>将已有表格、BED、PAF 或 GenBank 转成 ggchord 输入，保留记录身份，并且只转换一次坐标约定。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。 涉及文件时，先解压[教学输入](../assets/tutorial-inputs.zip)，使工作目录包含 `examples/`。

{{prepare:setup,docstyle,plasmiddata}}

## 把已有表格转换成需要的列 {#adapters}
这里使用随包的实际序列、连接和基因，
只重命名列来演示角色映射；没有额外生成示意序列或比对。

{{code:adapters}}

转换返回普通数据框，并记录 `ggchord_conversion` 信息。标识原样保留；缺少完整
序列长度时应补充可靠来源，不能把最大注释端点当作全长。核酸与蛋白比对的最小列
和单位见[导入比对结果](alignments.html)。

## BED 和 PAF：坐标只转换一次 {#zero-based}
BED 例子把已知 pUC19 注释写为最小 BED3，再读回；PAF 直接读取公开教程中的真实
minimap2 输出，保留原始可选标签和值。显示筛选保留 mapq ≥30 的记录，绘图只接收六个端点角色。

{{code:formats}}
{{figure:paf}}

BED 的零起点半开区间由读取器转成一基包含坐标；不要再手工加一次 1。
BED3 没有链列，读回时链为未知；需要链方向时应使用含 strand 的原始格式。
PAF 的负相对方向表现为递减的 subject 端点，这是合法反向命中，不能一律排序去掉方向。
比对质量和 `pident` 来自生产工具报告的字段，不在绘图时重新计算碱基比对。

## GenBank：保留记录身份 {#genbank}
{{code:genbank}}

这个文件是 NCBI 的原始 L09137.2 记录，名称是 pUC19c；不要把它强行改成随包 pUC19。
返回值包括 `seq_data`、`features$annotation` 与 `report`，分段位置及记录/注释来源
保留。这里需要所选 Python 已安装 Biopython；函数不自动安装依赖。
`unsupported = "report"` 将不支持的位置记入报告，不把远端或复杂位置猜成区间。

## 可选的区间对象 {#ranges}
{{code:optionalranges}}
{{results:optional}}

这里选择 pUC19 中原本未知链的注释；转换后的 `*` 仍表示未知。
完整长度来自已知序列长度和 Seqinfo。`GInteractions` 可以转换成连接，但仍需明确
两端序列及长度。缺少可选依赖时结果会标记跳过，不把源码存在当作运行通过。
详见[数据格式说明](reference/ggchord_data_contract.html)。当前不直接解析 SAM/BAM 的
CIGAR；需先由相应工具正确解析，再通过表格角色映射转换，不能仅靠改列名处理剪接。

## 完成后检查 {#check-result}
先查看返回列与转换报告，再选绘图列。PAF 反向命中保留端点方向；GenBank 记录身份不能因为展示需要而改成另一条随包质粒。

<div class="next-step">继续：<a href="alignments.html">导入核酸和蛋白比对 →</a></div>
