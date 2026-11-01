## 用实际编码序列建立分组 {#membership}
先运行[共享准备代码](index.html#quick-start)。这里选择随包 pUC19 与 pETDuet-1 的
AmpR 注释，提取实际 DNA 并按链方向校正。只有编码 DNA 完全相同，才放入同一组，
不会因为注释都叫 AmpR 就直接认定同源。

{{code:homologs}}
{{figure:homologs}}

这个严格例子展示已知完全相同的 CDS；不代表所有同源基因都应完全一致。
通常的同源分组应来自可靠的外部分析。查看 `homologs$report$comparisons` 和
`$members`；未知成员和歧义 paralog 会被报告，默认不随意连接所有组合。
`paralogs = "all"` 才明确请求受界限约束的笛卡尔积。跨分段注释的外部包围区间
是锚点，不是外显子逐碱基比对。`find_ggchord_loci()` 由确切 ID 定位并按生物学方向扩展。

## 将有证据的蛋白区间投影到两端 CDS {#projection}
{{code:linkprojection}}
{{results:projection}}

这里投影完全相同 CDS 的前 20 个残基，因此 `ungapped = TRUE` 有序列依据。
普通 BLASTp HSP 不能仅凭两端跨度相等就作出这个断言。工具按 CDS 分段投影，并
保留源蛋白统计，不把蛋白相似度冒充核苷酸相似度。真实 BLASTp 查询范围的投影
见[导入核酸和蛋白比对](alignments.html#protein)，单侧注释的完整流程见[坐标投影](projection.html)。
