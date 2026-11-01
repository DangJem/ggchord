## 成员归属是输入证据 {#membership}
先运行[共享准备代码](index.html#quick-start)。已声明的合成 cluster 将每条序列一个
feature 相连，用于教学，不由相同的基因注释文本推断关系。

{{code:homologs}}
{{figure:homologs}}

检查 `homologs$report$comparisons` 和 `$members`。未知成员会遗漏并报告，含糊的
paralog 默认不连接；`paralogs = "all"` 显式请求有上限的笛卡尔积。
重复 feature ID 需要明确分段策略，外包范围是锚点，不是外显子比对。
`find_ggchord_loci()` 按精确 feature ID 选择，并沿生物学方向扩展；未知方向要求
对称扩展。下一步可做[局部聚焦](views.html)。

## 投影已声明的蛋白对应 {#projection}
{{code:linkprojection}}
{{results:projection}}

蛋白端点跨度相等不构成无 gap 证据，必须显式声明 `ungapped = TRUE`。
任一 CDS 跨编码片段边界时，辅助函数拆分对应片段并保留源蛋白元数据。
它不把蛋白相似度作为核苷酸相似度，也不生成核苷酸比对长度。
domain、不完整末端及 phase 策略见 [CDS 投影](projection.html)。
