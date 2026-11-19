{{code:courseblockan}}


示例使用随包表或随安装提供的小文件。每段代码列出所用的依赖。

<a id="adapters"></a>

## 把已有表格转换成需要的列 {#adapters}
这里使用随包的实际序列、连接和基因，
只重命名列来演示角色映射；没有额外生成示意序列或比对。

{{code:adapters}}

{{console:adapters}}

转换返回普通数据框，并记录 `ggchord_conversion` 信息。标识原样保留；缺少完整
序列长度时应补充可靠来源，不能把最大注释端点当作全长。核酸与蛋白比对的最小列
和单位见[导入比对结果](alignments.md)。

<a id="genbank"></a>

## GenBank：保留记录身份 {#genbank}
{{code:genbank}}

这个文件是 NCBI 的原始 L09137.2 记录，名称是 pUC19c；不要把它强行改成随包 pUC19。
返回值包括 `seq_data`、`features$annotation` 与 `report`，分段位置及记录/注释来源
保留。这里需要所选 Python 已安装 Biopython；函数不自动安装依赖。
`unsupported = "report"` 将不支持的位置记入报告，不把远端或复杂位置猜成区间。

本节需要可用的 Python Biopython。成功读取后，序列表含一条记录 L09137.2，完整长度为 2,686 bp。

显式 `python` 路径保留选定虚拟环境，包括解释器符号链接；不安装或改动环境。



<a id="files"></a>

## FASTA 长度与 GFF3 注释 {#files}

{{code:frontimportseq}}
{{output:frontimportseq}}



FASTA 给出完整长度。GFF3 返回区间、类型和原始链，绘图前可按需保留 CDS。

{{code:frontimportgff}}

{{output:frontimportgff}}

<a id="blast"></a>

## 读取 BLAST 表 {#blast}

{{code:frontimportblast}}

{{output:frontimportblast}}

每一行是一对区间。用于同一张图时，先核对两端 ID 与骨架 ID。这里的 pUC19 FASTA 用于演示长度读取，BLAST 表来自另一组序列。见[比对输入](alignments.md)。

<a id="alternatives"></a>

## 直接转换 Bioconductor 对象 {#alternatives}

已安装 GenomicRanges 时，可从一条区间转换成特征表。已知完整长度写在 `seqlengths` 中，可另用 `as_ggchord_seq(gr)` 提取。

{{code:frontimportgr}}

{{output:frontimportgr}}

GInteractions 的第一端是 query，第二端是 subject。下面的最小输入用两个区间构成一条连接。此示例另需 InteractionSet。

{{code:frontimportgi}}

{{output:frontimportgi}}

其他输入对象见[导入与转换](import-intervals.md#ranges)。用转换后的表[添加图层](multi-sequence-start.md#draw)。
