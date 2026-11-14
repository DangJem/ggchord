<a id="choose-format"></a>

## 手中的文件从哪里开始 {#choose-format}

| 手中的文件 | 下一步 |
| --- | --- |
| FASTA 序列 + GFF3 注释 + BLAST 局部比对 | [从源文件画比较图](real-data.md#read) |
| BED 区间或 PAF（Pairwise mApping Format）比对 | [核对零基坐标](import-intervals.md#zero-based) |
| GenBank 单序列及注释 | [保留原始记录与分段](import-sequences.md#genbank) |
| bedGraph 或 BigWig 连续信号 | [显式保留值和缺失区间](import-signals.md#signals) |
| VCF/BCF 变异和结构断点 | [保留原始变异证据](import-relations.md#variants) |
| 已有蛋白 signature 注释 | [保留 AA 位置及源分数](import-relations.md#protein-signatures) |
| 已有 feature 对及组 | [保留原关系与成员](import-relations.md#feature-relations) |
| BAM/CRAM 已比对 reads | [使用成熟覆盖度桥接](import-signals.md#aligned-reads) |
| 已有 delta 结果 | [区分原生统计和坐标单位](import-relations.md#delta) |
| 已有 MAF 或 chain 基因组比对 | [保留块及源坐标](import-relations.md#genome-alignments) |
| 已有 R 表格 | [映射必需字段](import-sequences.md#adapters) |
| GRanges 或 GInteractions | [直接转换区间对象](import-intervals.md#ranges) |

<a id="files"></a>

## 三文件到首图 {#files}

[FASTA / GFF3 / BLAST](import-sequences.md#files)
