<a id="choose-format"></a>

## 手中的文件怎样读入？ {#choose-format}

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

## 从文件到图 {#files}

[读取 FASTA、GFF3 和 BLAST](import-sequences.md#files)，查看[完整的序列比较示例](real-data.md#read)。如果数据已经在 R 中，可直接核对[必需字段和坐标](data.md#tables)。

## 导入结果怎样接到图层 {#reader-results}

| 格式与返回对象 | 坐标与下游入口 | 额外解码依赖 |
| --- | --- | --- |
| FASTA 长度表，GFF3 注释表，BLAST 连接表 | FASTA 声明 bp 或 aa。GFF3 使用一基闭区间，绘方向时核对来源链。BLAST 保留两端顺序与原单位。长度表作为 `ggchord()` 的参考。注释和比对各自传给对应图层的 `data`。 | 包内读取 |
| BED、BEDPE、PAF、bedGraph 均返回表 | 原始零基半开区间转为一基闭区间。依次用于特征、成对连接、局部比对、数值轨道。 | 包内读取 |
| BigWig 数值表 | 解码结果已是一基闭区间，`value` 用于数值轨道。 | `rtracklayer`，Windows 可使用已有 bedGraph |
| GenBank 列表 | 全长表和 `features$annotation` 可交给准备函数。一基特征保留分段来源。 | Python 与 Biopython |
| 对齐 FASTA 长表，MAF 或 chain 列表 | 前者保留比对列和来源位置。后两者提供全长与配对区间，MAF 另含字符。报告保留解码信息和原文，大文件导入时这些对象会同时占用内存。 | 对齐 FASTA 由包内读取，MAF 与 chain 使用 Python、Biopython |
| delta 列表 | 全长与连接按解码器坐标和声明单位交给准备函数。 | 已安装的 `show-coords` |
| VCF 或 BCF 列表 | 使用一基参考等位区间 `features$variant`，显式解析的断点和连接单独取用。 | `VariantAnnotation`，BCF 另需 pysam，结构变异扩展依赖见[变异输入](import-relations.md#variants) |
| signature 或 feature 关系列表 | 前者按蛋白 aa 坐标给出 `features$signature`。后者保留原始 pair 或组成员，只有完整定位的原始 pair 才有可绘连接。 | signature JSON 需 `jsonlite`，本地表格可直接读取 |

## 导入后要处理什么？ {#data-operations}

- [校验标识、全长和单位](data.md#tables)
- [筛选、聚焦或整理比较数据](workflows.md)
- [在核酸与蛋白坐标之间投影](projection.md)
- [检查来源和布局结果](output.md#layout-results)
- [查看旗舰图的原始文件与构建过程](flagship-data.md)
