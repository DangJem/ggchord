先有序列与分析结果，再由ggchord读取；没有比对表时不应为了连线伪造同源关系。

## 从FASTA与GFF3开始 {#files}
解压[教学输入包](../assets/tutorial-inputs.zip)，使当前目录包含 `examples/`。FASTA提供完整长度，GFF3提供区间和链；保留原始记录，在显示表中选择必要字段。

{{prepare:sequencepalette,docstyle}}
{{code:showcasedata}}

`show_alignments` 保留32条原始HSP，`show_links` 明确筛为12条；`show_genes` 保留全部266条CDS，`show_features` 按长度至少500 bp、序列内注释顺序分四组各选一个，并补充有来源的负链首尾注释用于展示。筛选不是重跑比对，也不证明功能同源。

坐标排查时可用tidyr将四个端点列转成长表，原始连接表保持原样：

{{code:inspectends}}

## 用BLAST产生自己的结果 {#blast}
以下命令是分析者在安装BLAST后对自己的文件执行的生产步骤，网站构建不运行它们。

```sh
blastn -query query.fna -subject subject.fna -outfmt 6 -out result.tsv
```

```r
library(ggchord)
blast_source <- read_blast(file = "result.tsv")
```

无自定义字段列表的 `-outfmt 6` 为标准12列表，7另带注释；完全不写outfmt是其他格式。query/subject ID由读取器映射到qaccver/saccver。参见[NCBI字段说明](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/table/appendices.T.options_common_to_all_blast/)和[核酸/蛋白比对教程](alignments.html)。随站历史HSP只是格式转换，没有冒充本轮新BLAST结果。

## 已有其他格式时 {#alternatives}
| 来源 | 入口 | 注意 |
| --- | --- | --- |
| 已读表格 | as_ggchord_seq/feature/link | 显式mapping映射角色，保留原表 |
| BED、PAF | read_bed、read_paf | 读取器只转换一次零基半开坐标 |
| GenBank | read_genbank | 需要可用Biopython；保留记录与分段来源 |
| GRanges、Seqinfo | as_ggchord_feature、as_ggchord_seq | 完整长度来自Seqinfo，未知链不转成正链 |
| GInteractions | as_ggchord_link、as_ggchord_seq | 第一锚点为query；仍需可靠完整长度 |

这些适配器对应[用户issue #2](https://github.com/DangJem/ggchord/issues/2)提出的需求。当前代码已有入口；issue中的历史“开发版未经测试”不等于现阶段所有组合均已验证。本地缺失可选包的路径会明确跳过。[格式专题](formats.html#ranges)给出可选对象的完整代码和限制。

继续：[1.3 从骨架逐层画起](chord-build.html)。
