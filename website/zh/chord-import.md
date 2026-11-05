解压[教学输入包](../assets/tutorial-inputs.zip)，使工作目录包含 `examples/`。每次先读一种文件并检查返回的几行，再按序列 ID 对齐输入。

## FASTA 长度与 GFF3 注释 {#files}

{{code:frontimportseq}}

{{output:frontimportseq}}

FASTA 给出完整长度。GFF3 返回区间、类型和原始链，绘图前可按需保留 CDS。

{{code:frontimportgff}}

{{output:frontimportgff}}

## 读取 BLAST 表 {#blast}

{{code:frontimportblast}}

{{output:frontimportblast}}

每一行是一对区间。保留两端坐标顺序，核对两列序列 ID 能在 FASTA 中找到。如何生成标准 BLAST 文件见[比对输入](alignments.html)。

## 直接转换 Bioconductor 对象 {#alternatives}

已安装 GenomicRanges 时，可从一条区间转换成特征表。已知完整长度写在 `seqlengths` 中，可另用 `as_ggchord_seq(gr)` 提取。

{{code:frontimportgr}}

{{output:frontimportgr}}

GInteractions 的第一端是 query，第二端是 subject。下面的最小输入用两个区间构成一条连接。此示例另需 InteractionSet。

{{code:frontimportgi}}

{{output:frontimportgi}}

其他输入对象见[导入与转换](formats.html#ranges)。用转换后的表[添加图层](chord-build.html)。
