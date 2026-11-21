<a id="gff3"></a>

## 从 GFF3 读取注释区间 {#gff3}

随包 GFF3 记录了 MT108731.1 的 CDS 注释。位置、链向和原始行号应一起保留，绘图后可以追溯到单条记录。遇到链向未知的记录时，`unstranded = "keep"` 会保留其原始状态。

{{example:installed-gff}}

<a id="reference"></a>

## 补齐完整参考 {#reference}

用 `read_fasta_lengths()` 读取 FASTA 的 ID 和全长。一条序列的最后一个注释可能远在序列终点之前。对于此例，可以从 `seq_data_example` 选择 MT108731.1 得到完整参考，再把 GFF3 表交给 `geom_feature(data = features)`。换成自己的文件时，应核对两个表的 ID 拼写是否完全一致。

GenBank 记录同时包含序列和注释。安装 Biopython 后，`read_genbank()` 会在同一集合中返回 `sequences` 和 `features$annotation`。用编码区做[坐标投影](projection.md)时，请保留拼接片段及原始链向。对于 BED 区间或 VCF 变异，见 [read_bed()](../../man/read_bed.Rd) 和 [read_vcf()](../../man/read_vcf.Rd) 的坐标及返回表说明。
