<a id="gff3"></a>

## 从 GFF3 读取注释区间 {#gff3}

随包 GFF3 记录了 MT108731.1 的 CDS 注释。位置、链向和原始行号应一起保留，绘图后可以追溯到单条记录。遇到链向未知的记录时，`unstranded = "keep"` 会保留其原始状态。

{{example:installed-gff}}

<a id="reference"></a>

## 补齐完整参考 {#reference}

用 `read_fasta_lengths()` 读取 FASTA 的 ID 和全长。一条序列的最后一个注释可能远在序列终点之前。对于此例，可以从 `seq_data_example` 选择 MT108731.1 得到完整参考，再把 GFF3 表交给 `geom_feature(data = features)`。换成自己的文件时，应核对两个表的 ID 拼写是否完全一致。

GenBank 记录同时带有序列和注释。`read_genbank()` 返回的 `seq_data` 是参考表，`features$annotation` 是特征行。用编码特征进行[坐标投影](projection.md)时，保留拼接片段及其来源链。
