<a id="alignments"></a>

## 导入 HSP 端点 {#alignments}

读取 MT108731.1 对 MT118296.1 的随包 BLAST 表。每条 HSP 是一个局部比对片段。四个端点指明两侧区间，`pident` 表示一致性百分数。比对文件若缺少序列全长，需要从序列表取得。

{{example:installed-blast}}

若输入为 PAF，`read_paf()` 还能返回 query 与 subject 的全长。读取时会将文件边界换算为一基闭端点，并保留反向匹配。这些区间可直接交给[连接带图层](multi-sequence-links.md)。如果输入是 [BEDPE](../../man/read_bedpe.Rd)、[chain](../../man/read_chain.Rd) 或 [delta](../../man/read_delta.Rd)，可查阅相应接口，确认关系表的返回位置和解码依赖。

<a id="members"></a>

## 读取已知分组的成员 {#members}

分组文件记录哪些特征属于同一组。本例首个具名组包含四个植物特征。后续若要匹配 GFF3 注释，请保留 ID 和原始行号。绘制两两关系时，需要已有的配对，或明确指定组内要比较的成员。

{{example:installed-members}}
