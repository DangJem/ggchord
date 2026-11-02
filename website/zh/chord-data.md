先认识表格，再决定使用哪些图层。弦图是ggchord最早的核心视图；早期工作流围绕BLAST比对表展开，因此保留了紧凑的BLAST风格字段。这里的“多序列比对弦图”展示多条序列间已有的区间比对，绘图本身不执行多序列比对算法。

## 三张表的最小形状 {#tables}
| 表 | 必需列 | 每行的含义 |
| --- | --- | --- |
| seq_data | accver, length | 一条序列及其完整长度，bp |
| ribbon_data | qaccver, saccver, qstart, qend, sstart, send | 两条序列上的一对已有比对区间 |
| gene_data | accver, start, end, strand | 一个有已知正负链的基因区间 |

只画骨架时，后两张表都不需要。通用注释可保留未知链，请使用 `geom_feature()`，不要冒充有方向的基因。

## 名称与可选列 {#names}
`accver` 源于 accession.version，但值可为任意稳定字符串，例如自己的contig名称；不会自动添加或删除版本号。`q` 是query，`s` 是subject，所以 `qstart` 是query端起点。两个标识都必须能在序列表中找到。

`pident` 是百分比相似度，仅映射或筛选相似度时需要；比对表的 `length` 是比对长度，序列表的 `length` 是完整序列长度，二者不能混用。`evalue`、`bitscore`、`qlen`、`slen` 按分析需求保留。基因的 `anno` 是可选注释文字；要画标签时提供它。其他类别列可供 `aes()` 映射，不是强制输入。

## 查看真实表格 {#inspect}
```r
library(ggchord)
library(dplyr)
seq_data_example |> select(accver, length) |> slice_head(n = 3)
ribbon_data_example |> select(qaccver, saccver, qstart, qend, sstart, send) |> slice_head(n = 3)
gene_data_example |> transmute(accver, start, end, strand = source_strand, anno) |> slice_head(n = 3)
```

这些行来自随包真实记录；原始生物学链在该历史教学对象的 `source_strand` 中，因此显式恢复。自己的GFF3读入结果直接使用其 `strand`。

## 坐标与方向 {#coordinates}
使用一基闭区间bp坐标。基因 `start <= end`，方向另存于 `strand`。连接端点次序可表达反向比对；不要一律排序。未知完整长度不能由最大特征终点推断。详细契约见[数据与坐标](data.html#coordinates)。

继续：[1.2 从分析结果到输入表](chord-import.html)。
