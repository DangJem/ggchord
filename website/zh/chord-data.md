`ggchord()` 从序列表确定有哪些骨架。连接和特征按需加入。下列几行来自随包数据，运行代码可以看到同样的结果。

## 序列表 `seq_data` {#tables}

{{code:frontseqtable}}

{{output:frontseqtable}}

| 必需列 | 含义 |
| --- | --- |
| `accver` | 唯一序列 ID。可以用自己的 contig 名称，其他表按该值匹配。 |
| `length` | 完整序列长度，单位为 bp。 |

**可选列**可用于显示或分组，例如另存一个简短的显示名称。

## 连接表 `ribbon_data` {#inspect}

{{code:frontlinktable}}

{{output:frontlinktable}}

| 必需列 | 含义 |
| --- | --- |
| `qaccver`、`saccver` | 查询端（query）和目标端（subject）的序列 ID，都须出现在序列表中。 |
| `qstart`、`qend` | 查询序列上的两个区间端点。 |
| `sstart`、`send` | 目标序列上的两个区间端点。 |

**可选列** `pident` 是比对一致性百分比，可用于筛选或配色。例如，首行约为 98.6%。比对表也可以有 `length`，它表示比对长度。

## 基因与通用特征 {#names}

{{code:frontgenetable}}

{{output:frontgenetable}}

| 列 | 用途 |
| --- | --- |
| `accver` | 必需。该区间所在的序列 ID。 |
| `start`、`end` | 必需。区间的起点和终点，`start <= end`。 |
| `strand` | 基因箭头需要 `+` 或 `-`。无方向或方向未知的注释可用通用特征图层。 |
| `anno` | 可选。可以作为标签文字。 |
| `type` | 可选。通用特征的类别，例如 CDS、启动子或复制起点。 |

这里从随包示例的 `source_strand` 取出原始注释链。自己导入的 GFF3 数据可直接使用其中的 `strand`。

## 坐标怎么读 {#coordinates}

输入使用从 1 开始、两端均计入的 bp 坐标。例如，10 到 12 包含 3 个碱基。基因区间用 `strand` 表达方向，连接区间则可能通过端点的升降顺序表达方向，导入后应保留两端各自的顺序。

继续阅读[绘图基本概念](modes.html)。需要从文件读取时可查阅[输入格式与转换](formats.html)。

字段采用 BLAST 风格，是 ggchord 从比对可视化发展而来的命名沿革。query/subject 对应 q/s 前缀；自定义 ID 也可以使用。实际 BLAST 表见[比对输入](alignments.html)。

## 蛋白序列的单位 {#protein-unit}

核酸表可省略单位列，默认使用 bp。完整蛋白长度表应声明 `sequence_unit = "aa"`，同一图使用一个单位。[蛋白比对示例](alignments.html#protein)展示完整输入和配图。
