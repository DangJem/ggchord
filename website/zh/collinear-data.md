共线性图是ggchord从弦图扩展出的直线视图。设计目标是保留共同数据和图层，尽可能通过替换坐标系完成转换。

## 相同的三张表 {#tables}
| 输入 | 必需列 | 共线性含义 |
| --- | --- | --- |
| seq_data | accver, length | 每条序列占一行，共享bp横向比例 |
| ribbon_data | qaccver, saccver, qstart, qend, sstart, send | 两行之间的已有区间对应 |
| gene_data | accver, start, end, strand | 行上基因，正负链由Position上下分列 |

可选列和BLAST命名来历与[1.1](chord-data.html)相同；相似度需要pident，标签需要anno。局部窗口的source字段用于保留原始坐标，不能将显示间隙解释为缺失DNA长度。

## 数据如何获得 {#import}
[1.2](chord-import.html)的FASTA/GFF3/BLAST流程可以直接复用；已有PAF、GRanges、GInteractions通过[适配器](formats.html)读入。蛋白比对必须先[投影到已知CDS坐标](projection.html)，不能把aa直接画成bp。

## 当前转换边界 {#contract}
[下一页](collinear-build.html)给出只替换coord的实际例子。显式圆形半径、曲率、ring和径向序列标签位置在当前共线性合同中会被拒绝，而不是悄悄忽略。将这些坐标特有样式留给各自展示图；“任意组合完美转换”仍是长期目标，不是当前保证。

共线性默认只连接相邻显示行；要保持弦图的连接集合，需要在共同图层明确 `adjacent_only = FALSE`。显示方向可手动指定或交给coord自动判断，两种入口不能同时使用。
