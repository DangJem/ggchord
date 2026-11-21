<a id="inputs"></a>

## 按手头的数据找入口 {#inputs}

| 手中的数据 | 需要得到的结果 | 从哪里开始 |
| --- | --- | --- |
| 完整序列长度、特征区间或双端记录 | 标识可匹配、单位明确的绘图表格 | [参考表与图层表](data.md) |
| FASTA、GenBank 或 GFF3 | 长度表或注释区间 | [读取序列与注释文件](import-sequences.md) |
| BLAST、PAF 或特征分组文件 | 比对端点或明确记录的组成员 | [读取关系](import-relations.md) |
| bedGraph 或 BigWig | 区间及其测量数值 | [读取连续信号](import-signals.md) |
| 氨基酸区间及其父 CDS | 对应的基因组碱基与链向 | [将残基投影到 CDS](projection.md) |

<span id="reader-results"></span>

<a id="results"></a>

## 从源记录追到图上 {#results}

用 `ggchord()` 登记完整参考表，再把观察数据显式传给各图层。用 `ggchord_data()` 将相关的参考表、特征表和连接表组成一个集合。筛选记录时保留标识及原始端点。用 `validate_ggchord_data()` 检查输入角色，用 `data_report()` 查看读取器或数据操作保留的来源证据。[布局与来源查询](diagnostics.md)可核对图层的输入，[预览与导出](output.md)说明如何指定图形的最终尺寸。
