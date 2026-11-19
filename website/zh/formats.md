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

从完整参考表开始构图。各层显式提供所用观测表，或用 `prepare_ggchord_data()` 组装准备对象。为了展示而筛选记录时，保留 accession 与原端点。需要核对某个图层或统计窗口对应的源记录时，见[查询布局与来源](diagnostics.md)。设置最终图片尺寸见[预览与导出](output.md)。
