<a id="fields"></a>

## 把各类记录匹配到参考序列 {#fields}

| 表 | 起步字段 | 含义 |
| --- | --- | --- |
| 参考表 | `accver, length` | 每个来源 ID 一行，记录完整序列长度。蛋白参考请声明 `sequence_unit = "aa"`。 |
| 特征表 | `accver, start, end` | 一个区间位于一条参考上。保留已知的 `strand`。稳定的 `feature_id` 让子集标签找到相应特征。 |
| 连接带表 | `qaccver, saccver, qstart, qend, sstart, send` | 位于两条具名参考上的两个区间。端点顺序带有比对方向信息。 |
| 信号表 | `accver, start, end` 加数值列 | 根据所选几何，把数值映射到高度或填色。 |

位置从 1 开始，区间包含两端。1–100 在 bp 参考上跨越 100 个碱基，在 aa 参考上跨越 100 个残基。保留源文件中的标识，包括 accession 的版本号。

<a id="prepared"></a>

## 把参考、特征和连接交给对应组件 {#prepared}

从 `ggchord(reference)` 开始，第二个参数是 `mapping`。其他表交给负责绘制它们的层，例如 `geom_feature(data = features)` 或 `geom_link_ribbon(data = hits)`。也可以先组装各表。例中取回准备好的特征表，运行 `ggchord(prepared) + geom_seq() + geom_feature()` 就能把这些特征画在参考上。

{{example:prepared-tables}}

教学注释表在 `source_strand` 保留来源链，示例明确选择这一列。如果区间按蛋白残基计数、目标参考按基因组碱基计数，需使用对应的 [CDS 坐标模型](projection.md)。
