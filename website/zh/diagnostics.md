<a id="source"></a>

## 取回一个具名图层的输入 {#source}

给连接带图层设置稳定 ID，再请求 `data_stage = "source"`。结果中的三条记录保留进入这一层的 accession 和原始端点。要同时查询来源和连接带几何，可再指定 `include = "ribbon"`。

{{example:inspect-source}}

<a id="checks"></a>

## 根据现象选择核对项 {#checks}

| 问题 | 查看什么 |
| --- | --- |
| 一条连接带去了哪里？ | 比较来源 ID、端点、已选择的序列和邻接规则。重复显示一条参考时，可能需要指定具名副本。 |
| 哪些区间贡献了窗口覆盖率？ | 查看轨道的 `track_statistics` 表和 `source_rows`。[逐窗示例](tracks-statistics.md#members)同时列出实际分母。 |
| 标签为何隐藏或拥挤？ | 给 `diagnose_ggchord()` 指定输出设备和最终宽高，结合 `label_diagnostics`、测量状态和实际图形核对。 |
| 同一来源位置投影到了哪个面板？ | 分面之间的投影有差异时，可在 `map_ggchord_position()` 指定目标 `panel`。 |

查询统计阶段时，先看 `layer_stage_status`。统计成功但结果为零行时状态为 `available`。计算失败时状态为 `failed` 并给出原因。两种情况对应的后续处理不同。
