<a id="source"></a>

## 取回一个具名图层的输入 {#source}

用 `layers = "observed"` 选中具名图层，再请求 `stages = "source"`。来源表包含进入这一层的三条记录。要将其与显示几何对照，可请求 `stages = c("source", "geometry")` 并提供目标设备和尺寸。

{{example:inspect-source}}

<a id="checks"></a>

## 根据现象选择核对项 {#checks}

| 问题 | 查看什么 |
| --- | --- |
| 一条连接带去了哪里？ | 核对来源标识、端点和显示窗口边界。参考序列重复显示时，检查或指定两端的显示实例。 |
| 哪些区间贡献了窗口覆盖率？ | 导出 `statistic` 阶段，查看 `source_rows`。[逐窗示例](tracks-statistics.md#members)同时列出实际分母。 |
| 标签为何隐藏或拥挤？ | 给 `diagnose_ggchord()` 指定输出设备和最终宽高，结合返回的 `issues`、测量状态和实际图形核对。 |
| 同一来源位置投影到了哪个面板？ | 分面之间的投影有差异时，可在 `map_ggchord_position()` 指定目标 `panel`。 |

查看结果表时，同时核对 `stage_status`。可用阶段也可能有零行；失败或不可用的阶段会给出原因。比较不同尺寸下的位置之前，先核对 `measurement`。
