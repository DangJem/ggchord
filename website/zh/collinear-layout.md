先学会[第一张共线性图](collinear-start.html)，再用同一组真实输入观察对齐、锚点和方向。无需学习弦图。

## 输入与默认左对齐 {#baseline}

本页独立加载随包四条完整序列及其真实比对。序列表必须有 `accver`、`length`；连接表必须有两个 ID 及 query/subject 的起止位置，单位为 bp，一基闭区间。保留反向端点，不先排序。

{{code:alignbase}}
{{figure:alignleft}}

`base_align` 保留数据与图层，后续仅替换坐标系。默认 `align = "left"`、`anchor = NULL`、`orientation = "manual"`。横向长度使用共同尺度，短序列不会被拉伸至相同长度。

## 只改变对齐方式 {#alignment}

{{code:aligncenter}}
{{figure:aligncenter}}

四条完整序列的中点对齐。每条序列的显示平移不同；比对端点在其骨架上的位置没有改变。

{{code:alignright}}
{{figure:alignright}}

右端对齐。三幅图均为 9 × 6 英寸，数据、图层、行顺序与方向相同。比较的是横向平移，不是新的比对结果。

## 给每条序列一个锚点 {#anchor}

{{code:alignanchor}}
{{figure:alignanchor}}

与左对齐基线相比，仅增加 `anchor`。这个命名向量以 `accver` 为键，数值是源序列中点，因而在本例与完整序列居中有相同效果。它只是几何示范；科研中应从已确认的共同位点构造锚点。

也可给 `accver, position` 或 `accver, start, end` 数据框。每条显示序列恰有一个锚点，缺失或重复会报错。区间锚点由 `align` 选择显示左端、中点或右端；这是应用显示方向后的边缘。多拷贝基因须先明确选择哪一份。

## 自动选择显示方向 {#orientation}

{{code:alignauto}}
{{figure:alignauto}}

与左对齐基线相比仅启用 `orientation = "auto"`，`min_support = 0.1` 显式重述默认阈值。算法按已有连接对前面显示行的支持选择方向，第一条序列固定参照。阈值为正反支持差的绝对值除以支持总量，范围 0–1；平局、支持不足或无连接时保留初始方向。若图与手动基线相同，可能是当前证据支持初始方向，不代表自动设置没有执行。检查 `export_ggchord_layout()` 的诊断，不能由图看起来整齐就认定方向正确。

| 概念 | 数据或设置 | 意义 |
| --- | --- | --- |
| 原始序列方向 | FASTA 的记录顺序 | 输入坐标的参照 |
| 基因生物学链 | 注释的 `strand` | 相对原序列的转录方向 |
| 比对方向 | query/subject 端点递增或递减 | 两区间的相对匹配方向 |
| 显示方向 | `seq_orientation` 或自动方向 | 图中如何摆放序列 |

显示翻转不改写原始链，也不能单独作为真实基因组倒位的结论。手动设置 `geom_seq(seq_orientation = ...)` 时不要再启用自动方向。

## 参数选择与组合边界 {#parameters}

| 参数 | 默认值 | 类型、作用与限制 |
| --- | --- | --- |
| `align` | `"left"` | left/center/right；对齐完整序列或锚区间 |
| `anchor` | `NULL` | 命名数值向量或表；不裁剪数据 |
| `orientation` | `"manual"` | manual/auto；自动与显式序列方向冲突 |
| `min_support` | `0.1` | 0–1 的数值；控制方向证据阈值，不是 identity 门槛 |
| `segments` | `NULL` | 显式窗口表；不能与 anchor 或 auto 同用 |
| `segment_gap` | `100` | 正数，按 bp 显示尺度预留断点空间，不是缺失长度 |

这些设置属于 `coord_collinear()`，不传给其他坐标模式。分段时由窗口表掌管方向，详见[局部窗口](views.html#segments)。圆弧半径与曲率不能用于此模式。

## 练习和自己的数据 {#adapt}

先只把 `align` 从 left 改成 right，观察短序列与最长序列的相对位置；再给每条序列选择一个有业务依据的锚点，解释为何可以比较。使用自己的 FASTA 长度和 BLAST/PAF 端点时先核对 ID、单位、完整长度及反向端点，见[数据导入](real-data.html#read)。

没有连接时先查端点和相邻行筛选，见[跨行连接](collinear-routing.html)。接下来学习[局部窗口与分段](views.html)，完整参数见 [`coord_collinear()`](reference/coord_collinear.html)。
