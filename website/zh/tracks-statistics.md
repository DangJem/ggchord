<span id="interval-coverage"></span>

<a id="coverage"></a>

## 按组比较各窗口的覆盖率 {#coverage}

以 MT108731.1 为共同的 query 参考，按 subject accession 拆分实测 HSP。每个面板显示该 subject 在 query 连续 10,000 bp 窗口内覆盖的碱基比例。重叠的 HSP 取并集，每个碱基计一次。三个面板共用 0–1 高度轴，可以直接比较。

{{example:coverage-facets}}

<span id="count-and-col"></span>

<a id="members"></a>

## 读取实际分母与来源行 {#members}

完整窗口的分母为 10,000 bp。末窗在第 64,323 个碱基截止，分母为 4,323。对于 subject MT118296.1，20,001–30,000 bp 窗口内有 8,415 个碱基被覆盖。导出的 `source_rows` 列给出参与该窗统计的输入行。

{{example:window-members}}

对单位置观测，`geom_track_bar()` 统计各位置的记录数。如果输入已包含计数列，使用 `geom_track_col()` 将该列映射到高度。

<a id="domains"></a>

## 选择共用或局部高度域 {#domains}

若观测来自不同生物学样本，在区间表保留样本 ID，再按该列分面。默认共用高度域，方便跨样本比较。若要看各面板内的局部波动，在轨道数值标度设 `facet_domain = "panel"`，让数值上下限自动训练，保留可见的局部数值轴。两张独立图比较同一种度量时，可用 `shared_track_limits()` 算出范围，再传给两边的标度。

平均深度会分别计入重叠区间。要和覆盖率一同展示，可在第二条轨道映射 `aes(y = after_stat(mean_depth))`，为两种度量分别标出单位。
