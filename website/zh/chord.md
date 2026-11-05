<div class="lesson-intro">把四条基因组放在同一张关系图中，组合直线与弧形骨架、保留真实链方向的基因、连续相似度标度和精选标签，同时保留数据来源。</div>

## 综合结果 {#overview}
{{figure:showcasechord}}

先比较深蓝色 `MT118296.1` 和青绿色 `OQ646790.1`，它们的半径分别为 0.82 和 1.18。橙色 `OR222515.1` 和红色 `MT108731.1` 使用 1.02 和 1.08。半径用来拉开骨架层次，序列长度仍按 bp 刻度读取。[单独比较半径效果](chord-geometry.html#shape)。

每条序列按已有的正、负链分别抽样，无 CDS 长度门槛。两条 MT 序列没有负链 CDS，OQ646790.1 唯一的负链 CDS 为 207 bp。本图用显式 Position 把正链箭头放在内侧、负链箭头放在外侧。箭头方向、颜色和数据仍使用原始链。正链抽样优先保留与所显示比对端点重叠的记录，以展示连接带局部避障。每条箭头使用原始注释作为标签，包括 hypothetical protein。

## 准备原始记录 {#prepare}
将[教程输入包](../assets/tutorial-inputs.zip)解压到工作目录，使其包含 `examples/`。使用 `install.packages(c("dplyr", "purrr", "tidyr"))` 安装示例可选依赖；完整流程用它们整理表格和计算 GC 轨道。

准备脚本包含原始数据读取、筛选和显示设置。将脚本放在当前工作目录。

{{recipe:chord}}

输入包含四条未改写的 FASTA、对应 GFF3 注释和 32 条原始 HSP。展示保留长度至少 900 bp、相似度至少 75% 的 HSP，每对序列最多三条，以减少拥挤；没有合并或编造比对。

## 组合图层 {#compose}
端点重叠抽样使用 dplyr ≥ 1.1.0。

{{code:showcasechord}}

| 设置 | 在这张图里的作用 |
| --- | --- |
| `seq_orientation`、`seq_curvature`、`seq_radius` | 按序列调整显示方向、弯曲程度和半径。曲率为 0 时画直线。 |
| `seq_gap = 0.06`、`rotation = 50` | 增加骨架间距，并旋转整体布局以摆放图中元素。 |
| `link_avoid = "smooth"`、`alpha = 0.46` | 让比对带在箭头附近平滑留白，并以半透明显示交叠。 |
| `position = chord_position`、`width = 0.025` | 按真实链把箭头放在骨架两侧，偏移为 ±0.030。较窄的箭头贴近骨架，标签使用同一位置。 |
| `label_layout = "radial"`、`label_wrap = 24` | 沿径向排列基因名称，长名称可换行。 |
| `pident` 与默认分段标度 | 将比对一致性百分比映射到带的颜色，并显示 默认分段色标。 |

## 阅读与调整 {#read}
先看序列颜色和基因方向，再沿较大的比对区间阅读连接。连接交叉本身不能证明重排。色条使用 ggchord 默认分段标度，按实际显示的比对一致性读数。

需要减少连接时查看[连接筛选](ribbons.html)，调整文字密度时查看[标签](labels.html)，改变半径和曲率时查看[逐序列参数](parameters.html)。从文件读取开始学习可先看[真实序列比较](real-data.html)。

## 按预定尺寸保存 {#save}
```r
ggsave("four-genome-chord.png", p_show_chord, width = 12, height = 10, units = "in", dpi = 144)
```
