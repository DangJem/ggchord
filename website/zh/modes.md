<div class="task-intro"><p>为你的问题选择布局，而不是为了圆形效果改变问题。下面用同一批真实随包记录比较三种表达。</p></div>

## 怎么选择？ {#choose}
| 场景 | 更合适的图 | 应避免的误读 |
| --- | --- | --- |
| 几条序列之间的总体对应关系 | 弦图 | 连线交叉本身不等于倒位 |
| 一条质粒或环形序列上的注释 | 圆图 | 多条序列不能塞进单序列坐标 |
| 比较位置、基因顺序与局部邻域 | 共线性图 | 显示对齐或翻转不等于生物学结论 |

三种图都消费已有数据，不计算新的 BLAST 比对。先尝试下面的最小组合，再进入[原始记录流程](real-data.html)、[完整质粒图](plasmids.html)或[局部同源比较](homology.html)。

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。

{{prepare:setup,sequencepalette,featurepalette,docstyle,plasmiddata}}

下列对象都可按普通 ggplot 组合，
打印或使用 `ggsave()` 保存。

## 弦图：比较多条序列 {#chord}
{{code:chord}}
{{figure:chord}}

连接按 query 序列标识配色，与骨架使用同一类别色板；这不表示相似度或同源基因分组。按真实相似度配色见[原始记录教程](real-data.html#chord)。
`coord_chord()` 控制全局旋转和拟合。逐序列半径、曲率、间隙和方向属于 `geom_seq()`。
连接控制点为 NULL 时，弦图保留历史的原点路径，并不表示各模式统一使用直线。
区间连接默认避障为 `uniform`，点连接默认避障为 `none`。

## 共线性图：对齐局部比较 {#collinear}
{{code:collinear}}
{{figure:collinear}}
本例分别选每条序列最长的已知基因作显示锚点，不将它们当作同源基因。
`anchor` 对齐指定区间，`orientation = "auto"` 按现有证据选择支持的显示方向。
应查看方向诊断，不能直接把显示翻转解释为生物学倒位。区间连接默认连接相邻显示行；
使用 `adjacent_only = FALSE` 或明确 `pairs` 保留其他连接。非相邻边框默认虚线，
显式样式优先。这些规则与弦图不同。[局部视图与分段](views.html)为裁剪或局部翻转
保留源坐标。

## 圆图：注释单条序列 {#circular}
{{code:circular}}
{{figure:circular}}
本例只选 CDS 与复制起点显示，保留完整序列的 bp 刻度。
`coord_circular()` 要求单条序列。通用注释可带未知链。
`geom_feature_plasmid()` 是有独立默认值的便利预设，应在需要这些默认值时显式选择。
引物和酶切图层见[识别注释](recognition.html)。

## 完成后检查 {#check-result}
先看序列数：单序列注释用圆图；多序列总体连接用弦图；位置、距离或邻域比较用共线性图。

<div class="next-step">继续：<a href="installation.html">开始使用：绘制第一张图 →</a></div>
