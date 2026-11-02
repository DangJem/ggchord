<div class="lesson-intro">模块二的综合结果：四条真实基因组、分链特征和向外绕行的跨行连接。</div>

## 综合结果 {#overview}
{{figure:showcasecollinear}}

共线性来自弦图的坐标表达扩展；复用数据和图层是设计目标。先读[数据](collinear-data.html)、[坐标替换](collinear-build.html)、[路由与比例尺](collinear-routing.html)，再运行下方完整代码。沿用弦图的精选候选与相似度标度；线性图省略两条拥挤的短候选，展示16个CDS，每个特征传入标签层。正负链在骨架上下分列。

## 准备原始记录 {#prepare}
解压[教学输入包](../assets/tutorial-inputs.zip)，安装dplyr、purrr、tidyr后运行：

{{prepare:sequencepalette,docstyle}}
{{code:showcasedata}}

保留四条完整FASTA、266条CDS及32条原始HSP。相邻行使用明确筛选的观测区间；非相邻每对仅保留最长一条，减轻长连接重叠。

## 为连接和标签保留空间 {#gc}
首页不再添加GC轨道，序列组成分析留在[轨道专题](tracks.html)。避免将两种不同的阅读任务挤在同一宣传图中。

## 组合图层 {#compose}
{{code:showcasecollinear}}

虚线外框标识非相邻显示行，较明显的线宽便于辨认。控制点沿右侧分流，局部smooth避障在箭头附近退让。identity左侧色条与弦图一致，底部链图例使用箭头。比例尺水平offset为0，长度仍按共享bp转换计算。

## 调整与限制 {#read}
修改seq_gap、行顺序、方向或连接筛选后，应同步重建控制点。局部箭头避障不是绕过中间整条序列的自动路由器；跨行ribbon之间仍可合法相交。需要局部观察时使用[focus与分段](views.html)。当前显式圆形半径/曲率不能直接传入共线性坐标，见[转换边界](collinear-data.html#contract)。

## 按预定尺寸保存 {#save}
```r
ggsave("four-genome-collinear.png", p_show_collinear, width = 12, height = 8, units = "in", dpi = 144)
```
