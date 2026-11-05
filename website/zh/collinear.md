<div class="lesson-intro">模块二的综合结果：四条真实基因组、分链特征和向外绕行的跨行连接。</div>

## 综合结果 {#overview}

{{figure:showcasecollinear}}

四条骨架按碱基位置比较，精选的 16 条 CDS 放在各行两侧。先沿一条横行读基因位置，再跟随带状连接看匹配区间。实线轮廓连接相邻行，虚线轮廓提示中间还有其他显示行。若刚开始，先做[四条骨架示例](collinear-start.html)，再完成[添加比对和 CDS](collinear-build.html)。

## 准备原始记录 {#prepare}
解压[教学输入包](../assets/tutorial-inputs.zip)，安装dplyr、purrr、tidyr后运行：

准备脚本包含原始数据读取、筛选和显示设置。将脚本放在当前工作目录。

{{recipe:collinear}}

保留四条完整FASTA、266条CDS及32条原始HSP。相邻行使用明确筛选的观测区间；非相邻每对仅保留最长一条，减轻长连接重叠。

## 为连接和标签保留空间 {#gc}
如需比较序列组成，可在[轨道专题](tracks.html)学习添加 GC 轨道。

## 组合图层 {#compose}
{{code:showcasecollinear}}

虚线外框标识非相邻显示行，较明显的线宽便于辨认。控制点沿右侧分流，局部smooth避障在箭头附近退让。identity左侧色条与弦图一致，底部链图例使用箭头。比例尺水平offset为0，长度仍按共享bp转换计算。

## 调整与限制 {#read}
修改seq_gap、行顺序、方向或连接筛选后，应同步重建控制点。局部箭头避障不是绕过中间整条序列的自动路由器；跨行ribbon之间仍可合法相交。需要局部观察时使用[focus与分段](views.html)。当前显式圆形半径/曲率不能直接传入共线性坐标，见[转换边界](collinear-data.html#contract)。

## 按预定尺寸保存 {#save}
```r
ggsave("four-genome-collinear.png", p_show_collinear, width = 12, height = 8, units = "in", dpi = 144)
```
