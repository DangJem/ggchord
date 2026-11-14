<div class="lesson-intro">模块二的综合结果：四条真实基因组、分链特征和向外绕行的跨行连接。</div>


## 研究问题、对齐与取舍 {#question}

这张图回答：在同一 bp 尺度下，四条完整序列的已选择 CDS 与比对区间如何相对排列？源表仍是四份 FASTA、GFF3 CDS 和原始 BLAST HSP，与[弦图案例](chord.html#question)相同。自动方向是按连接支持选择的显示方向；center 对齐整条序列，不声称各序列中心同源。

相邻行展示经过明确阈值筛选的记录；跨行每对只保留最长 HSP。省略两个拥挤特征是图中取舍，完整 source 表保留。右侧控制点由当前布局推导，所以改变顺序、方向或行距后必须重新计算。局部避障让开箭头，并不保证整条长路径绕过所有中间行。

## 综合结果 {#overview}

{{figure:showcasecollinear}}

四条骨架按碱基位置比较，精选的 16 条 CDS 放在各行两侧。先沿一条横行读基因位置，再跟随带状连接看匹配区间。实线轮廓连接相邻行，虚线轮廓提示中间还有其他显示行。若刚开始，先做[四条骨架示例](multi-sequence-start.html#collinear)，再完成[添加比对和 CDS](multi-sequence-start.html#collinear)。

## 准备原始记录 {#prepare}
解压[说明输入包](../assets/tutorial-inputs.zip)，安装dplyr、purrr、tidyr后运行：

准备脚本包含原始数据读取、筛选和显示设置。将脚本放在当前工作目录。

{{recipe:collinear}}

保留四条完整FASTA、266条CDS及32条原始HSP。相邻行使用明确筛选的观测区间；非相邻每对仅保留最长一条，减轻长连接重叠。

## 为连接和标签保留空间 {#gc}
如需比较序列组成，可在[轨道专题](tracks.html)了解添加 GC 轨道。

## 组合图层 {#compose}
{{code:showcasecollinear}}

虚线外框标识非相邻显示行，较明显的线宽便于辨认。控制点沿右侧分流，局部smooth避障在箭头附近退让。identity左侧色条与弦图一致，底部链图例使用箭头。比例尺水平offset为0，长度仍按共享bp转换计算。

## 调整与限制 {#read}
修改seq_gap、行顺序、方向或连接筛选后，应同步重建控制点。局部箭头避障不是绕过中间整条序列的自动路由器；跨行ribbon之间仍可合法相交。需要局部观察时使用[focus与分段](views.html)。当前显式圆形半径/曲率不能直接传入共线性坐标，见[转换边界](multi-sequence-views.html#input)。

## 按预定尺寸保存 {#save}
{{code:courseblockf}}

## 保存前的证据检查 {#evidence-check}

按图中已声明的宽高调用 view_ggchord，再用相同尺寸 ggsave 保存。核查 ID、源链、比对端点方向、筛选报告、图例单位及标签隐藏原因；PNG 能渲染不代表字体在所有 PDF/SVG 阅读器中一致。保留完整源表、显示子集、参数、软件版本和输入哈希，使图可以追溯到记录。数据来源说明可下载自[资产与来源记录](../assets/ASSET-SOURCES.md)。

## 来源权利与复现 {#source-rights}

包代码采用 MIT 许可；原始数据库记录及参考目录保留各自使用条件，包许可不重新授权第三方数据。导出时应保留记录身份和来源。[资产与来源记录](../assets/ASSET-SOURCES.md)、输入归档哈希及本页准备脚本共同描述本例，计算识别结果仍需实验核实。
