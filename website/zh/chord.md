<div class="lesson-intro">把四条基因组放在同一张关系图中，组合直线与弧形骨架、保留真实链方向的基因、连续相似度标度和精选标签，同时保留数据来源。</div>

## 综合结果 {#overview}
{{figure:showcasechord}}

本页是模块一的终点。初学请按[输入表](chord-data.html) → [生成与导入](chord-import.html) → [逐层绘图](chord-build.html) → [参数与避障](chord-tuning.html)学习。四条真实序列上只展示18个精选CDS，每个均交给repel标签层，避免密集箭头淹没连接。

## 准备原始记录 {#prepare}
将[教程输入包](../assets/tutorial-inputs.zip)解压到工作目录，使其包含 `examples/`。使用 `install.packages(c("dplyr", "purrr", "tidyr"))` 安装示例可选依赖；完整流程用它们整理表格和计算 GC 轨道。

{{prepare:sequencepalette,docstyle}}
{{code:showcasedata}}

输入包含四条未改写的 FASTA、对应 GFF3 注释和 32 条原始 HSP。展示保留长度至少 900 bp、相似度至少 75% 的 HSP，每对序列最多三条，以减少拥挤；没有合并或编造比对。

## 组合图层 {#compose}
{{code:showcasechord}}

`geom_seq(seq_curvature = ...)` 通过命名向量将曲率绑定到稳定的序列标识，`0` 表示直线。曲率只改变展示几何，不改变序列长度、生物学链方向或比对端点。

## 阅读与调整 {#read}
先看序列颜色和基因方向，再沿较大的比对区间阅读连接。连接交叉本身不能证明重排。相似度色条与[共线性视图](collinear.html)使用相同范围和配色。

需要减少连接时查看[连接筛选](ribbons.html)，调整文字密度时查看[标签](labels.html)，改变半径和曲率时查看[逐序列参数](parameters.html)。从文件读取开始学习可先看[真实序列比较](real-data.html)。

## 按预定尺寸保存 {#save}
```r
ggsave("four-genome-chord.png", p_show_chord, width = 12, height = 10, units = "in", dpi = 144)
```
