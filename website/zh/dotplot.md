## 两条序列轴上的 HSP 关系 {#task}

综合矩阵保留四条完整序列之间的 32 条原始 HSP。两轴分别保留 accession 身份和局部位置，颜色表达 Identity，线宽表达原记录中的 HSP 长度。

{{manualfigure:flagdotplot}}

## 本模块的主要功能 {#path}

- 多条 query 与 subject 序列的拼接轴、显式顺序和原生单位间隙。
- 正向与反向 HSP 端点，两轴独立的显示方向。
- Identity 等已有指标的 link 映射与标度。
- 精确来源坐标转换、可见窗口和输入行追溯。

[数据、双轴与映射](dotplot-start.html) · [比对结果导入](alignments.html) · [蛋白坐标](protein-alignments.html)

## 输入要求 {#inputs}

提供完整长度、query/subject 标识和每条 HSP 的 qstart/qend/sstart/send。拼接轴不会合并生物学序列身份；bp 与 aa 需要独立坐标域或显式投影。

## 能力边界 {#support}

每个标记是连接 HSP 端点的线段，不是逐碱基重建。错配和缺口需要真实对齐字符或精细轨迹。骨架 gene/feature/track 图层不适用；Dotplot 目前不是 Fusion 的注册端点。
