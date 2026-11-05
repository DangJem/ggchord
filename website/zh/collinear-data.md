输入沿用[基本数据表](chord-data.html)。这里集中说明共线性图特有的方向和行关系。

## 从三张表到横行 {#tables}

序列表的一行对应一条骨架，连接表的一行对应两个区间。注释仍用原始 `strand`。用 Position 把正负链放到上下两侧时，保留原始链；局部窗口的显示间隙也不代表缺失 DNA 的长度。

## 对照一条实际比对 {#example}

{{code:frontlinktable}}

{{output:frontlinktable}}

首行的 query 从 26,298 到 51,139 bp，subject 从 7,121 到 31,959 bp。两端都递增。到[共线性图](collinear-build.html#links)沿第一、二行找这两个区间。换成递减的 subject 时，应保留递减端点，显示翻转由绘图设置单独决定。

## 数据如何获得 {#import}

用 FASTA 获取完整长度，GFF3 获取基因区间，BLAST 或 PAF 获取两个端点区间。ID 必须逐字匹配；基因区间采用一基闭区间，反向比对可以用递减的端点表达。拿到自己的文件后，直接跟做[从源文件到共线性图](real-data.html#read)，再到[线性绘图步骤](real-data.html#linear)。其他文件用[格式入口](formats.html)，蛋白坐标则先做[CDS 投影](projection.html)。

## 当前转换边界 {#contract}

骨架的显示方向可以手动设定，也可用坐标系按比对证据自动选择。每次选择其中一种方式。显示翻转保留原始链方向。相邻指图上行的顺序；若已有远隔行的比对，要明确启用非相邻连接。

已经学过其他布局时，可把通用图层复用于共线性图。逐序列圆弧半径、曲率、环与径向序列标签需在各自的坐标下配置。

继续查看[跨行连接](collinear-routing.html)。
