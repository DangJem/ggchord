<div class="task-intro"><p>把原始 FASTA、GFF3 和标准 BLAST 表格组合成一张可追溯的序列比较图。你将检查输入、选择显示子集，再对照弦图与全长概览。</p></div>

## 本页准备 {#prepare}
在新的 R 会话中先运行下方代码，再按页面顺序继续。配色和排版是可选的显示设置，输入只选择需要的列。 涉及文件时，先解压[教学输入](../assets/tutorial-inputs.zip)，使工作目录包含 `examples/`。

{{prepare:docstyle}}

这条完整工作流使用源码仓库已有的四条序列记录，不需要重新比对或预测基因。
使用匹配的开发 checkout，或[下载教学输入](../assets/tutorial-inputs.zip)，解压后让
工作目录下包含 `examples/`。保存图形时同时记录包版本与[构建清单](../build-manifest.json)。

## 1. 读取并检查原始记录 {#read}
{{code:realimport}}

FASTA 头部提供原样 accession version 与长度，GFF3 提供 CDS 区间和真实链方向，
BLAST 文件提供已有比对；读取器不运行 BLAST。继续前检查 `real_validation`：
连接端点递减可能是合法的反向命中。不能把未知注释链改成 `+` 来满足基因快捷层。

## 2. 明确显示筛选规则 {#selection}
{{code:realprepare}}
{{results:real}}

显示连接保留相似度 ≥80%、已报告比对长度 ≥500 bp 的记录。弦图每条序列选择最长
六条已导入 CDS，并在这些 CDS 中最多标记两条非 hypothetical protein 注释。这是显示规则，不是生物学结论，未显示的
基因不表示不存在。原始 `feature_real`、`link_real` 不变，`real_prepared$report`
和 `real_filtered$report` 保留准备与筛选证据。链方向不重编码；这与生成说明中
明确进行显示链重编码的小型随包教学 fixture 不同。

## 3. 构建弦图比较 {#chord}
{{code:realchord}}
{{figure:realchord}}

先读序列名称，再读真实链方向的注释箭头及区间连接。区间连接表达已有比对区间
和映射的相似度，连接交叉本身不构成重排证据。选定注释与更宽的比对区间分开阅读。
科学解释涉及显示子集之外的注释时，应使用完整输入表。

## 4. 在线性视图检查全部注释 {#linear}
{{code:reallinear}}
{{figure:reallinear}}

这张图绘制全部导入注释。区间连接默认只连接相邻显示序列，保留所有比较需要显式
图层设置。bp 参照条属于共线性坐标。裁剪、源坐标轴与独立窗口见[局部视图](views.html)，
按图集明确尺寸预览保存见[输出](output.html)。

## 5. 随图保留证据 {#provenance}
教学输入包包含未改写的 FASTA/GFF3/BLAST、SHA256 清单和共享示例脚本。
绘图清单记录代码 hash、PNG 设备、宽高和 dpi。accession 名称标识来源记录，
图片不替代原始注释或比对证据。

## 先选列，再传入函数 {#minimal-fields}
`*_source` 对象保留原字段与来源；`seq_real` 只有两列，`feature_real` 只保留几何、
类型和注释文字，`link_real` 因为本例实际筛选和映射相似度，才额外保留 `pident`
和已报告的 `length`。这里使用标准 12 列转换文件，不要求历史 17 列布局。
随包基因示例将真实链保存在 `source_strand`；简短例子恢复该方向后，只选四列传入
基因图层。格式和单位的区别见[比对输入教程](alignments.html)。

## 完成后检查 {#check-result}
检查序列 ID 在三张表中一致，查看筛选报告，再确认图例表达相似度与链方向。全长概览用于定位；若要研究具体 CDS，请继续聚焦邻域。

<div class="next-step">继续：<a href="homology.html">连接同源基因 →</a></div>

## 保存本页结果 {#save-result}
在工作目录中保存刚刚创建的图；预览时使用相同宽高。

```r
ggsave("phage-comparison.png", p_real_chord, width = 10, height = 8, units = "in", dpi = 144)
```
