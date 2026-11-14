<div class="task-intro"><p>把原始 FASTA、GFF3 和标准 BLAST 表格组合成一张可追溯的序列比较图。你将检查输入、选择显示子集，再对照弦图与全长概览。</p></div>


## 问题、数据来源与许可 {#research-question}

本案例以原始文件为起点，回答四条源基因组的已观测相似区间和 CDS 如何组合成全长比较。FASTA/GFF3 标识及来源对应[四序列案例](chord.html#question)。BLAST 标准表由保留原值的原始输出转换，转换规则和输入行来源在 examples/blastn-standard/provenance.csv 中。

包代码使用 MIT 许可；原始数据库记录及第三方参考目录的权利和使用条件按各来源保留，不把包许可自动套到所有输入。说明输入包提供未改写文件、SHA256 和来源说明。格式专题的 EMALE 记录另保留 gggenomes MIT 来源许可，见输入包中的 GGGENOMES-LICENSE.txt。

## 本页准备 {#prepare}
解压[说明输入](../assets/tutorial-inputs.zip)，使工作目录包含 `examples/`。

## 读取并检查原始记录 {#read}
{{code:realimport}}
{{console:realimport}}

FASTA 提供序列 ID 和长度，GFF3 提供注释范围与链，BLAST 提供已计算的比对。验证摘要为空，表示这批输入未报告问题。

## 明确显示筛选规则 {#selection}
{{code:realprepare}}
{{console:realprepare}}

保留一致性至少 80%、比对长度至少 500 bp 的连接。弦图每条序列最多画六条较长的 CDS，以便看清连接。

{{code:realfeatures}}
{{output:realfeatures}}
{{results:real}}

## 构建弦图比较 {#chord}
{{code:realchord}}
{{figure:realchord}}

先读序列名称，再读 CDS 和区间连接。区间连接表达已有比对区间
和映射的相似度，连接交叉本身不构成重排证据。选定注释与更宽的比对区间分开阅读。
科学解释涉及显示子集之外的注释时，应使用完整输入表。

## 在线性视图检查全部注释 {#linear}
{{code:reallinear}}
{{figure:reallinear}}

这张图绘制全部导入注释。区间连接默认只连接相邻显示序列，保留所有比较需要显式
图层设置。bp 参照条属于共线性坐标。裁剪、源坐标轴与独立窗口见[局部视图](views.html)，
按图集明确尺寸预览保存见[输出](output.html)。

## 随图保留证据 {#provenance}
说明输入包包含 FASTA、GFF3、BLAST 表格、SHA256 清单和共享示例脚本。
绘图清单记录代码 hash、PNG 设备、宽高和 dpi。accession 名称标识来源记录，
图片不替代原始注释或比对证据。

## 先选列，再传入函数 {#minimal-fields}
骨架使用 `accver` 与 `length`，连接使用两端 ID 和端点。这个例子还用 `pident` 和 `length` 筛选连接。参见[比对输入说明](alignments.html)。

## 检查要点 {#check-result}
检查序列 ID 在三张表中一致，查看筛选报告，再确认图例表达相似度与链方向。全长概览用于定位；若要研究具体 CDS，请继续聚焦邻域。

<div class="next-step">继续：<a href="modes.html">绘图基本概念 →</a></div>

## 保存本页结果 {#save-result}
在工作目录中保存刚刚创建的图；预览时使用相同宽高。

{{code:courseblockc}}

## 保存前的证据检查 {#evidence-check}

按图中已声明的宽高调用 view_ggchord，再用相同尺寸 ggsave 保存。核查 ID、源链、比对端点方向、筛选报告、图例单位及标签隐藏原因；PNG 能渲染不代表字体在所有 PDF/SVG 阅读器中一致。保留完整源表、显示子集、参数、软件版本和输入哈希，使图可以追溯到记录。数据来源说明可下载自[资产与来源记录](../assets/ASSET-SOURCES.md)。
