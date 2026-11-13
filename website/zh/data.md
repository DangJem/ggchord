<div class="task-intro"><p>把文件中的观测整理成绘图角色。先检查最小字段与坐标，再使用校验报告决定是否清理；不要为了出图猜测缺失的长度或方向。</p></div>

## 本页准备 {#prepare}
先选出校验所需的列。运行文件示例前解压[教学输入](../assets/tutorial-inputs.zip)，使工作目录包含 `examples/`。

{{prepare:setup}}

## 输入表 {#tables}
[基本数据表](chord-data.html#tables)展示序列、连接和特征的字段与实际数据。接下来用这些字段校验导入结果。

## 坐标契约 {#coordinates}
基因组坐标采用 **一基、闭区间 bp**。BED/PAF 读取器已转换源坐标约定，不要重复转换。
基因起止坐标有序，链单独存储；连接端点的次序可表达方向。没有方向证据的单 bp
连接保持未知。通用 feature 保留未知链。兼容例外是 `read_gff3()` 默认
`unstranded = "plus"` 会转为 `+`，但保留 `source_strand`。要保留未知链，读取时
显式使用 `unstranded = "keep"`，绘制时使用 `geom_feature()`。

蛋白骨架使用原生 aa 坐标，见[蛋白比对](protein-alignments.html)。显示在核酸骨架上时需要[CDS 投影](projection.html)。笛卡尔布局单位、基因组 bp、角度、
物理 mm 和轨道数值单位是不同的域。详见完整的[数据格式说明](reference/ggchord_data_contract.html)。

## 导入原始记录 {#import}
在仓库根目录运行这个仅用于源码仓库的例子，不改写原始文件。

{{code:imports}}

{{console:imports}}

`read_genbank()` 保留记录/注释来源和分段位置。`read_bed()`、`read_paf()` 保留
坐标转换元数据。`as_ggchord_seq()`、`as_ggchord_feature()`、`as_ggchord_link()`
适配支持的 S4 对象；对象缺少完整长度或标识时必须显式提供。
`file` 与 `files` 均支持路径向量或 glob，重复解析路径按首次出现顺序只读一次。
一个工作流选择一个清晰入口。缺失的序列或同源关系应由分析者
提供，绘图层不生成这些业务数据。

## 校验、清理与准备 {#validation}


{{code:validation}}

{{console:validation}}

校验报告问题；清理执行明确修复并返回报告；准备函数组合工作流步骤，不推断生物学
归属。绘图前检查 `validation`、`cleaned$report` 和 `prepared$report`。
几何表的 `source_row` 指图层已解析输入行，`.source_line` 才是物理文件行；GenBank
另有记录/注释来源。聚合可带 `source_rows` 和 `.length_estimated`，筛选合并时应保留。

## 检查问题输入及修复报告 {#repairs}
{{code:cleaning}}
{{results:cleaning}}

例子只修改原始记录的副本，故意加入未知 ID 和超过 0–100 的相似度。明确修复不表示应该静默
截断科学观测。检查 `problem_report`、`repaired$report`，保留原表，并按实际证据
选择策略。合法反向比对可以保持端点顺序，而非排序。`prepare_ggchord_data()`
因此有其工作流相关的默认策略，不能假定独立清理和准备函数完全相同。

## 声明 genome 与 contig 键 {#hierarchy}
{{code:hierarchy}}

{{console:hierarchy}}

`hierarchy$contigs` 汇总已声明归属，不拼接 contig，也不由 accession 前缀推断归属。
不同样本重复局部 contig 名时须显式提供全局唯一标识。
[文件格式与适配器](formats.html)说明格式及可选对象，[真实记录教程](real-data.html)
把这些报告带入完整绘图流程。

## 最小输入与蛋白比对 {#minimal}
{{results:minimal}}

构造函数的第二个位置参数是 `ribbon_data`，不是 `gene_data`。建议像共享例子一样
写明参数名，避免把基因表误当成连接表。相似度、比对长度和覆盖率字段只在相应
映射、筛选或统计中才需要。详见[标准 BLAST 与蛋白输入](alignments.html)和
[按序列、按链设置参数](parameters.html)。

## 完成后检查 {#check-result}
输入 ID 必须能匹配到完整序列表；原始记录与修复副本分开保存。反向比对端点递减可以合法，不能把所有警告都当作需要删除的数据。

<div class="next-step">继续：<a href="formats.html">导入和转换数据 →</a></div>
