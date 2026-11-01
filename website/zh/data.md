## 输入表 {#tables}
| 输入 | 必需字段 | 含义 |
| --- | --- | --- |
| 序列 | `accver`, `length` | 唯一标识；正值的完整序列长度，单位 bp |
| 区间连接 | `qaccver`, `saccver`, `qstart`, `qend`, `sstart`, `send` | query 与 subject 锚点 |
| 基因 | `accver`, `start`, `end`, `strand` | 链为已知 `+`/`-`；`anno` 可选 |
| 通用 feature | `accver`, `start`, `end`, `strand` 或已声明方向 | 可保留未知链；`type`, `label`, `feature_id` 可选 |

标识可以是任意字符串，accession version 原样保留。只有映射、筛选或统计需要时，
才要求 `pident` 和比对 `length`。不能由 feature 的范围推断完整序列长度。

## 坐标契约 {#coordinates}
基因组坐标采用 **一基、闭区间 bp**。BED/PAF 读取器已转换源坐标约定，不要重复转换。
基因起止坐标有序，链单独存储；连接端点的次序可表达方向。没有方向证据的单 bp
连接保持未知。未知链不会被静默改成 `+`；应使用 `geom_feature()`，而非严格的基因层。

氨基酸坐标需要[显式 CDS 投影](projection.html)。笛卡尔布局单位、基因组 bp、角度、
物理 mm 和轨道数值单位是不同的域。详见完整的[数据契约](reference/ggchord_data_contract.html)。

## 导入原始记录 {#import}
在仓库根目录运行这个仅用于源码仓库的例子，不改写原始文件。

{{code:imports}}

`read_genbank()` 保留记录/feature 来源和分段位置。`read_bed()`、`read_paf()` 保留
坐标转换元数据。`as_ggchord_seq()`、`as_ggchord_feature()`、`as_ggchord_link()`
适配支持的 S4 对象；对象缺少完整长度或标识时必须显式提供。
`file` 与 `files` 均支持路径向量或 glob，重复解析路径按首次出现顺序只读一次。
一个工作流选择一个清晰入口。缺失的序列或同源关系应由分析者
提供，绘图层不生成这些业务数据。

## 校验、清理与准备 {#validation}
先运行[共享准备代码](index.html#quick-start)：

{{code:validation}}

校验报告问题；清理执行明确修复并返回报告；准备函数组合工作流步骤，不推断生物学
归属。绘图前检查 `validation`、`cleaned$report` 和 `prepared$report`。
几何表的 `source_row` 指图层已解析输入行，`.source_line` 才是物理文件行；GenBank
另有记录/feature 来源。聚合可带 `source_rows` 和 `.length_estimated`，筛选合并时应保留。

## 检查问题输入及修复报告 {#repairs}
{{code:cleaning}}
{{results:cleaning}}

演示输入有意包含未知 ID、反向区间及超过 0–100 的相似度。明确修复不表示应该静默
截断科学观测。检查 `problem_report`、`repaired$report`，保留原表，并按实际证据
选择策略。合法反向比对可以保持端点顺序，而非排序。`prepare_ggchord_data()`
因此有其工作流相关的默认策略，不能假定独立清理和准备函数完全相同。

## 声明 genome 与 contig 键 {#hierarchy}
{{code:hierarchy}}

`hierarchy$contigs` 汇总已声明归属，不拼接 contig，也不由 accession 前缀推断归属。
不同样本重复局部 contig 名时须显式提供全局唯一标识。
[文件格式与适配器](formats.html)说明格式及可选对象，[真实记录教程](real-data.html)
把这些报告带入完整绘图流程。
