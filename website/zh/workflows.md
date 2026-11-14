## 按输入与目标选择 {#tools}

先确定手中的原始输入和要得到的图。每条路径都保留来源 ID、完整长度和坐标单位；读取函数不运行比对，也不生成缺失的关系或注释。

| 手中的输入 | 最短处理路径 | 产物与绘图入口 |
| --- | --- | --- |
| 已有 R 表 | 按[基本字段](data.html#tables)整理序列表，核对连接两端和注释 ID；用[校验报告](data.html#validation)定位问题 | 序列长度、HSP、注释三类表 → [共用数据与图层](multi-sequence-start.html) |
| FASTA、GFF3 与 BLAST/PAF | [读取序列及注释](import-sequences.html)，再[读取比对](alignments.html)；核对全长、单位和 ID，保留反向端点 | 真实来源表 → [序列关系图](multi-sequence-start.html)或 [HSP 点图](dotplot-start.html#dotplot) |
| GenBank 单序列 | 在[序列与注释导入](import-sequences.html)使用 `read_genbank()`；查看返回的长度表、注释层级和原始分段 | 单序列及注释 → [圆形骨架与已有注释](circular-start.html) |
| 蛋白 BLASTp 与全长 | [声明 aa 单位并读取比对](protein-alignments.html)，核对两端蛋白 ID；需要核酸落点时再提供真实 CDS/phase | 蛋白坐标 → [蛋白关系图](protein-alignments.html)；核酸落点 → [坐标投影](projection.html) |
| BED、BEDPE 等区间文件 | [读取区间或配对区间](import-intervals.html)，核对读取器完成的一基转换，避免重复换算 | 注释区间或成对链接 → [特征](annotations.html)与[连接](multi-sequence-links.html) |
| bedGraph/BigWig 等数值信号 | [读取信号](import-signals.html)，保留数值单位、缺失及来源区间；不把未知区间补零 | 连续信号 → [数据轨道](tracks.html)；计数与分箱见[轨道统计](tracks-statistics.html) |
| 已有字符对齐 | [读取对齐字符](aligned-sequences.html)，检查 gap、来源坐标和字母类型 | 真实逐列字符 → [局部放大与字符细节](fusion-details.html) |

## 绘图前的检查 {#checks}

1. 序列表的 `accver` 唯一且长度来自完整来源；不能用最后一条注释的终点猜全长。
2. 连接的 `qaccver`/`saccver` 和注释的 `accver` 都能匹配到序列表。不同样本同名 contig 须显式区分，读取器不猜归属。
3. 核酸使用 bp，蛋白使用 aa；BED/PAF 已由读取器转换。连接端点的升降序可能携带方向，不能统一排序。
4. 先画最小图，再按需要加入图层、[局部配置](parameters.html)和[目标尺寸输出](output.html)。问题按[症状排查](diagnostics.html#symptoms)。

## 完整来源实例 {#examples}

[从原始 FASTA/GFF3/BLAST 到序列比较图](real-data.html)保留文件路径、显示筛选、输入与结果。原始文件可从[输入归档](formats.html)下载；下载脚本在解压后、以包含 `examples/` 的目录为工作目录运行。
