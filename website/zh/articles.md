<p class="lesson-intro">从一张能读懂的图开始，再根据自己的数据和问题选择下一步。</p>

## 第一次使用：按顺序完成 {#start}
1. [绘制第一张图](installation.html)：安装、最小输入、绘图、读图与保存。
2. [选择图形](modes.html)：判断应该看总体连接、单序列注释还是局部位置关系。
3. [理解数据与坐标](data.html)：确认完整长度、一基闭区间与正负链。

## 完成一个真实任务 {#workflows}
| 目标 | 数据准备 | 完成后你会得到 |
| --- | --- | --- |
| [比较四条序列](real-data.html) | 教学输入中的 FASTA、GFF3、标准 BLAST 表格 | 有筛选依据的弦图与全长概览 |
| [绘制完整质粒图](plasmids.html) | 随包 pUC19 DNA | 分层注释、引物与选定酶位点的圆图 |
| [比较对应 CDS 邻域](homology.html) | 随包两条质粒及实际编码 DNA | 有证据的分组、局部窗口和源坐标图 |
| [计算并添加 GC 轨道](tracks.html) | 随包 pUC19 DNA | 有单位、窗口及值域的柱、线和类别带 |

每页给出本页的准备代码。使用原始文件的页面须先[下载教学输入](../assets/tutorial-inputs.zip)，解压后保持 `examples/` 路径；只用随包数据的页面无需下载。

## 已有图，解决一个具体问题 {#tasks}
| 遇到的问题 | 去哪里 |
| --- | --- |
| 文件列名或坐标约定不同 | [导入和转换](formats.html)、[核酸与蛋白比对](alignments.html) |
| 注释重叠或文字拥挤 | [注释堆叠](annotations.html)、[标签布局](labels.html) |
| 只查看局部、分段或改变显示方向 | [局部视图](views.html)、[逐序列参数](parameters.html) |
| 连接太多，想筛选或汇总 | [整理连接](ribbons.html) |
| 配色、图例或层级不清 | [样式映射](styling.html)、[常用调整](adjustments.html) |
| 蛋白范围需要画到基因组上 | [CDS 投影](projection.html) |
| 要保存或解释诊断 | [输出与预览](output.html)、[常见问题](faq.html) |

## 查阅概念与完整参数 {#reference}
[函数文档](reference.html)用于查精确契约；[术语说明](glossary.html)用于理解概念；[图集](gallery.html)用于按结果找教程。使用旧 API 时阅读[迁移指南](migration.html)与[更新日志](news.html)。

## 功能索引 {#coverage}
{{coverage}}

这是查阅索引，不要求按表格顺序学习。精确返回值和限制以各函数参考为准。
