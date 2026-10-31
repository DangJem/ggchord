# v0.16.2 本地开发完成记录

日期：2026-10-31。基线 `4510d3686fcb18ac9e691a37b17c1ac08fc211f6`；实现提交
`19000f1`。状态：已完成本版全部核查、实施与相关轻量验证；按授权流程提交并整合
本地 main。未执行 GitHub 正式发布验收，未 push、tag 或发布，也未启动下一版。

## 实际结果

- 34个公开数据入口、195个当前参数，逐项记录输入、NULL/NA、返回、方向、来源和
  失败合同；保留旧审阅出处，不把旧快照当作当前接口。见
  [完整契约](API-DATA-CONTRACTS-v0.16.2.md)及 `audit-generated/data-*-v0.16.2.csv`。
- F03 保留 clean sort / prepare keep 的任务相关默认；F16 保留 file/files 双入口并
  统一展开、顺序与规范化路径去重；F17/F20 保留可选成熟解析器与对象适配职责。
  没有删除、改名或新增公开函数参数，没有待用户选择的本版兼容冲突。
- F42 修正真正互惠重叠与反向端点；去重保留不同相对方向，validator不把圆图跨原点
  bounds猜成线性重叠。修正缺失ID和gene重复的清洗建议。
- F43 保留未合并行alignment length；pair始终独立，默认不合并关系顺序交叉区块；
  合并包络标估计，报告保存源length/pident，不把包络统计冒称新alignment。
- F44 linear甲基化上下文不再连接末首；F45 primer标签及metadata使用factor数值标签。
- F46 明确UTF-8、字段数/路径/参数错误，默认保留001和字面NA等ID；新增物理源行；
  GenBank显式UTF-8打开；不下载、不安装、不改原始业务数据。非法filter/protein阈值和
  fractional primer annealing length明确拒绝，正常单值调用不变。
- F47 单bp端点不再否定已有reverse或误作same；clean保留裁剪前证据，merge/bundle
  让单bp端点行独立，dedup/density区分已有方向与unknown，缺失方向不被推断。
- 聚合、投影、cluster、primer、restriction和common-feature的来源、信息损失与
  置信边界已写入新帮助页、roxygen/Rd、NEWS、迁移、决策和版本文档。

## 实际验证

环境：R 4.6.0、ggplot2 4.0.3、roxygen2 8.0.0。真实可选类来自已存在的
`/private/tmp/ggchord-bioc-lib`：GenomicRanges 1.64.0、InteractionSet 1.40.0、
Biostrings 2.80.2。离线GenBank使用现有Python 3.12和Biopython 1.86源码目录，
通过PYTHONPATH引用现有NumPy缓存；没有安装或修改用户环境。

| 验证 | 已观察结果与范围 |
| --- | --- |
| `Rscript tools/validate-data-contracts.R` | PASS：包含/降序/不同方向去重，单行与合并来源，跨pair/交叉block，单bp clean→merge→bundle，density staged count与unknown分组，F31 keep，file/files/glob/ID/源行/畸形字段/UTF-8/参数错误，BED单bp/零宽，linear/circular上下文，primer factor；标准码能力边界 |
| `Rscript tools/validate-clean-intervals.R` | PASS：单bp inclusive和clip保留/移除边界，F01/F02/F31另由定向脚本和既有公开行为覆盖 |
| `Rscript tools/validate-data-formats.R` | PASS：BED/PAF坐标、segments、tags、未知链、mapq来源、零宽/异常块错误 |
| `R_LIBS=/private/tmp/ggchord-bioc-lib Rscript tools/validate-data-convert.R` | PASS：真实GRanges/GInteractions/XStringSet/XString、metadata冲突、名字/单位/未知链/端点顺序和绘图；未将缺依赖路径写成类验证通过 |
| Seqinfo直接脚本 | PASS：未使用seqlevel仍保留在conversion属性，isCircular/genome/source metadata保留，冲突长度明确报错 |
| `Rscript tools/validate-genbank.R`（上述现有解析环境） | PASS：离线多record、joined/反链/partial/源文、codon_start/declared length、origin join、remote report、含汉字与空格路径；该脚本默认R库下Biostrings检查跳过，但另一个真实类脚本已通过 |
| GenBank失败直接脚本 | PASS：缺解释器与选中Python缺Biopython给具体失败信息；未改环境、未尝试联网 |
| `Rscript tools/validate-data-prepare.R` | PASS：显式步骤、原输入/命名集合、源行组合、方向保护、GFF Parent与源strand；代表性collinear build |
| `Rscript tools/validate-data-project.R` | PASS：CDS phase、正反链、protein/genome往返、分段与ungapped links、partial/异常输入边界；不执行翻译 |
| `Rscript tools/validate-data-relations.R` | PASS：hierarchy、cluster显式成员、paralog、未知链、分段anchor、展开限额与来源 |
| `Rscript tools/validate-focus-coordinates.R` | PASS：nested segments/cleavage、最早源坐标、反向连接、linear length估计与重复focus |
| `Rscript tools/validate-segments.R` | PASS：独立窗口、cross_gap、局部翻转、点来源、map逆向与gap拒绝；代表性build |
| 基线函数直接对照 | F42包含误报/降序漏检、F43交叉block从1行→2行、F44 linear blocked→unblocked、F45 factor 1/2→100/300、F46 ID 1/NA→001/字面NA、F47旧prepare拒绝→保留reverse均有运行证据 |
| 文档/静态记录 | 当前R源码与Rd无警告解析；只执行roxygen Rd/namespace生成；`audit-data-contracts.R`保存34入口/195参数/返回表及源码MD5；diff无空白错误 |

以上是各脚本的代表输入，不是每个参数、任意尺寸或等价语法的穷举。新增默认testthat
只保留严重缺陷的必要回归，导入成功案例复用原案例；本轮没有执行testthat。

## 已明确的能力边界

- common-feature protein模式只作标准码、线性六框架搜索；source geneticCode不是
  可选翻译表。非标准码、剪接或跨原点蛋白需DNA模式或外部明确分析；无命中不证明
  生物学不存在。confidence不是概率或实验验证，primer/restriction同样只是搜索与注释。
- legacy validate/clean/raw focus接收有限numeric坐标；严格正整数由as适配器或特定
  工作流保证。单bp是合法区间；零宽半开插入和不安全位置有明确拒绝/报告路径。
- focus/merge/bundle的统计保留或估计含义不同；使用原始alignment进行科研统计。
  query/subject、原ID、未知链与display orientation不被混淆。
- 保留历史chord锚点的连续显示公式差别，没有在数据版本改变默认几何。
- 正式发布仍须用户明确授权后运行完整testthat、R CMD check、英语/汉语vignette与
  PNG/PDF/SVG验收；本轮没有运行这些流程，也没有重建冻结docs网站。引擎状态与稳定
  导出字段继续归v0.16.3，网站/教程全量重构归v0.16.5，平台收敛归v0.16.6。

本版已选范围没有未完成开发项；本地完成与正式发布验收、已发布三种状态严格分开。
