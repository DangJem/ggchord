# v0.16.2 数据契约与证据

日期：2026-10-31。开发基线 `4510d3686fcb18ac9e691a37b17c1ac08fc211f6`。
本文件描述当前源码合同；有限运行范围见完成记录和 `tools/validate-data-contracts.R`。
旧审计 CSV 保留历史，不作为新版运行通过证明。用户帮助入口为
`help("ggchord_data_contract", package = "ggchord")`。

## 入口、返回与失败

下表每一行对应一个公开入口；所有输入均保留原对象，不在绘图时执行外部分析。
NULL/NA 的含义按职责区分，不强制制造一个全包统一含义。

| 入口 | 输入与 NULL/NA | 返回、来源与错误边界 |
| --- | --- | --- |
| `read_blast` | 本地 UTF-8；file/files 均可向量或 glob；至少提供一个 | 六 canonical 端点、可选统计列、物理源行、import 属性；字段数异常停止；数值转换失败警告并保留 NA；默认不把 `001` 或字面 ID `NA` 改写 |
| `read_gff3` | 本地 GFF3；feature_types 非空；未知链由 unstranded 显式决定 | 9列记录选择后的 feature 表，ID/Parent 解码、原 attributes/source_strand/物理源行；无目标类型停止，数值解析失败警告；legacy plus 默认保留，推荐 keep 用于未知链 |
| `read_fasta_lengths` | 本地 FASTA；header_delim NULL 取首个空白前 ID | ID、实际非空白字符数、header 物理源行；无 header、空 ID、header 前孤立文本停止；重复 ID 警告，不推测拼接、单位或业务 ID |
| `read_bed` | BED3–12；零宽拒绝；可无 strand | start0+1/end0、原零基字段、segments、物理源行；字段宽度/块越界/重叠异常停止；未知链 `.`/none |
| `read_paf` | 12个必需字段及可选 tags；零宽和越界拒绝 | query 增序，反向 target 降序；原坐标、tags、物理源行；length/pident 来自报告计数，非重新比对；mapq255→NA 且保留 source_mapq |
| `read_genbank` | 本地 UTF-8 GenBank/GBFF；现有 Python/Biopython 与 jsonlite；未知长度需显式 lengths | seq/features/report；记录与 feature 索引、qualifiers、原文本、位置描述、parser版本；unsupported error/report 明确分流；不安装、不联网、不猜长度或亲缘关系 |
| `as_ggchord_seq` | 表/Seqinfo/GRanges/GInteractions/XString；未知全长需 lengths；单位保留 | 普通表与 conversion；已知长度冲突停止；Seqinfo 不从端点推测全长；metadata/seqlevels/名字/来源链保留；非 bp 表不能直接当基因组绘图 |
| `as_ggchord_feature` | 表/GRanges；ID与坐标不可缺失；未给链允许未知 | 普通 feature 表；*→.，source_strand 保留；正整数坐标；metadata 冲突加前缀并记录 column_map；不擅自排序反向端点 |
| `as_ggchord_link` | 表/GInteractions；六端点必需，无统计值可不提供 | 首 anchor=query、次 anchor=subject，负链端点降序；未知链 orientation_known=FALSE，direction=NA；两端 metadata、Seqinfo、object metadata 保留，不补 length/pident |
| `validate_ggchord_data` | constructor 三表；NULL 可选表不检查；非法类型形成结构问题 | validation(errors/warnings/summary/invalid_rows/cleanable)；strict 收集后决定停止；参数逻辑错误立即停止；NA/未知 ID 与越界有源行，校验不改输入 |
| `validate_ggchord_hierarchy` | 显式 genome/contig 表；主外键不可缺失 | 原表、issues、contigs；重复/缺失/未知外键为错误，孤立 genome 为警告；不推断成员或拼接 contig |
| `clean_ggchord_data` | data.frame 三表；NULL 保留缺席；参考 ID 必须有效 | 复制后的三表、逐变更报告；缺失/空/未知 ID 按 drop/error/keep；keep 不是可绘图承诺；反向默认 sort 不变 |
| `prepare_ggchord_data` | 适配器输入与命名 features/links 集合；mapping/steps 为命名 list | 三表、独立集合与报告；先适配/验证，仅按请求 clean→filter→deduplicate→merge→bundle→optimize；clean 默认 keep，方向丢失会阻止后续方向操作；来源行集合逐步组合 |
| `find_ggchord_loci` | 明确 feature_id/targets、完整长度；未知链只允许对称扩展 | 独立 loci、原 target 范围、源行集合、clip 标记；不按注释猜 ID；反链上游向大坐标；跨原点需显式转换 |
| `focus_ggchord_data` | seq、可选连接/feature、显式 loci；legacy finite numeric | 窗口相对三表/loci/report；原 ID、最早源坐标与截断标记；反向 loci 按 min/max 不代表跨原点；部分连接仅做线性显示投影，length 标估计、pident 不重算 |
| `focus_ggchord_track` | point/interval 表、参考和窗口；value 可 NA，坐标不可 NA | 保留数值/metadata、不插值；源坐标与 output-input 行、遗漏行；不跨窗口接线；反向区间要求先分段 |
| `segment_ggchord_data` | 正整数有序窗口、可选三表/点连接、显式 orientation | 独立窗口与 segments/transforms/report；不拼接源序列；split/drop 明确；候选复制数预先限额；遗漏源行有报告 |
| `map_ggchord_position` | collinear 图、显示窗口 ID、有限位置；source/inverse 为逻辑单值 | position/x/y 表；可逆局部或源 bp；分段 gap/越窗拒绝；非 segmented 旧普通线性映射允许外推，不推广为所有输入皆严格拒绝越界 |
| `build_ggchord_cds` | 明确 CDS、phase/codon_start、链与可选 order/partial；异常翻译声明拒绝 | models/segments/source_data/expanded_data/provenance；只跳首 phase，后续验证连续；显式跨原点 order 可用；genetic_code 为保留元数据，不执行翻译 |
| `project_ggchord_features` | 明确 parent 与 protein/genome 区间；跨 codon 由 error/expand | segment 化普通表、源坐标/row、excluded_bases/codon_expanded；不把 intron 包络当 coding；未知 parent、冲突链/坐标或无有效 coding overlap 停止 |
| `project_ggchord_links` | 明确 ungapped=TRUE、两个 CDS 与等长残基区间；未知 parent 拒绝 | 按两侧 CDS 边界拆分的基因组连接，源行/piece/方向；protein length/pident 只留 source_ 元数据，不伪造 nucleotide alignment |
| `filter_ggchord_ribbons` | pair 表；可选统计列按启用过滤器要求；阈值单个有限数 | 完整保留行、source_rows、逐行理由；排序另记来源；统计 NA 在对应过滤器下删除，不把缺值当0；非法阈值/缺列明确停止 |
| `deduplicate_ggchord_ribbons` | 六端点及 keep 所用统计列；必用数值必须有限 | 代表行和 removed→最终代表报告；overlap 两侧均为真正互惠比例，降序可用，保留相反相对方向；不合并统计量 |
| `merge_ggchord_ribbons` | 六端点/length/pident、有效 pair；group_by 是附加分组 | 单行保留原 length；合并包络标 .length_estimated；report 保留源行/length/pident；默认拒绝不合相对方向的 block 次序；extra 冲突默认 typed NA；不是新 alignment |
| `bundle_ggchord_ribbons` | 六端点和显式权重所需统计；参考长度、bins/min_bundle | 按 pair/相对方向/bin 聚合；包络沿首行绝对方向；length 为源 length 加权均值取整，pident 为 length 加权均值，冲突额外列 typed NA；报告源行集合 |
| `optimize_ggchord_layout` | seq/连接/权重、可选 order/orientation；NA 参数拒绝 | 显示排序/翻转建议与 objective/report；原生物链/query-subject 不变；approximate 只影响目标评估，不声称全局最优 |
| `stat_ribbon_bundle` | 标准 layer 与 bundle 参数 | 对 bundle 结果作几何，after_stat 暴露 bundle_n/weight/density；来源聚合含义与 utility 一致；不是绘图时重新比对 |
| `stat_ribbon_density` | 标准 layer、bins/weight/group_by | 每输入连接保留一条，报告与计算变量；pair 内 bin 权重相对最大值，非 coverage/概率；默认/显式固定样式规则沿用 v0.16.1 |
| `link_ggchord_clusters` | 明确 membership/pairs，paralog/segment 策略 | 连接与 members/comparisons/report；未知 feature 报告遗漏，缺失/冲突 ID 停止；未知链无方向证据；展开前限额；同源关系来自用户声明 |
| `primer_catalog` | set/aliases 显式选择，非法逻辑拒绝 | distinct oligo 或全部原 alias 行、数据库 provenance 属性；source direction 只为提示 |
| `find_primer_bindings` | DNA、可选自定义 oligo/alias；相同 ID 不同 DNA 拒绝 | site 表/aliases；exact 或最长精确3′ suffix，mismatch=0；Wallace 为简单估计；reference 模式保留参考注释而非新实验结论 |
| `find_restriction_sites` | DNA、源数据库/自定义 motif、circular/methylation | 识别范围/未折返及显示 cut/来源行/版本与上下文标记；未知 cut 保持 NA；linear 上下文不接末首；仅有限敏感性规则，不模拟完整消化 |
| `filter_restriction_sites` | site 表、preferred/commercial/count 等显式策略 | 返回原表子集，保留 match_id/pattern_source_row 与属性；无独立 removed-row 报告；显示筛选不改变识别序列，不替代酶反应选择或实验验证 |
| `find_common_features` | DNA、库与过滤/mode/阈值 | method/identity/coverage/heuristic confidence、segments 与库来源；标准码六框架不支持非标准翻译、剪接或跨原点蛋白；source geneticCode 不被解释为算法参数 |

## 已作出的取舍

- 单bp端点不带升降序，已有strand/direction可作证据，缺失时保持unknown；
  merge/bundle对此类行不造聚合包络，dedup/density不将unknown当作same。
- F03 保留任务相关默认：独立 clean 的 sort 与 prepare 的 keep 都不改，方向敏感
  工作流显式使用 keep；不是尚待用户统一默认。
- F16 保留 file/files，不为名称统一删除工作入口；两个参数共用同一展开与去重路径。
  现在同一规范化路径仅读一次，file 中 glob 与 files 一致。没有参数删除版本。
- F17 保留现有成熟可选解析器；显式 UTF-8 与失败信息，无联网/安装。unsupported
  政策保留；不是把所有不支持位置变成 guessed envelope。
- F20 保留适配器/未知链/metadata/Seqinfo 合同，真实对象验证与缺依赖路径分别记录。
- F30 不扩展搜索算法；明确其标准遗传码、线性蛋白框架、甲基化目录和启发式置信边界。
  非标准翻译的最小教学输入证实它不在当前能力范围；不虚称支持任意 geneticCode。
- bundle 的相对方向分组保留：同时翻转两端不改变相对关系；首行绝对方向只是包络表示。
  需要分开来源方向时传 `group_by = c("qstrand", "sstrand")` 等实际存在的列。
- 不在本版统一历史 chord feature/axis 与 point link 的显示锚点公式，不修改存储坐标。
  后续内部输出阶段可继续核对，但不得把已有显示差别直接改成生物学坐标迁移。

## 来源与外部规范

坐标定义对照 [UCSC BED 规范](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)、
[PAF 规范](https://github.com/lh3/miniasm/blob/master/PAF.md)、
[GFF3 规范](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)；
位置类型对照 [Biopython SeqFeature](https://biopython.org/docs/latest/api/Bio.SeqFeature.html)，
遗传码边界对照 [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)。
仅引用规范与公开解析 API，未复制第三方实现。

内置数据来自既有 `data-raw/generate_*_database.R`，保持原始工作簿、序列与源记录：
common-feature 的库版本/原文件和工作簿 SHA、primer 的 Addgene URL/源 sheet/row/SHA、
restriction 的 REBASE/EMBOSS version/原 motif/source row 均可追踪。本版未替换、下载或
重新生成这些资产；既有来源记录不等于本版重新认证全部第三方资产的再分发许可。

## 验证与未运行边界

机读参数表使用当前 formals，复用已有逐参数语义并以本版结论覆盖变化处；旧 review
SHA 单独保留，源码散列另存。运行证据是有限的代表输入，不能将表中每个参数都标成
运行通过。真实可选依赖检查覆盖安装在临时库中的类对象；未安装这些依赖时公开入口
仍保持明确错误。没有执行完整 testthat、包构建/check、vignette 或发布设备矩阵。
