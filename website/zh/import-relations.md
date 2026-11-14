{{code:courseblockaz}}

{{prepare:setup}}

先解压[真实输入包](../assets/tutorial-inputs.zip)，工作目录须含 examples/；每节单列可选依赖。

## 变异记录与结构断点 {#variants}

VCF 需要可选 **VariantAnnotation**，显式结构变异解码另外需要 **StructuralVariantAnnotation**。示例读取依赖包未改动的有来源输入：chr22 的 1000 Genomes 子集，以及 StructuralVariantAnnotation 分发的真实 caller VCF。数据保留在这些包中，ggchord 不复制第三方解析器实现或输入。`check-variant-import.R` 记录确切依赖版本和源哈希。

{{code:courseblockba}}

每条可放置记录保留一个 feature 行及 ALT 列表，多等位调用也不拆分。Feature 表示包括 padding 碱基的参考等位基因跨度，不是推断的符号 SV footprint。完整 INFO、样本/基因型数组及表头在 `report$records[[1]]$parsed` 中，应在那里检查实际 caller 字段；缺失支持度或 Identity 不变成零，也不制造计数。未知完整长度保持 NA，绘图前应提供可靠参考长度。端粒和越界跨度保留在 `unplaced_features` 及源对象中，不裁切。

结构解码需显式启用，保留置信区间，关闭缺失 mate 推断；每个互惠对生成一条连接，包含两个原始记录索引。警告标明解码器未配对的记录，报告保留全部记录。`.plottable = FALSE` 标记超出已知参考边界的断点，相关配对保留在 `unplaced_links`。精确单碱基断点可显式选取位置后用于 `geom_link_line()`；不确定区间适合区间几何。`q_breakend_orientation`、`s_breakend_orientation` 描述衍生序列，不是相对同源链。不能制造比对 Identity 为这些连接配色。

现代 BCF（`.bcf` 或 `.bcf.gz`）通过 `python` 指定已安装可选 **pysam** 的现有解释器，例如 `read_vcf("calls.bcf", python = "/path/to/environment/bin/python")`。仅 BCF 会调用它，经 HTSlib 转临时 VCF，并记录解码器、版本和消息。虚拟环境解释器路径保留其环境。原记录顺序、数值及基因型保留；转码文本不是原始 BCF 物理行。不会安装任何软件。本读取器将选定本地文件加载到内存，大数据应在上游选取或子集化。

## 基因组对齐块：MAF 与 chain {#genome-alignments}

使用已有 Python 环境，其中须安装支持 MAF/chain 的可选 **Biopython Bio.Align**，已直接验证 1.85。章节可通过 `ALIGNMENT_PYTHON` 选择该解释器。源码 checkout/输入包须包含 `examples/genome-alignments/`；两个未改动的固定上游版本文件有哈希、许可归属及验证脚本，不运行 aligner。可选解析器不可用意味着跳过示例，不是成功导入。

{{code:courseblockbb}}

MAF 保留实际带 gap 字符和声明的完整源长度，负链分量映射至递减的**正向源坐标**。第一分量锚定关系、源块/文件及原始坐标矩阵均显式保留，不自动扩展成所有可能配对。Gap 单元的源位置缺失，完整长度不能独立验证字符。对齐列在每块重新开始，因此 `map_alignment_position()` 或字符面板前必须筛选一个文件/块。原质量、上下文及空分量证据保留在源与报告中，第一锚为空时不会静默换锚。

Chain 示例显式选择源块 1，即有来源的 NR_046654.1/chr3 关系；上游 `_modified` 格式对照仍在未改动输入中，不冒充新生物学观测。三个 exon 对齐块保留中间缺口。Chain 不提供错配或字符单元，chain 分数不是 Identity。Query 映射 qaccver，target 映射 saccver，原 ID 独立可追溯。默认 bp 不适用时应显式声明 `sequence_unit`。

两套组装都叫 `chr1` 时须声明身份，例如 `query_ids = c(chr1 = "assemblyB.chr1")`、`subject_ids = c(chr1 = "assemblyA.chr1")`。默认保留 ID，并拒绝跨角色歧义，不自动加前缀或合并组装。完整长度冲突明确失败；最终共享 ID 要求两侧显式且一致地声明。这是源身份，与显示实例分开。两种格式都不会仅由端点推断 pident、支持计数或完整源字符串。

## 已有 delta 比对 {#delta}

`read_delta()` 需要已有 **show-coords** 程序；本章节可通过 `DELTA_DECODER` 选择其可执行文件。读取的是已有 NUCMER/PROMER 结果，不调用 aligner，也不打开序言中的历史 FASTA 路径。未改动的上游 Drosophila 示例、固定版本、原始许可及哈希保留在 `examples/genome-alignments/`。程序缺失时跳过可选示例，不代表成功解码。

{{code:courseblockbc}}

这里显式选取的 query 仍是字面字符串 `"3210101"`。Query/reference 映射 qaccver/saccver，原一基有方向端点不再次平移。每条原生输出先与实际源表头和端点核对，再关联错误计数、有符号 delta 脚本和物理行号。完整长度来自源声明，不取最后比对端点。跨角色同名源遵循 `read_chain()` 的显式 ID 映射，不猜组装身份。

`pident` 为工具报告的百分比及其打印精度，不由两个端点重建。PROMER 中它是 **AA 比对统计**，而区间坐标和完整源长度仍是 **bp**。原生 similarity/stop 百分比及 frame 独立保留。qspan/sspan 是基因组跨度，不是对齐列数；不制造 BLAST `length`、read support 或字符。原文、原生输出及消息保留在 `report`。若程序没有可靠版本选项，读取报告中的版本保持未知，应另行记录准备环境的来源。定向检查记录了实际固定工具源码 tag 和可执行文件校验和。ggchord 不重新分发第三方程序或实现。

来源输入包还包含 `arabidopsis-first-pair.delta`：保留原序言及首个完整拟南芥序列对块的全部 39 条记录，字节为原文件的未改动前缀。`arabidopsis-source.json` 记录固定提交的原始压缩文件、解压哈希、完整块规则及产物哈希，并保留作者的 MIT 许可原文。它提供独立于翻译比对示例的真实核酸比对验证，不重新运行比对，也不认证剩余源块。复用同一个已有解码器：

{{code:courseblockbd}}

## 已有蛋白 signature {#protein-signatures}

安装包内的输入为从 InterPro API 缓存的两个普通人类蛋白的已有 Pfam 注释。原 JSON 字节、获取日期、URL、CC0 归属及哈希在 `inst/extdata/signature-sources.json` 中。`data-raw/prepare-signature-inputs.py` 将相同已观测片段显式重编码为 TSV，供坐标往返核对；不运行 InterProScan，缺失的 MD5、状态、运行日期及整合描述保持 `-`。TSV 是派生格式示例，不是新分析结果。

{{code:courseblockbe}}

完整长度为声明的 AA 长度，不取 feature 最大端点。API 返回的小写 ID（`p00517`、`p51587`）、accession version 及重复模型命中均不改动。模型/signature ID 不是唯一 feature 实例 ID，因此 `feature_id` 保持缺失；与其他文件的源身份协调须显式完成。

每个缓存 JSON location 为一个 feature，实际片段存入 `segments`，start/end 仅包围其 footprint。现有 Feature 几何保留坐标空隙，并维持一个源命中。代表性标志是元数据，不自动筛选。原 entry/protein/location 索引及文件内 `.source_match` 可追溯，多文件须结合 `.source_input`。TSV 没有 JSON 分组证据，不按模型 ID 自动合并。分数保留原生列表列：JSON 数字、TSV 字符串和缺失值保留各自来源类型。不同库的 score/e-value 不能当作一种通用置信度比较，也不制造 identity/read support。若要使用数值标度，应先选择所需成员库并明确其分数含义。

JSON 需要可选 **jsonlite**，只读提供的本地缓存。后续页 URL 和声明计数留在报告中，绘图时不获取数据。质量/上下文元数据、原始字段及完整源对象可供核对。这些位置注释不提供完整蛋白字符，也不构成独立字符核验。[官方格式说明](https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html)及[分数解释](https://interpro-documentation.readthedocs.io/en/latest/faq.html#why-are-there-no-e-values-associated-with-interpro-entries)说明其底层约定。

## 已有 feature 关系 {#feature-relations}

`read_feature_relations()` 离线读取已有成对关系或组成员。安装包中的植物子集保留固定源提交中 MCScanX 数据的完整原始行。许可、原 URL、哈希及抽样物理行号见 `inst/extdata/feature-relation-sources.json`；`data-raw/prepare-feature-relations.py` 可重现输入，不运行分析。

{{code:courseblockbf}}

调用者声明升序一基闭合连续锚点，文件后缀不能确立坐标约定。配套表没有完整染色体长度，因此本例检查关系，不伪造 `seq_data` 绘制整条染色体。绘制 `ribbon_data` 区间关系前，应独立取得完整长度；这些关系不是逐位点 HSP。只有明确分数含义后才映射分数；原第三列保持不透明字符串，不当作 Identity 或 e-value。

`pairs` 接受两个 ID 及可选分数；`named_groups` 将第一列保留为原生字符组 ID，包括数字状 ID；`groups` 将全部字段作为成员，仅记录文件内组索引。成员关系不会展开为全体两两连接。重复、单成员组、成员顺序、原字段及物理行均保留。

可选 `feature_data` 按唯一 `feature_id` 精确匹配，不修改大小写、版本或前缀。成对示例与组示例有不同原始命名空间及配套位置表。缺失 ID 默认报错，`unmatched = "report"` 保留原行并标记不完整映射；重复 feature ID 的歧义始终报错。没有位置时仍返回全部关系，不生成映射连接。注释链向单独保留，不据此赋予比对方向；分段注释需要显式连续锚点表。不推断完整长度、Identity、支持度或源字符。

## 图路径：探索性子集 {#graph-paths}

GFA 图边表达邻接/重叠，不自动提供整基因组间 HSP。已有显式路径步骤在坐标构建有据时，可用普通 sequence/feature 表展示。节点 ID 可在同一路径重复，须保留路径、源步骤索引、出现次序及遍历方向，不能连接所有出现位置，也不能把节点/边数当作直系同源、Identity 或 read support。

真实探索例使用固定提交中未改动的 ODGI 人类 LPA 图，由可选 **gfapy** 解码。5,195 条原边都是零重叠；原型核对每个原始路径转移和节点长度，投影 13 条给定路径的 202,806 个步骤。展示子集显式保留两条命名路径各前 100 个步骤，原遍历方向及 ID 不改；给定完整路径长度与染色体长度严格区分，不生成 pair 连接。代码、源 URL/哈希、原 MIT 许可、输出表及实际图见源码 checkout 的 `tools/long-term-audit/OPTIONAL-GRAPH-BRIDGE-v0.27.0.md` 及其原型资产；`data-raw/prepare-graph-prototype.py` 仅在显式调用时取得两个固定原文件。

这不是通用 `read_gfa()` 或 graph Coord。非零/未知重叠、源字符缺失、GFA W/J/GFA2 及 ODGI 原生二进制不在此原型范围。用成熟上游工具导出明确路径子集并保留日志/来源映射，不推断缺失坐标或成对证据；文档中的上游导出路线也不等同于此处已验证二进制解码。
