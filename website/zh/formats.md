<div class="task-intro"><p>先根据手中的数据选择入口。读入文件后，核对返回的 ID、长度和端点，再交给绘图层。</p></div>

## 手中的文件从哪里开始 {#choose-format}

| 手中的文件 | 下一步 |
| --- | --- |
| FASTA 序列 + GFF3 注释 + BLAST 局部比对 | [从源文件画比较图](real-data.html#read) |
| BED 区间或 PAF（Pairwise mApping Format）比对 | [核对零基坐标](#zero-based) |
| GenBank 单序列及注释 | [保留原始记录与分段](#genbank) |
| bedGraph 或 BigWig 连续信号 | [显式保留值和缺失区间](#signals) |
| VCF/BCF 变异和结构断点 | [保留原始变异证据](#variants) |
| 已有蛋白 signature 注释 | [保留 AA 位置及源分数](#protein-signatures) |
| 已有 feature 对及组 | [保留原关系与成员](#feature-relations) |
| BAM/CRAM 已比对 reads | [使用成熟覆盖度桥接](#aligned-reads) |
| 已有 delta 结果 | [区分原生统计和坐标单位](#delta) |
| 已有 MAF 或 chain 基因组比对 | [保留块及源坐标](#genome-alignments) |
| 已有 R 表格 | [映射必需字段](#adapters) |
| GRanges 或 GInteractions | [直接转换区间对象](#ranges) |

## 本页准备 {#prepare}
文件示例使用仓库中的小型真实记录。先解压[教学输入](../assets/tutorial-inputs.zip)，使工作目录包含 `examples/`。

{{prepare:setup}}

## 把已有表格转换成需要的列 {#adapters}
这里使用随包的实际序列、连接和基因，
只重命名列来演示角色映射；没有额外生成示意序列或比对。

{{code:adapters}}

{{console:adapters}}

转换返回普通数据框，并记录 `ggchord_conversion` 信息。标识原样保留；缺少完整
序列长度时应补充可靠来源，不能把最大注释端点当作全长。核酸与蛋白比对的最小列
和单位见[导入比对结果](alignments.html)。

## BED 和 PAF：坐标只转换一次 {#zero-based}
BED 的起点从 0 开始，终点不包含在区间内。下例读入一条 pUC19 复制起点记录，比较读入前后的坐标。

{{code:formats}}
{{output:formats}}

PAF 提供两端全长和区间。先选出 mapq ≥ 30 的比对，再由两端长度生成骨架表。

{{code:formatpaf}}
{{output:formatpaf}}

{{code:drawpaf}}
{{figure:paf}}

反向比对用递减的 subject 端点表达。绘图使用读取器转换后的坐标。

## GenBank：保留记录身份 {#genbank}
{{code:genbank}}
{{console:genbank}}

这个文件是 NCBI 的原始 L09137.2 记录，名称是 pUC19c；不要把它强行改成随包 pUC19。
返回值包括 `seq_data`、`features$annotation` 与 `report`，分段位置及记录/注释来源
保留。这里需要所选 Python 已安装 Biopython；函数不自动安装依赖。
`unsupported = "report"` 将不支持的位置记入报告，不把远端或复杂位置猜成区间。

本例已在启用 Python Biopython 1.85 的环境重跑，实际返回 L09137.2、2686 bp 及一行注释。可在独立文档环境使用对应 Python 的 `python3 -m pip install biopython` 安装，并在启动 R 前将其可执行文件加入 PATH。可选依赖导致跳过或返回 NULL 时，不应记作导入成功。

显式 `python` 路径保留选定虚拟环境，包括解释器符号链接；不安装或改动环境。

## 可选的区间对象 {#ranges}

第一段把序列 a 上的 10–30 bp 转成特征表。第二段把 a 上的 10–30 bp 与 b 上的 40–60 bp 连成一条记录。GRanges 表示一个区间，GInteractions 表示一对区间。先安装 GenomicRanges 和 InteractionSet，完整序列长度由 seqlengths 提供。

{{code:frontimportgr}}

{{output:frontimportgr}}

{{code:frontimportgi}}

{{output:frontimportgi}}

## BED双端区间 {#bedpe}

`read_bedpe()`读取六个必填的制表符分隔字段，以及可选name、score、两条链和额外元数据。第一端保持query、第二端保持subject；score保留为字符串，不推断Identity或测序支持度。格式依据[BEDPE说明](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)。

下面将随包已观测比对的区间转写为BEDPE，用于坐标往返验证，不将它们当作配对测序证据。六列输入没有链证据，因此导入后的锚点链和相对方向保持未知。有真实链信息时才提供对应列，不由坐标升序猜测链。

```r
data(ribbon_data_example)
paired_rows <- ribbon_data_example |>
  transmute(
    chrom1 = qaccver, start1 = pmin(qstart, qend) - 1, end1 = pmax(qstart, qend),
    chrom2 = saccver, start2 = pmin(sstart, send) - 1, end2 = pmax(sstart, send)
  )
paired_file <- tempfile(fileext = ".bedpe")
write.table(paired_rows, paired_file,
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
paired_links <- read_bedpe(paired_file)
unlink(paired_file)
paired_links |>
  select(qaccver, saccver, qstart, qend, sstart, send, q_mapped, s_mapped)
```

负链锚点采用降序canonical端点，与配对区间适配器一致；同时保留原始start/end和链字段。name不视为唯一ID。额外列在`bedpe_extra`中保留，按`column_11`起命名。未知染色体（`.`）或边界（`-1`）默认报错；显式`unmapped="keep"`为审计保留不完整行，适配或绘图前须筛选`q_mapped & s_mapped`。非法已知坐标和零宽区间始终报错。区间Ribbon没有Identity时，显式设置固定填充，不制造pident。

## 完成后检查 {#check-result}
先查看返回列与转换报告，再选绘图列。PAF 反向命中保留端点方向；GenBank 记录身份不能因为展示需要而改成另一条随包质粒。

<div class="next-step">继续：<a href="alignments.html">导入核酸和蛋白比对 →</a></div>


## 保留真实 FASTA 字符 {#characters}

通过 `keep_sequence = TRUE` 保留真实字符、头部物理行号和可选文件来源。仅读取时去除空白，保留大小写；显式声明 bp 或 aa，不猜字母表。聚焦同时裁剪字符串与区间，下例局部位置 1 对应原始位置 1001。降序区域仍选择同一区间，不执行反向互补。只有长度不能恢复字符，只有 HSP 端点不能恢复 mismatch 或 gap 列。有 gap 的 FASTA 会将 gap 计入长度，基因组用途需要单独的比对列到源位置映射。

```r
character_source <- read_fasta_lengths(
  "examples/plasmid/pUC19.fna",
  sequence_unit = "bp", keep_sequence = TRUE, source_file = TRUE
)
character_locus <- focus_ggchord_data(
  character_source,
  loci = data.frame(accver = character_source$accver, start = 1001, end = 1060)
)
character_locus$seq_data |>
  select(accver, length, sequence, .source_origin)
```


## 已有字符对齐：区分对齐列与源位置 {#aligned-characters}

以下示例需要源码 checkout 中的 `examples/alignment/`。该目录保留 Biopython 原始 DNA/蛋白比对输入、上游许可、哈希及可重复的 FASTA 格式转写规则。不运行比对工具，不改动完整记录 ID 或原始字符；显式选择行与列。

`read_aligned_fasta()` 读取一个已有对齐，须声明 DNA、RNA 或 AA；等长本身不能证明同源性。`alignment_column` 包括 `-` 和 `.` gap，这些列的 `source_position` 为 NA。默认源位置在去 gap 输入片段中从 1 计数，不是完整基因组坐标。小写字符也消耗位置，因此不能把此读取器当作 A3M 插入语义解析器。

```r
aligned_dna <- read_aligned_fasta(
  "examples/alignment/opuntia-aligned.fasta", "DNA", source_file = TRUE
)
selected_rows <- unique(aligned_dna$accver)[c(1, 3, 7)]
local_dna <- aligned_dna |>
  filter(accver %in% selected_rows, between(alignment_column, 321, 380))

map_alignment_position(local_dna, selected_rows[1], c(342, 343, 344))
map_alignment_position(local_dna, selected_rows[1], 343, inverse = TRUE)

p_characters <- ggplot2::ggplot(
  local_dna, aes(alignment_column, accver, label = symbol)
) +
  geom_residue(aes(fill = symbol), size = 3, overflow = "error") +
  ggplot2::scale_fill_manual(
    values = c(
      A = "#B8D8BA", C = "#A6CBE3", G = "#F2D28C", T = "#EAA8A8",
      N = "#D8D8D8", "-" = "#F5F5F5", "." = "#F5F5F5"
    ),
    guide = "none"
  ) +
  ggplot2::labs(x = "Alignment column", y = NULL) +
  ggplot2::theme_minimal()
p_characters
```

首条选中记录在对齐列 343 有 gap；按列查找返回 `matched = TRUE`、`is_gap = TRUE` 和缺失源位置。反向查找源位置 343 则定位到其实际对齐列，不做插值。被筛选删除的位置返回 `matched = FALSE`；筛选不重编号保留列。即使另行声明了不同源 ID，`map_alignment_position()` 仍按对齐行 ID 选择记录。

`geom_residue()` 是普通 Cartesian ggplot 图层，样式和用户标度遵循 ggplot2。绘制时按实际单元宽高测量字形：`overflow = "error"` 要求文字可读，`"tile"` 对放不下的字母保留色格，不缩小字号。显示字母时可减少列数或增大物理宽度，例如用普通 `ggsave(width = 12, height = 3.4, units = "in")` 保存 `p_characters`；仅增加 DPI 不会增加物理空间。支持离散行、identity/reverse 位置标度与 Cartesian 坐标；该层不直接投射到骨架 Coord。具有明确来源映射的 Cartesian 字符图可经 `register_ggchord_alignment()` 注册，用于 [Fusion 字符细节](fusion.html#character-detail)。

### 蛋白局部片段 {#aligned-protein}

同一数据契约可接收真实蛋白字符。该窗口保留上游 gap 与残基差异；fill 区分原始字符，不代表推断的功能域。

```r
aligned_protein <- read_aligned_fasta(
  "examples/alignment/protein-aligned.fasta", "AA"
)
protein_rows <- unique(aligned_protein$accver)[2:4]
local_protein <- aligned_protein |>
  filter(accver %in% protein_rows, between(alignment_column, 221, 260))

ggplot2::ggplot(local_protein, aes(alignment_column, accver, label = symbol)) +
  geom_residue(aes(fill = symbol), size = 3, overflow = "tile") +
  ggplot2::scale_fill_hue(l = 85, c = 35, guide = "none") +
  ggplot2::labs(x = "Alignment column", y = NULL) +
  ggplot2::theme_minimal()
```

### 显式源坐标 {#aligned-source}

`source_data` 的每个对齐记录对应一行，包含 `accver`、`start` 和 `strand`。`start` 是首个非 gap 源位置；负链位置递减，字符串必须已处于比对方向。可选 `source_accver` 标识原分子。若提供完整正向 `sequence`，逐一核对非 gap 字符：负链核酸使用互补，AA 仅反序。未提供完整源字符串时标为未核对，不猜起点或链。

下例只检查同一对齐去 gap 后的片段一致性，**不构成独立基因组验证**。实际基因组投影应提供有独立来源的序列和区间元数据。

```r
original_strings <- read_fasta_lengths(
  "examples/alignment/opuntia-aligned.fasta", keep_sequence = TRUE
)
fragment_sources <- original_strings |>
  transmute(
    accver, source_accver = accver, start = 1, strand = "+",
    sequence = gsub("[-.]", "", sequence)
  )
checked_cells <- as_ggchord_alignment(
  original_strings, "DNA", source_data = fragment_sources
)
```

仅凭六个 HSP 端点不能恢复内部 mismatch/gap。此工作流接收已有对齐列，不解析 CIGAR/A3M/Clustal/Stockholm，不推断相位，不计算共识、Logo 或新比对。显式CDS翻译见下节；精确宏观/细节注册见Fusion课程。

## 有来源的 CDS 翻译 {#cds-translation}

`build_ggchord_cds()` 定义已观测编码坐标；`translate_ggchord_cds()` 还需要真实的完整正向 DNA/RNA 字符串及可选依赖 **Biostrings**。本例真实 EMALE 负链 CDS 的源 phase 为0、起始密码子为ATG。表11为显式选择，并与保留的上游蛋白 FASTA 核对；输入来源见 `examples/format-inputs/README.md`。

```r
library(dplyr)
# Requires the source checkout, dplyr and optional Biostrings.
genome <- read_fasta_lengths(
  "examples/format-inputs/emales-BVI_023A.fna", keep_sequence = TRUE
)
coding <- read_gff3("examples/format-inputs/emales-BVI_023A.gff3") |>
  filter(feature_id == "BVI_023A_0009")
cds <- build_ggchord_cds(coding, genetic_code = 11)
observed <- translate_ggchord_cds(cds, genome, initiation = "start")
expected <- read_fasta_lengths(
  "examples/format-inputs/protein-pair.faa", keep_sequence = TRUE
) |>
  filter(accver == "BVI_023A_0009")
stopifnot(sub("[*]$", "", observed$models$protein_sequence) == expected$sequence)
observed$residues |>
  select(alignment_column, symbol, codon, source_positions) |>
  head()
```

`models` 报告编码/蛋白字符串、翻译表、起始及模糊碱基策略、不完整末端和末端余数。`nucleotides` 保留每个基因组碱基；`residues$source_positions` 按生物学顺序列出每个 codon 的三个基因组位置，也保留跨 exon 的 codon。残基位置使用 CDS feature ID 的已观测蛋白坐标域，不推断外部蛋白 accession。内部及末端终止符均保留为 `*`；本例仅在对照时删除一个末端终止符，以符合上游 FASTA 约定。

默认 `initiation = "literal"` 不改写替代起始密码子；`"start"` 要求完整5′端。`ambiguous = "error"` 拒绝模糊 codon，显式 `"solve"` 交由 Biostrings 解析，不能唯一确定时输出X。未知翻译表、缺失/过短源序列及已声明的翻译例外均报错。`partial = "trim"` 在构建 CDS map 时声明：只排除首段 phase 与末端不完整 codon，后续 phase 用于校验连续性，不删除跨段 codon 的碱基。缺失序列不会被补造。RNA 必须显式声明；有序 `segments` 沿用既有分段/跨原点契约。

完整的碱基、codon 和宏观关联可运行示例在 `examples/m2-detail.R`，负链氨基酸标在其三个编码碱基中心。真实正/负链分段对照、源hash和原始预期翻译保存在 `examples/cds/`。[注册字符图](fusion.html#character-detail)后可连接精确源边界。

## 连续信号：bedGraph 与 BigWig {#signals}

这个安装包内的示例由未改动的随包 pUC19 序列按 100 bp 窗口计算真实 GC 百分比，末窗为 86 bp；它是序列组成，不是实验测序覆盖度。`data-raw/prepare-signal-inputs.R` 保留生成规则，`inst/extdata/signal-sources.json` 保留哈希及显式的原始表头映射。读取器保留标识。

```r
signal_file <- system.file("extdata", "pUC19-gc-100.bedGraph", package = "ggchord")
signal <- read_bedgraph(signal_file, source_file = TRUE)
data(plasmid_example_pUC19)
gc_band <- position_track("GC", side = "inside", width = unit(7, "mm"))
p_signal <- ggchord(plasmid_example_pUC19) +
  geom_seq(show.legend = FALSE) +
  geom_track(
    aes(xmin = start - 1, xmax = end, ymin = 0, ymax = value, fill = value),
    data = signal, geom = "rect", track = gc_band
  ) +
  scale_track_value_continuous("GC", limits = c(0, 100), unit = "GC (%)") +
  scale_track_fill_gradientn(
    "GC", colours = c("#F1E51D", "#31849B"), limits = c(0, 100), name = "GC (%)"
  ) +
  coord_circular()
p_signal

if (.Platform$OS.type != "windows" && requireNamespace("rtracklayer", quietly = TRUE)) {
  indexed_file <- system.file("extdata", "pUC19-gc-100.bw", package = "ggchord")
  selected_signal <- read_bigwig(
    indexed_file,
    which = data.frame(accver = "pUC19", start = c(50, 100), end = c(150, 250)),
    source_file = TRUE
  )
  selected_signal |> select(accver, start, end, value)
}
```

bedGraph 为四列制表符分隔数据，使用零基半开坐标及有限的有符号数值。导入仅转换一次坐标；原始端点、值、物理行号和不解释样式的显示表头均可追溯。每文件一条信号，多个 track 定义应拆分。字面序列 ID `track`、`browser` 合法，表头不设置图形样式。重叠行保留，不隐式采用末行覆盖规则。

矩形以 `start - 1` 和 `end` 为单元边界，不在缺失区间插值连线。BigWig 需要可选 **rtracklayer**，其解码器不支持 Windows；该平台使用已有 bedGraph 文件。二进制导入的 GRanges 已是一基闭区间。区段值保留二进制浮点精度，区段数量不代表原始测量次数。完整长度来自文件表头。查询采用无链区间并集，结果裁切至并集，空选区返回零行。缺失区间仍缺失，与存储的零不同；不制造 coverage、identity、support、单位或裁切前观测 ID。各 Track 的映射、数值域、单位和 Guide 仍显式配置。

## 变异记录与结构断点 {#variants}

VCF 需要可选 **VariantAnnotation**，显式结构变异解码另外需要 **StructuralVariantAnnotation**。示例读取依赖包未改动的有来源输入：chr22 的 1000 Genomes 子集，以及 StructuralVariantAnnotation 分发的真实 caller VCF。数据保留在这些包中，ggchord 不复制第三方解析器实现或输入。`check-variant-import.R` 记录确切依赖版本和源哈希。

```r
if (requireNamespace("VariantAnnotation", quietly = TRUE)) {
  variant_file <- system.file("extdata", "chr22.vcf.gz", package = "VariantAnnotation")
  imported_variants <- read_vcf(variant_file, source_file = TRUE)
  imported_variants$features$variant |>
    select(accver, start, end, feature_id, REF, ALT) |>
    head()
}

if (requireNamespace("StructuralVariantAnnotation", quietly = TRUE)) {
  structural_file <- system.file("extdata", "gridss.vcf", package = "StructuralVariantAnnotation")
  structural <- read_vcf(structural_file, structural = TRUE)
  structural$links |>
    select(qaccver, saccver, qstart, qend, sstart, send,
      q_breakend_orientation, s_breakend_orientation)
  structural$report$records[[1]]$warnings
}
```

每条可放置记录保留一个 feature 行及 ALT 列表，多等位调用也不拆分。Feature 表示包括 padding 碱基的参考等位基因跨度，不是推断的符号 SV footprint。完整 INFO、样本/基因型数组及表头在 `report$records[[1]]$parsed` 中，应在那里检查实际 caller 字段；缺失支持度或 Identity 不变成零，也不制造计数。未知完整长度保持 NA，绘图前应提供可靠参考长度。端粒和越界跨度保留在 `unplaced_features` 及源对象中，不裁切。

结构解码需显式启用，保留置信区间，关闭缺失 mate 推断；每个互惠对生成一条连接，包含两个原始记录索引。警告标明解码器未配对的记录，报告保留全部记录。`.plottable = FALSE` 标记超出已知参考边界的断点，相关配对保留在 `unplaced_links`。精确单碱基断点可显式选取位置后用于 `geom_link_line()`；不确定区间适合区间几何。`q_breakend_orientation`、`s_breakend_orientation` 描述衍生序列，不是相对同源链。不能制造比对 Identity 为这些连接配色。

现代 BCF（`.bcf` 或 `.bcf.gz`）通过 `python` 指定已安装可选 **pysam** 的现有解释器，例如 `read_vcf("calls.bcf", python = "/path/to/environment/bin/python")`。仅 BCF 会调用它，经 HTSlib 转临时 VCF，并记录解码器、版本和消息。虚拟环境解释器路径保留其环境。原记录顺序、数值及基因型保留；转码文本不是原始 BCF 物理行。不会安装任何软件。本读取器将选定本地文件加载到内存，大数据应在上游选取或子集化。

## 基因组对齐块：MAF 与 chain {#genome-alignments}

使用已有 Python 环境，其中须安装支持 MAF/chain 的可选 **Biopython Bio.Align**，已直接验证 1.85。课程可通过 `ALIGNMENT_PYTHON` 选择该解释器。源码 checkout/输入包须包含 `examples/genome-alignments/`；两个未改动的固定上游版本文件有哈希、许可归属及验证脚本，不运行 aligner。可选解析器不可用意味着跳过示例，不是成功导入。

```r
alignment_python <- Sys.getenv("ALIGNMENT_PYTHON", "python3")
alignment_available <- nzchar(Sys.which(alignment_python)) &&
  identical(suppressWarnings(system2(
    alignment_python,
    c("-c", shQuote("from Bio import Align; assert hasattr(Align, 'parse')")),
    stdout = FALSE, stderr = FALSE
  )), 0L)
if (alignment_available) {
  maf <- read_maf(
    "examples/genome-alignments/ucsc_mm9_chr10.maf",
    alphabet = "DNA", python = alignment_python, source_file = TRUE
  )
  local_maf <- maf$alignment |>
    filter(.source_input == 1, source_block == 1, between(alignment_column, 61, 100))
  ggplot2::ggplot(local_maf, aes(alignment_column, accver, label = symbol)) +
    geom_residue(aes(fill = symbol), size = 3, overflow = "tile") +
    ggplot2::scale_fill_manual(
      values = c(
        A = "#B8D8BA", C = "#A6CBE3", G = "#F2D28C", T = "#EAA8A8",
        N = "#D8D8D8", "-" = "#F5F5F5"
      ), guide = "none"
    ) +
    ggplot2::labs(x = "MAF alignment column", y = NULL) +
    ggplot2::theme_minimal()

  chains <- read_chain(
    "examples/genome-alignments/dna_rna.chain", python = alignment_python
  )
  observed_chain <- chains$ribbon_data |>
    filter(source_block == 1)
  observed_chain |>
    select(qaccver, saccver, qstart, qend, sstart, send)
}
```

MAF 保留实际带 gap 字符和声明的完整源长度，负链分量映射至递减的**正向源坐标**。第一分量锚定关系、源块/文件及原始坐标矩阵均显式保留，不自动扩展成所有可能配对。Gap 单元的源位置缺失，完整长度不能独立验证字符。对齐列在每块重新开始，因此 `map_alignment_position()` 或字符面板前必须筛选一个文件/块。原质量、上下文及空分量证据保留在源与报告中，第一锚为空时不会静默换锚。

Chain 示例显式选择源块 1，即有来源的 NR_046654.1/chr3 关系；上游 `_modified` 格式对照仍在未改动输入中，不冒充新生物学观测。三个 exon 对齐块保留中间缺口。Chain 不提供错配或字符单元，chain 分数不是 Identity。Query 映射 qaccver，target 映射 saccver，原 ID 独立可追溯。默认 bp 不适用时应显式声明 `sequence_unit`。

两套组装都叫 `chr1` 时须声明身份，例如 `query_ids = c(chr1 = "assemblyB.chr1")`、`subject_ids = c(chr1 = "assemblyA.chr1")`。默认保留 ID，并拒绝跨角色歧义，不自动加前缀或合并组装。完整长度冲突明确失败；最终共享 ID 要求两侧显式且一致地声明。这是源身份，与显示实例分开。两种格式都不会仅由端点推断 pident、支持计数或完整源字符串。

## 已有 delta 比对 {#delta}

`read_delta()` 需要已有 **show-coords** 程序；本课程可通过 `DELTA_DECODER` 选择其可执行文件。读取的是已有 NUCMER/PROMER 结果，不调用 aligner，也不打开序言中的历史 FASTA 路径。未改动的上游 Drosophila 示例、固定版本、原始许可及哈希保留在 `examples/genome-alignments/`。程序缺失时跳过可选示例，不代表成功解码。

```r
delta_decoder <- Sys.getenv("DELTA_DECODER", "show-coords")
if (nzchar(Sys.which(delta_decoder))) {
  delta <- read_delta(
    "examples/genome-alignments/promer.delta",
    executable = delta_decoder, source_file = TRUE
  )
  observed_intervals <- delta$ribbon_data |>
    filter(qaccver == "3210101")
  displayed_ids <- unique(c(observed_intervals$qaccver, observed_intervals$saccver))
  displayed_sources <- delta$seq_data |>
    filter(accver %in% displayed_ids)

  ggchord(displayed_sources) +
    geom_seq(show.legend = FALSE) +
    geom_seq_label() +
    geom_link_ribbon(data = observed_intervals) +
    coord_collinear() +
    labs(
      title = "Existing translated alignment intervals",
      subtitle = "Drosophila source coordinates in bp; native identity in AA"
    )
}
```

这里显式选取的 query 仍是字面字符串 `"3210101"`。Query/reference 映射 qaccver/saccver，原一基有方向端点不再次平移。每条原生输出先与实际源表头和端点核对，再关联错误计数、有符号 delta 脚本和物理行号。完整长度来自源声明，不取最后比对端点。跨角色同名源遵循 `read_chain()` 的显式 ID 映射，不猜组装身份。

`pident` 为工具报告的百分比及其打印精度，不由两个端点重建。PROMER 中它是 **AA 比对统计**，而区间坐标和完整源长度仍是 **bp**。原生 similarity/stop 百分比及 frame 独立保留。qspan/sspan 是基因组跨度，不是对齐列数；不制造 BLAST `length`、read support 或字符。原文、原生输出及消息保留在 `report`。若程序没有可靠版本选项，读取报告中的版本保持未知，应另行记录准备环境的来源。定向检查记录了实际固定工具源码 tag 和可执行文件校验和。ggchord 不重新分发第三方程序或实现。

来源输入包还包含 `arabidopsis-first-pair.delta`：保留原序言及首个完整拟南芥序列对块的全部 39 条记录，字节为原文件的未改动前缀。`arabidopsis-source.json` 记录固定提交的原始压缩文件、解压哈希、完整块规则及产物哈希，并保留作者的 MIT 许可原文。它提供独立于翻译比对示例的真实核酸比对验证，不重新运行比对，也不认证剩余源块。复用同一个已有解码器：

```r
if (nzchar(Sys.which(delta_decoder))) {
  plant_delta <- read_delta(
    "examples/genome-alignments/arabidopsis-first-pair.delta",
    executable = delta_decoder, source_file = TRUE
  )
  plant_delta$ribbon_data |>
    select(qaccver, qstart, qend, saccver, sstart, send, pident) |>
    head()
}
```

## 已有蛋白 signature {#protein-signatures}

安装包内的输入为从 InterPro API 缓存的两个普通人类蛋白的已有 Pfam 注释。原 JSON 字节、获取日期、URL、CC0 归属及哈希在 `inst/extdata/signature-sources.json` 中。`data-raw/prepare-signature-inputs.py` 将相同已观测片段显式重编码为 TSV，供坐标往返核对；不运行 InterProScan，缺失的 MD5、状态、运行日期及整合描述保持 `-`。TSV 是派生格式示例，不是新分析结果。

```r
signature_file <- system.file("extdata", "P00517-signatures.tsv", package = "ggchord")
signatures <- read_signature_matches(signature_file, source_file = TRUE)
signatures$features$signature |>
  select(accver, start, end, signature_id, analysis, source_score)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  cached_file <- system.file("extdata", "P51587-pfam.json", package = "ggchord")
  cached <- read_signature_matches(cached_file, format = "entry_json")
  ggchord(cached$seq_data) +
    geom_seq(show.legend = FALSE) +
    geom_seq_label() +
    geom_feature(
      aes(feature_fill = signature_id), data = cached$features$signature,
      shape = "block", position = position_feature_stack()
    ) +
    coord_collinear()
}
```

完整长度为声明的 AA 长度，不取 feature 最大端点。API 返回的小写 ID（`p00517`、`p51587`）、accession version 及重复模型命中均不改动。模型/signature ID 不是唯一 feature 实例 ID，因此 `feature_id` 保持缺失；与其他文件的源身份协调须显式完成。

每个缓存 JSON location 为一个 feature，实际片段存入 `segments`，start/end 仅包围其 footprint。现有 Feature 几何保留坐标空隙，并维持一个源命中。代表性标志是元数据，不自动筛选。原 entry/protein/location 索引及文件内 `.source_match` 可追溯，多文件须结合 `.source_input`。TSV 没有 JSON 分组证据，不按模型 ID 自动合并。分数保留原生列表列：JSON 数字、TSV 字符串和缺失值保留各自来源类型。不同库的 score/e-value 不能当作一种通用置信度比较，也不制造 identity/read support。若要使用数值标度，应先选择所需成员库并明确其分数含义。

JSON 需要可选 **jsonlite**，只读提供的本地缓存。后续页 URL 和声明计数留在报告中，绘图时不获取数据。质量/上下文元数据、原始字段及完整源对象可供核对。这些位置注释不提供完整蛋白字符，也不构成独立字符核验。[官方格式说明](https://interproscan-docs.readthedocs.io/en/v5/OutputFormats.html)及[分数解释](https://interpro-documentation.readthedocs.io/en/latest/faq.html#why-are-there-no-e-values-associated-with-interpro-entries)说明其底层约定。

## 已有 feature 关系 {#feature-relations}

`read_feature_relations()` 离线读取已有成对关系或组成员。安装包中的植物子集保留固定源提交中 MCScanX 数据的完整原始行。许可、原 URL、哈希及抽样物理行号见 `inst/extdata/feature-relation-sources.json`；`data-raw/prepare-feature-relations.py` 可重现输入，不运行分析。

```r
pair_file <- system.file("extdata", "plant-feature-pairs.homology", package = "ggchord")
group_file <- system.file("extdata", "plant-feature-groups.txt", package = "ggchord")
position_file <- system.file("extdata", "plant-pair-features.gff", package = "ggchord")

# 此工具的四列位置表不是 GFF3；显式声明源锚点。
anchors <- read.table(
  position_file, header = FALSE, na.strings = character(),
  col.names = c("accver", "feature_id", "start", "end"),
  colClasses = c("character", "character", "numeric", "numeric")
)
relations <- read_feature_relations(
  pair_file, feature_data = anchors, source_file = TRUE
)
relations$ribbon_data |>
  select(qfeature_id, sfeature_id, qaccver, qstart, qend, source_score)

groups <- read_feature_relations(group_file, format = "named_groups")
groups$membership |>
  select(group_id, feature_id, .source_group, .source_member)
```

调用者声明升序一基闭合连续锚点，文件后缀不能确立坐标约定。配套表没有完整染色体长度，因此本例检查关系，不伪造 `seq_data` 绘制整条染色体。绘制 `ribbon_data` 区间关系前，应独立取得完整长度；这些关系不是逐位点 HSP。只有明确分数含义后才映射分数；原第三列保持不透明字符串，不当作 Identity 或 e-value。

`pairs` 接受两个 ID 及可选分数；`named_groups` 将第一列保留为原生字符组 ID，包括数字状 ID；`groups` 将全部字段作为成员，仅记录文件内组索引。成员关系不会展开为全体两两连接。重复、单成员组、成员顺序、原字段及物理行均保留。

可选 `feature_data` 按唯一 `feature_id` 精确匹配，不修改大小写、版本或前缀。成对示例与组示例有不同原始命名空间及配套位置表。缺失 ID 默认报错，`unmatched = "report"` 保留原行并标记不完整映射；重复 feature ID 的歧义始终报错。没有位置时仍返回全部关系，不生成映射连接。注释链向单独保留，不据此赋予比对方向；分段注释需要显式连续锚点表。不推断完整长度、Identity、支持度或源字符。

## 已比对 reads：可选覆盖度桥接 {#aligned-reads}

使用可选 **GenomicAlignments**、**Rsamtools** 和 **GenomicRanges** 读取已有 BAM，显式选择过滤与计数口径，再由 ggchord 绘制信号。本例使用依赖包中未改动的果蝇 RNA-seq 子集，原 BAM header 保留 TopHat、参考及 SRR 来源；固定源 URL 和校验和见 `tools/long-term-audit/v0.27.0-read-coverage-results.json`。不运行 read 比对或变异检测。

```r
if (requireNamespace("GenomicAlignments", quietly = TRUE) &&
    requireNamespace("Rsamtools", quietly = TRUE)) {
  bam_file <- system.file("extdata", "sm_treated1.bam", package = "GenomicAlignments")
  read_filter <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(
      isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE,
      isNotPassingQualityControls = FALSE, isDuplicate = FALSE
    ), mapqFilter = 0L
  )
  aligned <- GenomicAlignments::readGAlignments(bam_file, param = read_filter)
  read_coverage <- GenomicAlignments::coverage(aligned, drop.D.ranges = TRUE)
  coverage_track <- bind_rows(lapply(names(read_coverage), function(id) {
    runs <- read_coverage[[id]]
    ends <- cumsum(S4Vectors::runLength(runs))
    data.frame(
      accver = id, start = c(1, head(ends, -1L) + 1), end = ends,
      value = S4Vectors::runValue(runs)
    )
  }))
  coverage_track |> head()
}
```

每条保留的**alignment 记录**在其已比对参考区块上贡献一次。跳过区（`N`）以及本例显式排除的删除区（`D`）贡献零；保留 supplementary alignment。不设碱基质量阈值，不去除重叠配对端或重复分子。这不是自动得到的 fragment、唯一分子、表达量、Identity 或变异支持度。零段是该明确计算中的已存储零，与实验信号缺失不同。大型输入可用依赖包的流式 `BamFile` coverage 方法，并沿用相同过滤声明。

将这些一基闭合区间交给 `geom_track(geom = "rect")`，设置 `xmin = start - 1`、`xmax = end`、`ymin = 0`、`ymax = value`，方法见[连续信号](#signals)。完整参考长度取 `scanBamHeader(bam_file)[[1]]$targets`，不取最后 read 端点。绘制小区域前显式选择源窗口。

CRAM 应由已有 HTSlib/pysam 或 samtools 环境解码，需要参考时提供**实际匹配的参考**，保留原 header、过滤和转换日志。不伪造参考、不自动下载未请求资源，也不宣称所有 CRAM 都可独立解码。定向原型仅用 `no_ref=1` 将相同真实 BAM 记录重编码，无外部参考访问，再核对解码后的全部记录字段/标签、参考字典及覆盖度。这是格式往返，不是原始 CRAM 观测，也不是所有参考型 CRAM 的认证。上游 `SO:sorted` header 和临时索引缺失产生的解码消息均保留；顺序读取成功，没有改写它们。

## 图路径：探索性子集 {#graph-paths}

GFA 图边表达邻接/重叠，不自动提供整基因组间 HSP。已有显式路径步骤在坐标构建有据时，可用普通 sequence/feature 表展示。节点 ID 可在同一路径重复，须保留路径、源步骤索引、出现次序及遍历方向，不能连接所有出现位置，也不能把节点/边数当作直系同源、Identity 或 read support。

真实探索例使用固定提交中未改动的 ODGI 人类 LPA 图，由可选 **gfapy** 解码。5,195 条原边都是零重叠；原型核对每个原始路径转移和节点长度，投影 13 条给定路径的 202,806 个步骤。展示子集显式保留两条命名路径各前 100 个步骤，原遍历方向及 ID 不改；给定完整路径长度与染色体长度严格区分，不生成 pair 连接。代码、源 URL/哈希、原 MIT 许可、输出表及实际图见源码 checkout 的 `tools/long-term-audit/OPTIONAL-GRAPH-BRIDGE-v0.27.0.md` 及其原型资产；`data-raw/prepare-graph-prototype.py` 仅在显式调用时取得两个固定原文件。

这不是通用 `read_gfa()` 或 graph Coord。非零/未知重叠、源字符缺失、GFA W/J/GFA2 及 ODGI 原生二进制不在此原型范围。用成熟上游工具导出明确路径子集并保留日志/来源映射，不推断缺失坐标或成对证据；文档中的上游导出路线也不等同于此处已验证二进制解码。
