{{code:courseblockaw}}

{{prepare:setup}}

先解压[真实输入包](../assets/tutorial-inputs.zip)，工作目录须含 examples/；每节单列可选依赖。

## 连续信号：bedGraph 与 BigWig {#signals}

这个安装包内的示例由未改动的随包 pUC19 序列按 100 bp 窗口计算真实 GC 百分比，末窗为 86 bp；它是序列组成，不是实验测序覆盖度。`data-raw/prepare-signal-inputs.R` 保留生成规则，`inst/extdata/signal-sources.json` 保留哈希及显式的原始表头映射。读取器保留标识。

{{code:courseblockax}}

bedGraph 为四列制表符分隔数据，使用零基半开坐标及有限的有符号数值。导入仅转换一次坐标；原始端点、值、物理行号和不解释样式的显示表头均可追溯。每文件一条信号，多个 track 定义应拆分。字面序列 ID `track`、`browser` 合法，表头不设置图形样式。重叠行保留，不隐式采用末行覆盖规则。

矩形以 `start - 1` 和 `end` 为单元边界，不在缺失区间插值连线。BigWig 需要可选 **rtracklayer**，其解码器不支持 Windows；该平台使用已有 bedGraph 文件。二进制导入的 GRanges 已是一基闭区间。区段值保留二进制浮点精度，区段数量不代表原始测量次数。完整长度来自文件表头。查询采用无链区间并集，结果裁切至并集，空选区返回零行。缺失区间仍缺失，与存储的零不同；不制造 coverage、identity、support、单位或裁切前观测 ID。各 Track 的映射、数值域、单位和 Guide 仍显式配置。

## 已比对 reads：可选覆盖度桥接 {#aligned-reads}

使用可选 **GenomicAlignments**、**Rsamtools** 和 **GenomicRanges** 读取已有 BAM，显式选择过滤与计数口径，再由 ggchord 绘制信号。本例使用依赖包中未改动的果蝇 RNA-seq 子集，原 BAM header 保留 TopHat、参考及 SRR 来源；固定源 URL 和校验和见 `tools/long-term-audit/v0.27.0-read-coverage-results.json`。不运行 read 比对或变异检测。

{{code:courseblockay}}

每条保留的**alignment 记录**在其已比对参考区块上贡献一次。跳过区（`N`）以及本例显式排除的删除区（`D`）贡献零；保留 supplementary alignment。不设碱基质量阈值，不去除重叠配对端或重复分子。这不是自动得到的 fragment、唯一分子、表达量、Identity 或变异支持度。零段是该明确计算中的已存储零，与实验信号缺失不同。大型输入可用依赖包的流式 `BamFile` coverage 方法，并沿用相同过滤声明。

将这些一基闭合区间交给 `geom_track(geom = "rect")`，设置 `xmin = start - 1`、`xmax = end`、`ymin = 0`、`ymax = value`，方法见[连续信号](#signals)。完整参考长度取 `scanBamHeader(bam_file)[[1]]$targets`，不取最后 read 端点。绘制小区域前显式选择源窗口。

CRAM 应由已有 HTSlib/pysam 或 samtools 环境解码，需要参考时提供**实际匹配的参考**，保留原 header、过滤和转换日志。不伪造参考、不自动下载未请求资源，也不宣称所有 CRAM 都可独立解码。定向原型仅用 `no_ref=1` 将相同真实 BAM 记录重编码，无外部参考访问，再核对解码后的全部记录字段/标签、参考字典及覆盖度。这是格式往返，不是原始 CRAM 观测，也不是所有参考型 CRAM 的认证。上游 `SO:sorted` header 和临时索引缺失产生的解码消息均保留；顺序读取成功，没有改写它们。
