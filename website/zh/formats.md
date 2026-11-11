<div class="task-intro"><p>先根据手中的数据选择入口。读入文件后，核对返回的 ID、长度和端点，再交给绘图层。</p></div>

## 手中的文件从哪里开始 {#choose-format}

| 手中的文件 | 下一步 |
| --- | --- |
| FASTA 序列 + GFF3 注释 + BLAST 局部比对 | [从源文件画比较图](real-data.html#read) |
| BED 区间或 PAF（Pairwise mApping Format）比对 | [核对零基坐标](#zero-based) |
| GenBank 单序列及注释 | [保留原始记录与分段](#genbank) |
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
