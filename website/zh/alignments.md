<div class="task-intro"><p>用真实比对记录学习 BLAST 的默认列，并选出绘图所需的六个端点字段。</p></div>

## 为什么从 BLAST 输入开始 {#prepare}

ggchord 最初为 BLASTn 比对结果设计，用弦图展示核酸序列之间的对应区间。因此，BLAST 表格是 ggchord 的默认推荐比对输入，连接表也延续 query/subject 的字段风格。随着绘图能力扩展，这套命名继续用于弦图、共线性图及数据转换。

绘图需要完整序列的长度和已计算的比对结果。FASTA 提供骨架长度，BLAST 提供区间对应关系；ggchord 将这些记录组织成图。下面先读懂默认 BLAST 表格，再分别绘制核酸与蛋白比对。

`-outfmt 6` 输出无表头的制表符分隔数据，`-outfmt 7` 在相同数据行之外添加以 `#` 开头的注释。省略 `-outfmt` 时，BLAST 输出供阅读的配对文本。两个命令分别演示核酸与蛋白比对。

```sh
blastn -query query.fna -subject subject.fna -outfmt 6 -out hits-dna.tsv
blastp -query proteins.faa -subject reference.faa -outfmt 7 -out hits-protein.tsv
```

## 默认表格列就够用 {#standard}

不指定字段列表时，两种表格格式都按以下顺序输出 12 列。

| 顺序 | 字段 | 读法 |
| --- | --- | --- |
| 1–2 | `qseqid`、`sseqid` | query 与 subject 的序列 ID，读入后名为 `qaccver`、`saccver`，值原样保留 |
| 3–6 | `pident`、`length`、`mismatch`、`gapopen` | 一致性百分比、比对列数、错配数、缺口开启次数 |
| 7–10 | `qstart`、`qend`、`sstart`、`send` | 两端的起止坐标，从 1 开始并包含两端。递减端点保留比对方向 |
| 11–12 | `evalue`、`bitscore` | 期望值和以 bit 为单位的得分，通常前者越小、后者越大，命中证据越强 |

解压[教学输入](../assets/tutorial-inputs.zip)后，从包含 `examples/` 的目录运行。示例文件保存了实际 HSP 的前 12 列。

{{code:standardblast}}
{{output:standardblast}}

读取器把绘图所需的六个端点字段排在表格前面。

## 选出绘图所需的六列 {#minimal}

{{code:blastendpoints}}
{{output:blastendpoints}}

一条带由两条序列的 ID 和各自的起止坐标确定。要按一致性着色或筛选时保留 `pident`，要按比对长度筛选时保留 `length`。用自定义字段输出无表头文件时，在 `read_blast()` 中用 `col_names` 按实际顺序声明列。已经读入的其他表可用[连接格式转换](formats.html#adapters)。

## 先选择坐标单位 {#protein}

| 想比较的对象 | 原生单位 | 比对输入 | 下一步 |
| --- | --- | --- | --- |
| 核酸序列 | bp | BLASTn | 将 FASTA 全长和连接表交给绘图层 |
| 蛋白序列 | aa | BLASTp | 用完整蛋白全长声明 sequence_unit = "aa" |
| 蛋白命中在基因组上的位置 | bp | BLASTp + CDS 模型 | 将蛋白区间投影到核酸坐标 |

一张图使用一个坐标单位。蛋白长度来自完整序列，命中端点仅描述参与比对的区间。BLASTx、tBLASTn 的两端单位不同，需要先选定统一的显示坐标。基因的正负链与蛋白的 N→C 顺序是不同的生物学含义。

[查看 BLASTp 原生蛋白比较图](protein-alignments.html#protein)，或继续[了解三种绘图布局](modes.html)。

<span id="projection"></span>要将蛋白命中定位到核酸基因组，可阅读[CDS 坐标投影](projection.html)。

## 来源与限制 {#source}
EMALE 原始输出来自 gggenomes 的公开 EMALE 数据。
本项目只复用数据记录，用 ggchord 自行绘图；没有复制上游图或绘图代码。
下载、原始成员、行号范围、SHA256 和筛选规则见教学输入包的 `sources.json`。
表格解析不重新运行比对，也不验证命中的生物学功能。

## 完成后检查 {#check-result}
先确认本次比较的对象是核酸还是蛋白，再检查表格列、完整长度与坐标单位。

<div class="next-step">继续：<a href="projection.html">将蛋白坐标转换为基因组坐标 →</a></div>
