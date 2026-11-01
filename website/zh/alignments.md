## 默认表格列就够用 {#standard}
不需要早期教程中的 17 列自定义字段。使用 BLAST 标准表格输出即可：

```sh
blastn -query query.fna -subject subject.fna -outfmt 6 -out result.tsv
blastp -query proteins.faa -subject reference.faa -outfmt 7 -out result.tsv
```

`6` 是无注释表格，`7` 加上注释头；不再列出自定义字段时，两者使用默认 12 列。
完全省略 `-outfmt` 得到的是 pairwise 文本，当前表格读取器不解析该格式。
默认列为 query ID、subject ID、pident、length、mismatch、gapopen、qstart、qend、
sstart、send、evalue、bitscore。读取器将两个 ID 的角色命名为 `qaccver`、`saccver`，
保留字符串原值，并不要求它们一定是 accession version。
[NCBI 格式说明](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/table/appendices.T.options_common_to_all_blast/)。

## 文件布局与绘图最小列是两回事 {#minimal}
文件里有 12 列，不表示绘图必须接收 12 列。先从[共享准备代码](index.html#quick-start)
加载包，然后运行真实记录的例子：

{{code:standardblast}}

`blast_source` 保留全部统计和来源信息；`blast_links` 只有六个端点角色。
`pident` 只有在按相似度配色/筛选时才需要，`length` 只有按比对长度筛选/加权时才需要。
`qlen`、`slen` 用于相应覆盖率计算，不能为了基础出图强制要求它们。
无表头的自定义列须明确 `col_names`，读取器不猜测任意六列的含义。
也可对任意已读入的表使用 `as_ggchord_link(x, mapping = ...)`，见[格式转换](formats.html#adapters)。

仓库 `examples/blastn-standard/` 的文件是历史原始 HSP 的前 12 列转换，保留行和值，
不是本轮重新计算的 BLAST 结果。旧 `.o7` 文件保持原样，转换脚本和 hash 在来源清单中。

## BLASTp 可以读取，但坐标单位不同 {#protein}
下面是真实 BLASTP 2.9.0+ 结果与对应查询 CDS 注释，来自公开 EMALE 教程。
先保留完整输入，再选择有完整 CDS 映射、足够长度和较小 evalue 的查询 HSP。

{{code:protein}}
{{results:protein}}

BLASTp 的 `qstart/qend/sstart/send` 是氨基酸坐标。这里将查询蛋白的 HSP 范围
转换到已知 CDS 的基因组坐标，保留负链和 phase；有缺口的 HSP 仍可投影查询范围，
但这只是范围，不能解释为逐个核苷酸的配对。原始弱命中、gapopen 和统计仍在原表中。

要得到两端的核苷酸连接，还须提供 subject 的 CDS 映射。不能从蛋白 accession、
最大端点或同名注释猜出来。`project_ggchord_links()` 的 `ungapped = TRUE` 必须有
证据，不能为了让默认 BLASTp 文件通过而随意打开。[投影教程](projection.html)
使用两条编码 DNA 完全一致的真实质粒 CDS，展示可证明无缺口的双端例子。

本版不直接绘制以 aa 为轴的蛋白骨架。`blastx`、`tblastn` 的两侧单位也可能不同，
必须逐侧核对和转换；读取成功不代表可以直接用 bp 坐标绘图。

## 来源与限制 {#source}
EMALE 原始输出来自 [gggenomes 的公开数据与教程](https://thackl.github.io/gggenomes/articles/emales.html)。
本项目只复用数据记录，用 ggchord 自行绘图；没有复制上游图或绘图代码。
下载、原始成员、行号范围、SHA256 和筛选规则见教学输入包的 `sources.json`。
表格解析不重新运行比对，也不验证命中的生物学功能。
