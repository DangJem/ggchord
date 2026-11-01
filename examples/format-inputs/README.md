# Source-backed format inputs / 有来源的格式输入

- `pUC19c-L09137.2.gb`: unchanged NCBI GenBank record L09137.2, retrieved through
  EFetch. Its definition is **pUC19c**. This is not silently renamed to the
  package's pUC19 fixture or its historical `plasmid_example_pUC19c` alias.
- EMALE inputs: original data records from Thomas Hackl's gggenomes tutorial
  archive, `https://github.com/thackl/gggenomes/raw/master/data-raw/emales.tgz`.
  The upstream package declares MIT + LICENSE (copyright 2018 Thomas Hackl).
  Dataset context: <https://thackl.github.io/gggenomes/articles/emales.html>,
  endogenous mavirus-like elements in *Cafeteria burkhardae* (Hackl et al. 2021).
  No upstream figure or plotting code is copied.
- `emales-pair.paf`: retains BVI_023A → Cflag_017A rows and optional tags.
- `emales-mavirus-blastp.tsv`: complete unchanged BLASTP 2.9.0+ output, with
  standard 12 fields; weak hits remain in the original and are filtered explicitly
  in tutorial code. Protein coordinates are amino-acid positions.
- `emales-BVI_023A.gff3`: original comments and BVI_023A annotations retained;
  CDS IDs match the query protein IDs. Prodigal prediction is source annotation,
  not experimental gene-function validation.
- `emales-BVI_023A.fna`: unchanged single-genome archive member.

`sources.json` records URLs, retrieval date, archive/member/file SHA256 and
physical source-line ranges (inclusive). To repeat selection from the downloaded
archive: `python3 data-raw/prepare-format-inputs.py /path/to/emales.tgz`.

这些文件保留真实记录，不为演示捏造 identity、长度、比对或链方向。原始文件与教程
显示筛选分开；教程代码只将必要列传入绘图函数。GenBank 解析仍使用已安装的可选
依赖；没有依赖时明确标记跳过，不自动安装。
