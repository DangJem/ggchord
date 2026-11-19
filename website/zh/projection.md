<a id="prepare"></a>

## 选择输入 {#prepare}

绘图前准备一张 CDS 表和一个已知的氨基酸区间。完整 CDS 提供来源坐标、链及阅读框。

<a id="cds"></a>

## 从氨基酸区间找到核酸坐标 {#cds}

随包 pUC19 参考注释给出 lacZα 的完整 CDS。第 2–5 个残基落在原始负链的 455–466 bp。

```r
record <- flagship_circular_example$references[["pUC19"]]
coding <- record$features |>
  dplyr::filter(anno == "lacZα", type == "CDS") |>
  dplyr::transmute(accver, start, end, strand, feature_id = "lacZ_alpha", phase = 0L)
cds <- build_ggchord_cds(coding, genetic_code = 11)
span <- data.frame(feature_id = "display", parent_id = "lacZ_alpha", start = 2, end = 5)
mapped <- project_ggchord_features(span, cds)
mapped[c("accver", "start", "end", "strand")]
```

<a id="boundaries"></a>

## 跨片段与 phase {#boundaries}

多个 CDS 片段可以共享一个父 ID。投影按转录方向排列片段，返回的 segments 列保留内含子间隔。单看外包区间会把这些间隔也计入。

<a id="links"></a>

## 投影两端的已知对应区间 {#links}

下面两条教学坐标轴各含十个密码子。示例明确给定三个残基的一一对应，同时展示正负链的端点方向。真实输入需要自己的对应证据。

```r
models <- data.frame(
  feature_id = c("A", "B"), accver = c("row1", "row2"),
  start = c(101, 301), end = c(130, 330), strand = c("+", "-"), phase = 0L
)
pair <- data.frame(qfeature_id = "A", sfeature_id = "B",
  qstart = 2, qend = 4, sstart = 5, send = 7)
project_ggchord_links(pair, build_ggchord_cds(models), ungapped = TRUE)
```

<a id="check-result"></a>

## 核对投影 {#check-result}

检查父 CDS 标识、原生单位和完整片段列表。实际 BLASTp 命中若包含 gap，先判断已有证据能支持单侧区间还是双端逐位对应。
