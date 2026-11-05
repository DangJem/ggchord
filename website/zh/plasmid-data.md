单序列圆图可从一张两列表开始。要识别特征、引物或酶切位点，则另外提供真实DNA；已有注释时可以直接画表格，无需再次识别。

## 输入表与可选信息 {#tables}
| 输入 | 必需内容 | 可选内容 |
| --- | --- | --- |
| 骨架 | accver, length，且只有一条序列 | label用于显示名称；sequence用于搜索 |
| 通用特征 | accver, start, end与方向声明 | type、anno/label、feature_id、形状类别 |
| 引物 | 已计算的结合区间、链与名称等对应geom要求 | 识别与来源信息 |
| 酶切位点 | 搜索结果中对应geom的位点与名称 | 酶目录、甲基化假设与筛选来源 |

最简单骨架只需要前两列，不能要求所有人都提供DNA。后两类优先使用本包搜索返回的表，字段细节查相应函数帮助。未知链特征由通用feature层表达，不伪造正链。

## 查看随包真实记录 {#inspect}
```r
library(dplyr)
plasmid_example_pUC19 |> select(accver, length)
find_common_features(plasmid_example_pUC19) |> select(accver, start, end, strand, type, anno) |> slice_head(n = 5)
```

单序列 ID 列叫 `accver`，可以使用自定义名称。完整DNA与长度必须一致。识别得到的是计算匹配，不是实验验证。

## 导入自己的质粒 {#import}
`read_fasta_lengths()`只读取FASTA标识与完整长度。识别还需要额外读取实际DNA，例如用可选Biostrings的 `readDNAStringSet()` 和 `as.character()`，保留与表格一致的ID；有GFF3/BED/GenBank或GRanges时按[格式专题](formats.html)导入注释，保留来源。跨圆周起点与分段记录按对应导入契约处理，不手工把所有片段合成一个虚假的连续区间。

## 案例如何循序渐进 {#path}

在[添加特征和标签](plasmid-build.html)逐步构建 pUC19，然后做[pUC19 与 pETDuet-1 完整案例](plasmids.html)。按需查阅[搜索条件](recognition.html)和[13 张参考图](plasmid-catalog.html)。
