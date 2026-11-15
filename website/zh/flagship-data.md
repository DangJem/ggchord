旗舰绘图页直接使用本页列出的随包数据。安装 ggchord 后即可加载；只有重新构建数据时才需要源码中的原始文件。完整来源表与展示子集分开保留。

## 数据入口 {#datasets}

| 数据集 | 对应旗舰 | 可直接使用的内容 |
| --- | --- | --- |
| `flagship_comparison_example` | Chord、Collinear、Radial、Sankey | 每个视图自己的序列表、已选择 HSP、CDS、标签顺序和显示方向；Collinear 另含固定布局的路由控制点 |
| `flagship_dotplot_example` | 四序列 Dotplot | 全部32条 HSP、完整序列长度、16个单元格的计数标签及固定轴间隙 |
| `flagship_fusion_example` | 六面板 Fusion | 总览 HSP/CDS、跨面板关系、源窗口、局部表、60个真实编码碱基及其源位置、放大区间 |
| `flagship_circular_example` | 密集圆图与13个完整参考圆图 | 每图的完整序列长度、特征、独立参考引物、实际采用的酶集合、候选位点数和集合选择说明 |

{{code:courseblockg}}

```r
data(flagship_comparison_example)
chord_data <- flagship_comparison_example$chord
names(chord_data)
```

`flagship_comparison_example` 的四个旗舰入口是 `chord`、`collinear`、`radial`、`sankey`。`flagship_circular_example$dense` 对应既有密集旗舰；`flagship_circular_example$references` 按原始序列标识命名。另有 `comparison` 入口，保留四坐标严格对照的原始教学子集。加载这些入口无需再运行识别、筛选、聚焦或字符提取。

## 原始来源与单位 {#sources}

四基因组来源为 MT108731.1、MT118296.1、OQ646790.1、OR222515.1 的原始 FASTA、GFF3 和标准 BLAST 表，保留在源码的 `examples/fasta/`、`examples/gff3/`、`examples/blastn-standard/`。完整输入仍可从 `genome_comparison_example` 追溯：4条完整序列、32条已观测 HSP、266条 CDS。原始输入与上游记录见[真实输入流程](real-data.md)。

圆图来自随包的13条完整质粒序列；识别数据库、参考引物和酶集合的来源及许可见源码 `inst/NOTICE` 和相应 `data-raw/` 生成规则。保留原始标识、序列字符、坐标、链和记录顺序，不为改善视觉效果修改生物学数据。

核酸长度与坐标单位为 bp，区间采用一基闭区间，Identity 为百分比。反向 HSP 的端点可降序，不能为了消除提示排序端点。显示 orientation 与生物学 strand 分别保存。

## 比较图的展示选择 {#comparison}

Chord 和 Radial 使用同一12条 HSP及20条 CDS：HSP 长度至少900 bp、Identity 至少75%，每个来源序列对最多取长度最大的3条；CDS 按真实链和源位置分组抽样，正链优先保留与显示端点相交的记录。原始短 CDS 和真实链保持不变。

Collinear 和 Sankey 保留16条 CDS及10条 HSP。两个拥挤候选只从这组展示中省略，仍保留在完整来源表中；相邻行保留选中 HSP，非相邻行每对仅选最长的一条。标签输入顺序也随数据保存。

显示方向和 Collinear 路由控制点是既有旗舰的明确预设，适用于原顺序、间隙及对齐配置。改变这些几何设置时，应重新计算相应预设。它们不是生物学注释，也不是设备测量或标签避让结果。

## Dotplot 与 Fusion 的派生表 {#composition}

Dotplot 保留全部32条 HSP，不补造反向、自比对或缺失记录。计数表包含四序列的16个组合，零值仅表示本次输入无记录。标签位置按原生序列长度和4,000 bp轴间隙构建；改变轴拼接方式时应重新构建标签位置。

Fusion 每个来源序列对选择最接近4,000 bp的一条真实 HSP，每条序列按来源注释顺序选择3条 CDS。局部窗口来自已选择关系的真实端点，两侧各扩展800 bp并裁切到源边界；完整来源、局部裁剪表和用于上下文的特征分别保存。

字符面板使用局部下方序列第二条显示 CDS 的前60个真实编码碱基；负链按编码方向取反向互补。每个字符保留真实源位置，注册字符图及放大连接时直接使用这些字段，不从宏观 HSP 猜测逐碱基对齐。

## 圆图的识别与酶集合 {#circular}

特征识别、参考引物匹配和酶切搜索只在数据构建时执行。密集旗舰沿用商业全集中的唯一六碱基及以上集合。完整参考图优先使用保存的参考酶集合；参考集合明确为 None 时，沿用商业非冗余集合中的唯一六碱基及以上位点，并在每条记录的 `policy` 中说明。

`sites` 是实际绘制的完整集合，`all_site_count` 记录未筛选搜索结果的行数。大量未采用的候选不重复随包保存，可从原序列及冻结的识别数据库重建。旧三切点图库的抽样另存为 `gallery_sites`，与完整参考图分开使用，不能把三切点抽样说成完整酶集合。识别和集合选择均保留现有规则。

例如，密集圆图的位点集合可这样重建；搜索结果仅用于构建，不是实验切割测量。

```r
all_sites <- find_restriction_sites(plasmid_example_pETDuet_1)
dense_sites <- filter_restriction_sites(
  all_sites, set = "unique_6plus", parent_set = "commercial_all"
)
dense_sites |>
  select(accver, enzyme, position) |>
  head()
```

## 重建与核对 {#rebuild}

在源码根目录、安装 pkgload、dplyr（≥1.1.0）、purrr 和 ggplot2 后执行：

```sh
Rscript data-raw/generate_website_data.R
Rscript data-raw/generate_flagship_data.R
```

完整代码见[来源表与教学数据生成脚本](../../data-raw/generate_website_data.R)和[旗舰数据生成脚本](../../data-raw/generate_flagship_data.R)。原始输入校验和见[来源文件清单](../../data-raw/website-data-sources.tsv)，旗舰生成依赖校验和见[旗舰输入清单](../../data-raw/flagship-data-sources.tsv)。

随包对象只保存来源/展示表和明确的固定预设；不保存 ggplot 对象、grob、设备测量或最终标签位置。绘图时继续按真实物理尺寸计算几何与标签避让。重建后应核对来源记录、展示顺序及成图效果，再更新随包数据。

## 返回绘图 {#plots}

[四种比较图](multi-sequence-views.md#flagship) · [Dotplot](dotplot-start.md#matrix) · [Fusion](fusion-layout.md#flagship) · [密集圆图](circular-dense-case.md) · [完整参考圆图](plasmid-catalog.md)
