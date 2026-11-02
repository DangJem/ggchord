<p align="center">
  <a href="https://dangjem.github.io/ggchord/zh/index.html"><img src="man/figures/logo.png" alt="ggchord" width="320"></a>
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# ggchord：弦图、圆图与共线性图

基于 ggplot2 的 R 包。用独立图层组合序列骨架、区间连接、基因、通用注释、标签和定量轨道。

本说明对应 **ggchord 0.16.5**。从旧版本迁移时，请先阅读迁移指南；正式发布资产见 GitHub Releases。

## 安装

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages(c("remotes", "dplyr"))
remotes::install_github("DangJem/ggchord@v0.16.5")
```

如已检出本版源码，在仓库根目录使用 `pkgload::load_all(".")`。不声称存在 CRAN 发布。 本例用 dplyr 和原生 `|>` 管道整理数据；dplyr 是示例依赖。

## 三种模式

下列代码先从随包数据中选出必需列，再绘图；原始数据与来源信息保留在原对象中。

<details>
<summary>先运行：弦图与共线性图共用的数据和样式</summary>

```r
library(ggchord)
library(dplyr)
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)
seqs <- seq_data_example |>
  select(accver, length)
links <- ribbon_data_example |>
  select(qaccver, saccver, qstart, qend, sstart, send)
# The fixture preserves biological strands separately from its old display strands.
genes <- gene_data_example |>
  transmute(accver, start, end, strand = source_strand)
```

```r
# Optional appearance choices; none of these columns or scales is required input.
library(ggchord)
seq_palette <- setNames(
  c("#3867D6", "#E58A28", "#219B83", "#9B59B6"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1")
)
gene_palette <- c("+" = "#3B6FD8", "-" = "#E1812C", "." = "#8291A6")
```

```r
library(ggchord)
doc_theme <- theme_ggchord_publication(
  axis.text = element_text(size = 11),
  axis.text.check.overlap = TRUE
) + theme(legend.position = "bottom")
doc_axis <- function(coordinates = "relative") {
  scale_seq_position_continuous(
    coordinates = coordinates, n.breaks = 3, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
  )
}
```

</details>

### 弦图：查看序列间的总体连接

展示区间之间的对应关系。

```r
p_chord <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .7) +
  scale_ribbon_fill_manual(values = seq_palette, guide = "none") +
  geom_seq(show.legend = FALSE) +
  geom_gene(position = position_strand(.1), width = .04) +
  geom_seq_label(size = 4, seq_label_radius = .84, seq_label_hjust = .5) +
  coord_chord() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  doc_theme
p_chord
```

![弦图：查看序列间的总体连接](man/figures/chord.png)

[查看完整教程 →](https://dangjem.github.io/ggchord/zh/real-data.html)

### 圆图：读取单条序列的注释

先准备 pUC19 注释，再组合骨架与特征层。

<details>
<summary>先运行：圆图的数据和配色</summary>

```r
library(dplyr)
data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19 |>
  select(accver, length)
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source |>
  select(accver, start, end, strand, type, anno)
```

```r
library(ggchord)
feature_palette <- c(
  CDS = "#3973CB", promoter = "#8B5EC4", rep_origin = "#E8A42D",
  primer_bind = "#168F86", misc_feature = "#A0ADBD", protein_bind = "#CC6485"
)
```

</details>

```r
circle_features <- features |>
  filter(type %in% c("CDS", "rep_origin"))
p_circular <- ggchord(seq_data = plasmid_seq) +
  geom_seq() +
  geom_feature(aes(feature_fill = type), data = circle_features, position = position_plasmid()) +
  geom_feature_label_repel(data = circle_features, position = position_plasmid(), size = 4) +
  geom_seq_center_label() +
  scale_feature_fill_manual(values = feature_palette, breaks = unique(circle_features$type)) +
  coord_circular() +
  theme_ggchord_publication(axis.text = element_text(size = 11))
p_circular
```

![圆图：读取单条序列的注释](man/figures/circular.png)

[查看完整教程 →](https://dangjem.github.io/ggchord/zh/plasmids.html)

### 共线性图：对齐并比较位置

复用第一节的数据和样式，按明确选择的基因锚点对齐。

```r
anchors <- genes |>
  group_by(accver) |>
  slice_max(end - start, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(accver, start, end)
p_collinear <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .7) +
  scale_ribbon_fill_manual(values = seq_palette, guide = "none") +
  geom_seq() +
  geom_gene(position = position_strand(.1)) +
  geom_seq_label(size = 4) +
  coord_collinear(anchor = anchors, align = "center", orientation = "auto") +
  doc_axis() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) +
  doc_theme
p_collinear
```

![共线性图：对齐并比较位置](man/figures/collinear.png)

[查看完整教程 →](https://dangjem.github.io/ggchord/zh/modes.html)

`aes()` 映射数据，`scale_*()` 定义映射和图例，geom 定义几何与固定样式，Position 负责局部放置，coord 负责全局坐标，theme 负责非数据外观。

## 预览与保存

RStudio Plots 会按面板尺寸重新排列标签。按相同显式宽高预览并保存：

```r
view_ggchord(p_collinear, width = 8, height = 5, units = "in")
ggsave("comparison.png", p_collinear, width = 8, height = 5, units = "in", dpi = 144)
```

## 继续阅读

- [开始使用](https://dangjem.github.io/ggchord/zh/installation.html)
- [按任务查找教程](https://dangjem.github.io/ggchord/zh/articles.html)
- [函数文档](https://dangjem.github.io/ggchord/zh/reference/index.html)
- [示例图集](https://dangjem.github.io/ggchord/zh/gallery.html)
- [迁移指南](https://dangjem.github.io/ggchord/zh/migration.html)
- [更新日志](https://dangjem.github.io/ggchord/zh/news.html)

[Issues](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
