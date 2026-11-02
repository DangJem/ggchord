<p align="center">
  <a href="https://dangjem.github.io/ggchord/zh/index.html"><img src="man/figures/logo.png" alt="ggchord" width="320"></a>
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# ggchord：弦图、共线性图与质粒图

基于 ggplot2 的 R 包。用独立图层组合序列骨架、区间连接、基因、通用注释、标签和定量轨道。

本版为 **ggchord 0.16.6**；以下展示由本版代码生成。旧代码迁移见迁移指南，正式资产见 GitHub Releases。

## 安装

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages(c("remotes", "dplyr", "purrr", "tidyr"))
remotes::install_github("DangJem/ggchord@v0.16.6")
```

如已检出本版源码，在仓库根目录使用 `pkgload::load_all(".")`。不声称存在 CRAN 发布。 示例使用 dplyr、purrr、tidyr 和原生 `|>`，这些是示例可选依赖。

## 三种模式

以下三张综合示例展示本版的组合能力。原始输入保持可追溯，完整准备代码和解释见对应教程；初次使用请从最小首图开始。

<details>
<summary>先运行：三张图共用的样式</summary>

```r
# Optional appearance choices; none of these columns or scales is required input.
library(ggchord)
seq_palette <- setNames(
  c("#E64B35", "#233E4A", "#269B82", "#D59B17"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1")
)
gene_palette <- c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")
```

```r
library(ggchord)
doc_background <- "#FFFCF5"
identity_palette <- c("#303B71", "#36537D", "#269B82", "#A5D52C", "#ECEB31")
identity_values <- c(0, 0.5, 0.8, 0.9, 1)
identity_breaks <- c(0, 50, 80, 90, 100)
doc_legend_theme <- theme(
  legend.title = element_text(size = 13, face = "bold"),
  legend.text = element_text(size = 11),
  legend.key.width = unit(9, "mm"), legend.key.height = unit(5, "mm"),
  legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
  legend.margin = margin(l = 4, r = 4, unit = "mm")
)
sequence_guide <- guide_ggchord_legend(nrow = 1, theme = doc_legend_theme)
strand_guide <- guide_ggchord_legend(
  nrow = 1, override.aes = list(strand = c("+", "-")), theme = doc_legend_theme
)
identity_guide <- guide_ggchord_colourbar(
  position = "left", direction = "vertical",
  theme = theme(
    legend.key.height = unit(75, "mm"), legend.key.width = unit(5, "mm"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11), legend.margin = margin(r = 7, unit = "mm")
  )
)
doc_theme <- theme_ggchord_publication(
  axis.text = element_text(size = 11),
  axis.text.check.overlap = TRUE
) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(9, "mm"), legend.key.height = unit(5, "mm"),
    legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
    legend.box.spacing = unit(6, "mm"),
    plot.background = element_rect(fill = doc_background, colour = NA),
    panel.background = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 10, colour = "#5C625A"),
    plot.margin = margin(14, 18, 14, 18)
  )
doc_axis <- function(coordinates = "relative") {
  scale_seq_position_continuous(
    coordinates = coordinates, n.breaks = 3, minor_breaks = NULL,
    labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE)
  )
}
```

</details>

### 弦图：四基因组、混合曲率与相似度

直线与弧形骨架、真实链方向基因、精选标签及筛选后的观测比对。

![弦图：四基因组、混合曲率与相似度](man/figures/showcasechord.png)

<details>
<summary>数据准备与绘图代码 · 展开运行</summary>

```r
library(ggchord)
library(dplyr)
library(purrr)
library(tidyr)
# Original repository records; extract the tutorial input bundle first.
show_sequences <- read_fasta_lengths(files = "examples/fasta/*.fna") |>
  select(accver, length)
show_annotations <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep")
show_genes <- show_annotations |>
  filter(type == "CDS") |>
  select(accver, start, end, strand, anno)
show_alignments <- read_blast(files = "examples/blastn-standard/*.tsv")
# At most three substantial observed HSPs per pair; no coordinates are invented.
show_links <- show_alignments |>
  filter(length >= 900, pident >= 75) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 3, with_ties = FALSE) |>
  ungroup()
# Select up to four source CDS per sequence, spread across annotation order.
# No feature is moved or invented; keep complete source tables separately.
show_features <- show_genes |>
  filter(end - start + 1 >= 500) |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  mutate(display_bin = ntile(row_number(), 4), informative = !is.na(anno) & anno != "hypothetical protein") |>
  group_by(accver, display_bin) |>
  arrange(desc(informative), desc(end - start), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(accver, start, end, strand, anno)
# Also retain well-separated negative-strand CDS, if not already represented.
show_reverse <- show_genes |>
  filter(strand == "-", end - start + 1 >= 500) |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  slice(c(1, n())) |>
  ungroup()
show_features <- bind_rows(show_features, show_reverse) |>
  distinct(accver, start, end, strand, .keep_all = TRUE) |>
  arrange(accver, start)
show_labels <- show_features
show_optimization <- optimize_ggchord_layout(
  show_sequences, show_links,
  optimize = "orientation", weight = "length"
)
```

```r
show_curvature <- c(
  "MT108731.1" = 0, "MT118296.1" = 0.85,
  "OQ646790.1" = 1.1, "OR222515.1" = -0.4
)
show_radius <- c(
  "MT108731.1" = 1, "MT118296.1" = 1.05,
  "OQ646790.1" = 0.95, "OR222515.1" = 1.08
)
p_show_chord <- ggchord(show_sequences, show_links, show_features) +
  geom_link_ribbon(aes(ribbon_fill = pident), link_avoid = "smooth", alpha = 0.58, linewidth = 0.18) +
  geom_seq(
    seq_orientation = show_optimization$seq_orientation,
    seq_curvature = show_curvature, seq_radius = show_radius,
    seq_gap = 0.06, linewidth = 1.1
  ) +
  geom_gene(position = position_strand(0.07), width = 0.035) +
  geom_gene_label_repel(
    data = show_labels, label_layout = "radial", label_wrap = 24,
    size = 3.1, label_max_lines = 3, position = position_strand(0.07)
  ) +
  geom_seq_label(size = 4, seq_label_radius = 0.84, seq_label_hjust = 0.5) +
  scale_seq_colour_manual(
    values = seq_palette, name = "Sequence",
    guide = sequence_guide
  ) +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  scale_ribbon_fill_gradientn(
    colours = identity_palette, values = identity_values,
    limits = c(0, 100), breaks = identity_breaks, name = "Identity (%)", guide = identity_guide
  ) +
  guides(gene_fill = strand_guide) +
  scale_seq_position_continuous(breaks = c(0, 20000, 40000), minor_breaks = NULL) +
  coord_chord(rotation = 50) +
  doc_theme +
  labs(
    title = "Four genomes, one connected view",
    subtitle = "Signed curvature and radius · selected strand-aware CDS · local ribbon clearance"
  )
p_show_chord
```

</details>

[查看完整教程 →](https://dangjem.github.io/ggchord/zh/chord.html)

### 共线性图：四序列与跨行连接

正负链分列、逐一标注的精选CDS、局部避障、虚线跨行ribbon与右侧控制点路由。

![共线性图：四序列与跨行连接](man/figures/showcasecollinear.png)

<details>
<summary>数据准备与绘图代码 · 展开运行</summary>

```r
library(ggchord)
library(dplyr)
library(purrr)
library(tidyr)
# Original repository records; extract the tutorial input bundle first.
show_sequences <- read_fasta_lengths(files = "examples/fasta/*.fna") |>
  select(accver, length)
show_annotations <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep")
show_genes <- show_annotations |>
  filter(type == "CDS") |>
  select(accver, start, end, strand, anno)
show_alignments <- read_blast(files = "examples/blastn-standard/*.tsv")
# At most three substantial observed HSPs per pair; no coordinates are invented.
show_links <- show_alignments |>
  filter(length >= 900, pident >= 75) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 3, with_ties = FALSE) |>
  ungroup()
# Select up to four source CDS per sequence, spread across annotation order.
# No feature is moved or invented; keep complete source tables separately.
show_features <- show_genes |>
  filter(end - start + 1 >= 500) |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  mutate(display_bin = ntile(row_number(), 4), informative = !is.na(anno) & anno != "hypothetical protein") |>
  group_by(accver, display_bin) |>
  arrange(desc(informative), desc(end - start), .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(accver, start, end, strand, anno)
# Also retain well-separated negative-strand CDS, if not already represented.
show_reverse <- show_genes |>
  filter(strand == "-", end - start + 1 >= 500) |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  slice(c(1, n())) |>
  ungroup()
show_features <- bind_rows(show_features, show_reverse) |>
  distinct(accver, start, end, strand, .keep_all = TRUE) |>
  arrange(accver, start)
show_labels <- show_features
show_optimization <- optimize_ggchord_layout(
  show_sequences, show_links,
  optimize = "orientation", weight = "length"
)
```

```r
# Omit two short, crowded candidates from this linear display only.
# Both source rows remain in show_features and in the chord display.
show_linear_features <- show_features |>
  filter(!anno %in% c("tail assembly protein", "gp4 protein"))
# Separate adjacent and cross-row links, retaining original endpoint values.
show_pair_rows <- show_links |>
  mutate(
    query_row = match(qaccver, show_sequences$accver),
    subject_row = match(saccver, show_sequences$accver),
    span = abs(query_row - subject_row)
  )
show_adjacent <- show_pair_rows |> filter(span == 1)
# One longest observed HSP per non-neighbour pair keeps the long routes legible.
show_cross <- show_pair_rows |>
  filter(span > 1) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 1, with_ties = FALSE) |>
  ungroup()
# Inspect a backbone-only layout using the same alignment and gap as the final view.
show_linear_coord <- coord_collinear(align = "center", orientation = "auto", ratio = 0.7)
show_routing_base <- ggchord(show_sequences, show_links) +
  geom_seq(seq_gap = 0.42) + show_linear_coord
show_reference <- get_chord_layout(show_routing_base)$sequence_reference
# Shared controls sit just beyond the intervening backbones, rather than far away.
show_cross_controls <- map2(show_cross$query_row, show_cross$subject_row, function(query_row, subject_row) {
  between <- seq(min(query_row, subject_row) + 1L, max(query_row, subject_row) - 1L)
  middle_ids <- show_sequences$accver[between]
  right_edge <- max(show_reference$x_offset[middle_ids] + show_reference$lens[middle_ids] * show_reference$x_scale)
  middle_y <- mean(show_reference$y[show_sequences$accver[c(query_row, subject_row)]])
  c(right_edge + if (abs(query_row - subject_row) == 3) 0.85 else 0.65, middle_y)
})
p_show_collinear <- ggchord(show_sequences, gene_data = show_linear_features) +
  geom_link_ribbon(
    aes(ribbon_fill = pident),
    data = show_adjacent,
    link_avoid = "smooth", alpha = 0.48, linewidth = 0.22
  ) +
  geom_link_ribbon(
    aes(ribbon_fill = pident),
    data = show_cross,
    adjacent_only = FALSE, link_ctrl_point = show_cross_controls,
    link_avoid = "smooth", alpha = 0.28, linewidth = 0.3
  ) +
  geom_seq(seq_gap = 0.42, linewidth = 1) +
  geom_gene(position = position_strand(0.09), width = 0.045) +
  geom_gene_label_repel(
    data = show_linear_features |> arrange(desc(nchar(anno)), accver, start),
    label_layout = "auto", label_wrap = 24,
    label_side = "auto", position = position_strand(0.09), size = 3.3
  ) +
  geom_seq_label(size = 4) +
  geom_seq_scale(length = 10000, offset = unit(c(0, 10), "mm")) +
  scale_seq_colour_manual(
    values = seq_palette, name = "Sequence",
    guide = sequence_guide
  ) +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  scale_ribbon_fill_gradientn(
    colours = identity_palette, values = identity_values,
    limits = c(0, 100), breaks = identity_breaks, name = "Identity (%)", guide = identity_guide
  ) +
  guides(gene_fill = strand_guide) +
  scale_seq_position_continuous(breaks = c(0, 20000, 40000, 60000, 80000), minor_breaks = NULL) +
  show_linear_coord +
  doc_theme +
  labs(
    title = "Sequence architecture, across four genomes",
    subtitle = "Selected labelled CDS · strand-separated arrows · routed cross-row comparisons"
  )
p_show_collinear
```

</details>

[查看完整教程 →](https://dangjem.github.io/ggchord/zh/collinear.html)

### 质粒图：pETDuet-1密集注释

双链骨架、20个特征、5个参考引物和58个合并切点；位置和名称来自真实输入。

![质粒图：pETDuet-1密集注释](man/figures/showcaseplasmid.png)

<details>
<summary>数据准备与绘图代码 · 展开运行</summary>

```r
library(ggchord)
library(dplyr)
data(plasmid_example_pETDuet_1)
show_plasmid <- plasmid_example_pETDuet_1 |>
  select(accver, length)
show_plasmid_features <- find_common_features(plasmid_example_pETDuet_1)
show_plasmid_primers <- find_primer_bindings(plasmid_example_pETDuet_1, set = "reference")
show_plasmid_all_sites <- find_restriction_sites(plasmid_example_pETDuet_1)
show_plasmid_sites <- filter_restriction_sites(
  show_plasmid_all_sites,
  set = "unique_6plus", parent_set = "commercial_all"
)
```

```r
show_plasmid_position <- position_feature_stack(
  spacing = 0.085, base_position = position_plasmid()
)
p_show_plasmid <- ggchord(show_plasmid, validate = "none") +
  geom_seq(seq_style = "double", colour = "#252525") +
  geom_feature_plasmid(data = show_plasmid_features, position = show_plasmid_position) +
  geom_feature_label_repel(
    data = show_plasmid_features, position = show_plasmid_position,
    external = TRUE, max_overlaps = 0,
    label_fit = "shrink", label_min_size = 2.5
  ) +
  geom_primer(data = show_plasmid_primers, colour = "#A020F0") +
  geom_primer_label_repel(
    data = show_plasmid_primers, colour = "#A020F0", max_overlaps = 0
  ) +
  geom_restriction_site(data = show_plasmid_sites, colour = "#252525") +
  geom_seq_center_label() +
  coord_circular(
    rotation = 90, fit = "manual",
    xlim = c(-3.4, 4.6), ylim = c(-1.6, 2.4)
  ) +
  theme_ggchord_plasmid() +
  theme(
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
    plot.background = element_rect(fill = doc_background, colour = NA),
    panel.background = element_blank()
  )
p_show_plasmid
```

</details>

[查看完整教程 →](https://dangjem.github.io/ggchord/zh/plasmids.html#complex)

`aes()` 映射数据，`scale_*()` 定义映射和图例，geom 定义几何与固定样式，Position 负责局部放置，coord 负责全局坐标，theme 负责非数据外观。

## 预览与保存

RStudio Plots 会按面板尺寸重新排列标签。按相同显式宽高预览并保存：

```r
view_ggchord(p_show_collinear, width = 12, height = 8, units = "in")
ggsave("comparison.png", p_show_collinear, width = 12, height = 8, units = "in", dpi = 144)
```

## 继续阅读

- [开始使用](https://dangjem.github.io/ggchord/zh/installation.html)
- [按任务查找教程](https://dangjem.github.io/ggchord/zh/articles.html)
- [函数文档](https://dangjem.github.io/ggchord/zh/reference/index.html)
- [示例图集](https://dangjem.github.io/ggchord/zh/gallery.html)
- [迁移指南](https://dangjem.github.io/ggchord/zh/migration.html)
- [更新日志](https://dangjem.github.io/ggchord/zh/news.html)

[Issues](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
