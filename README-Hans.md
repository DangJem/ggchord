<p align="center">
  <a href="https://dangjem.github.io/ggchord/zh/index.html"><img src="man/figures/logo.png" alt="ggchord" width="320"></a>
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# ggchord：弦图、共线性图与质粒图

基于 ggplot2 的 R 包。用独立图层组合序列骨架、区间连接、基因、通用注释、标签和定量轨道。

正式版为 **ggchord 0.16.6**；本文示例由本地 **v0.27.0** 开发源码生成。正式资产见 GitHub Releases。

## 安装

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages("ggchord")
```

如已检出本版源码，在仓库根目录使用 `pkgload::load_all(".")`。正式版可从 CRAN 安装。 示例使用 dplyr（≥ 1.1.0）、purrr、tidyr 和原生 `|>`，这些是示例可选依赖。

## 三种模式

以下三张综合示例展示本版的组合能力。原始输入保持可追溯，完整准备代码和解释见对应教程；初次使用请从最小首图开始。

<details>
<summary>先运行：三张图共用的样式</summary>

```r
# 以下为可选外观设置。
seq_palette <- setNames(
  c("#E64B35", "#233E4A", "#269B82", "#D59B17"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1")
)
gene_palette <- c("+" = "#CE1536", "-" = "#346B9D", "." = "#8291A6")
```

```r
doc_background <- "#FFFCF5"
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
identity_guide_theme <- theme(
  legend.key.height = unit(100, "mm"), legend.key.width = unit(5, "mm"),
  legend.title = element_text(size = 13, face = "bold"),
  legend.text = element_text(size = 11), legend.margin = margin(r = 7, unit = "mm")
)
identity_guide <- ggplot2::guide_coloursteps(
  available_aes = "ribbon_fill", even.steps = FALSE, show.limits = TRUE,
  position = "left", direction = "vertical", theme = identity_guide_theme
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

按每条序列的真实正负链分别抽样，无 CDS 长度门槛。两条 MT 序列的原始注释只有正链 CDS。正链置于内侧，连接带在重叠端点局部避让。

![弦图：四基因组、混合曲率与相似度](man/figures/showcasechord.png)

<details>
<summary>数据准备与绘图代码 · 展开运行</summary>

```r
library(dplyr)
library(purrr)
library(tidyr)
# 使用仓库原始记录，先解压教学输入包。
show_sequences <- read_fasta_lengths(files = "examples/fasta/*.fna") |>
  select(accver, length)
show_annotations <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep")
show_genes <- show_annotations |>
  filter(type == "CDS") |>
  select(accver, start, end, strand, anno)
show_alignments <- read_blast(files = "examples/blastn-standard/*.tsv")
# 每对序列最多保留三条较长的实际比对，使用原始坐标。
show_links <- show_alignments |>
  filter(length >= 900, pident >= 75) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 3, with_ties = FALSE) |>
  ungroup()
# 每条序列沿注释顺序选取最多四条原始 CDS。
# 特征保留原始位置，完整源表单独保存。
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
# 补充尚未入选且间隔较大的负链 CDS。
show_reverse <- show_genes |>
  filter(strand == "-", end - start + 1 >= 500) |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  dplyr::slice(c(1, n())) |>
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
# 沿源坐标对实际出现的正负链分别抽样。
# 保留短 CDS 及其原始链和坐标。
# 每个正链分组优先选择与已展示连接端点重叠的 CDS。
chord_link_ends <- bind_rows(
  transmute(show_links, accver = qaccver, lo = pmin(qstart, qend), hi = pmax(qstart, qend)),
  transmute(show_links, accver = saccver, lo = pmin(sstart, send), hi = pmax(sstart, send))
)
chord_endpoint_genes <- show_genes |>
  filter(strand == "+") |>
  semi_join(chord_link_ends, by = join_by(accver, overlaps(start, end, lo, hi)))
chord_genes <- show_genes |>
  filter(strand %in% c("+", "-")) |>
  left_join(mutate(chord_endpoint_genes, at_endpoint = TRUE),
    by = c("accver", "start", "end", "strand", "anno")) |>
  group_by(accver, strand) |>
  arrange(start, .by_group = TRUE) |>
  mutate(display_bin = ntile(row_number(), if (first(strand) == "+") 4L else 3L)) |>
  group_by(accver, strand, display_bin) |>
  arrange(desc(coalesce(at_endpoint, FALSE)), is.na(anno) | anno == "hypothetical protein",
    start, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(accver, start, end, strand, anno)
chord_labels <- chord_genes
chord_position <- position_strand(c("+" = -0.030, "-" = 0.030))
```

```r
show_curvature <- c(
  "MT108731.1" = 0, "MT118296.1" = 1,
  "OQ646790.1" = -0.65, "OR222515.1" = 0.75
)
show_radius <- c(
  "MT108731.1" = 1.08, "MT118296.1" = 0.82,
  "OQ646790.1" = 1.18, "OR222515.1" = 1.02
)
p_show_chord <- ggchord(seq_data = show_sequences, ribbon_data = show_links) +
  geom_link_ribbon(link_avoid = "smooth", alpha = 0.46, linewidth = 0.14) +
  geom_seq(
    seq_orientation = show_optimization$seq_orientation,
    seq_curvature = show_curvature, seq_radius = show_radius,
    seq_gap = 0.06, linewidth = 1.1
  ) +
  geom_gene(data = chord_genes, position = chord_position, width = 0.025) +
  geom_gene_label_repel(
    data = chord_labels, label_layout = "radial", label_wrap = 24,
    size = 2.8, label_max_lines = 3, position = chord_position
  ) +
  geom_seq_label(size = 4, seq_label_radius = 0.77, seq_label_hjust = 0.5) +
  scale_seq_colour_manual(
    values = seq_palette, name = "Sequence",
    guide = sequence_guide
  ) +
  scale_gene_fill_manual(values = gene_palette, breaks = c("+", "-"), name = "Strand") +
  guides(gene_fill = strand_guide, ribbon_fill = identity_guide) +
  doc_axis() +
  coord_chord(rotation = 50) +
  doc_theme +
  labs(
    title = "Four genomes, one connected view",
    subtitle = "Mixed backbone geometry · CDS sampled by source strand"
  )
view_ggchord(p_show_chord, width = 12, height = 10, units = "in")
```

</details>

[查看完整教程 →](https://dangjem.github.io/ggchord/zh/chord.html)

### 共线性图：四序列与跨行连接

正负链分列、逐一标注的精选CDS、局部避障、虚线跨行ribbon与右侧控制点路由。

![共线性图：四序列与跨行连接](man/figures/showcasecollinear.png)

<details>
<summary>数据准备与绘图代码 · 展开运行</summary>

```r
library(dplyr)
library(purrr)
library(tidyr)
# 使用仓库原始记录，先解压教学输入包。
show_sequences <- read_fasta_lengths(files = "examples/fasta/*.fna") |>
  select(accver, length)
show_annotations <- read_gff3(files = "examples/gff3/*.gff3", unstranded = "keep")
show_genes <- show_annotations |>
  filter(type == "CDS") |>
  select(accver, start, end, strand, anno)
show_alignments <- read_blast(files = "examples/blastn-standard/*.tsv")
# 每对序列最多保留三条较长的实际比对，使用原始坐标。
show_links <- show_alignments |>
  filter(length >= 900, pident >= 75) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 3, with_ties = FALSE) |>
  ungroup()
# 每条序列沿注释顺序选取最多四条原始 CDS。
# 特征保留原始位置，完整源表单独保存。
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
# 补充尚未入选且间隔较大的负链 CDS。
show_reverse <- show_genes |>
  filter(strand == "-", end - start + 1 >= 500) |>
  group_by(accver) |>
  arrange(start, .by_group = TRUE) |>
  dplyr::slice(c(1, n())) |>
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
# 本张线性图省略两个较短且拥挤的候选特征。
# 两条原始记录保留在 show_features 和弦图中。
show_linear_features <- show_features |>
  filter(!anno %in% c("tail assembly protein", "gp4 protein"))
# 分别绘制相邻和跨行连接，保留原始端点。
show_pair_rows <- show_links |>
  mutate(
    query_row = match(qaccver, show_sequences$accver),
    subject_row = match(saccver, show_sequences$accver),
    span = abs(query_row - subject_row)
  )
show_adjacent <- show_pair_rows |> filter(span == 1)
# 每对非相邻序列选一条最长的实际比对，便于辨认长连接。
show_cross <- show_pair_rows |>
  filter(span > 1) |>
  group_by(qaccver, saccver) |>
  slice_max(length, n = 1, with_ties = FALSE) |>
  ungroup()
# 用与最终图相同的对齐和间隙设置检查骨架布局。
show_linear_coord <- coord_collinear(align = "center", orientation = "auto", ratio = 0.7)
show_routing_base <- ggchord(show_sequences, show_links) +
  geom_seq(seq_gap = 0.42) + show_linear_coord
show_reference <- get_chord_layout(show_routing_base)$sequence_reference
# 共享控制点放在中间骨架稍外侧。
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
  guides(gene_fill = strand_guide, ribbon_fill = identity_guide) +
  doc_axis() +
  show_linear_coord +
  doc_theme +
  labs(
    title = "Sequence architecture, across four genomes",
    subtitle = "Selected labelled CDS · strand-separated arrows · routed cross-row comparisons"
  )
view_ggchord(p_show_collinear, width = 12, height = 8, units = "in")
```

</details>

[查看完整教程 →](https://dangjem.github.io/ggchord/zh/collinear.html)

### 质粒图：pETDuet-1密集注释

双链骨架、20个特征、5个参考引物和58个合并切点；位置和名称来自真实输入。

![质粒图：pETDuet-1密集注释](man/figures/showcaseplasmid.png)

<details>
<summary>数据准备与绘图代码 · 展开运行</summary>

```r
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
p_show_plasmid <- ggchord(show_plasmid, feature_data = show_plasmid_features, validate = "none") +
  geom_seq(seq_style = "double", colour = "#252525") +
  geom_feature_plasmid(position = show_plasmid_position) +
  geom_feature_label_repel(
    position = show_plasmid_position,
    external = TRUE, max_overlaps = 0,
    label_fit = "shrink", label_min_size = 2.5
  ) +
  geom_primer(data = show_plasmid_primers, colour = "#A020F0") +
  geom_primer_label_repel(
    data = show_plasmid_primers, colour = "#A020F0", max_overlaps = 0
  ) +
  geom_restriction_site(data = show_plasmid_sites, colour = "#252525") +
  geom_seq_center_label() +
  coord_circular(rotation = 0) +
  theme_ggchord_plasmid() +
  theme(
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.spacing.x = unit(9, "mm"), legend.spacing.y = unit(7, "mm"),
    plot.background = element_rect(fill = doc_background, colour = NA),
    panel.background = element_blank()
  )
view_ggchord(p_show_plasmid, width = 16.5, height = 8, units = "in")
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

## 未来方向

后续开发继续围绕清晰的静态科研图和一致的绘图语法展开。当前函数、支持输入及限制见教程与参考文档；以下列出进一步改进的方向。

- 改善实际成图尺寸下的密集注释、引导线和图例布局，空间不足时提供明确诊断。
- 完善整体与局部细节的组合，保留原始坐标、单位及每条显示关系的含义。
- 扩展已有科研结果的适配，优先使用成熟的可选工具，保留源记录与缺失信息。
- 完善 API 一致性、兼容性和可重现的教学示例，保持普通数据框输入简洁。

## 设计参考与致谢

ggchord 的图形语法、基因组表达和预览流程参考了以下开源项目与公共工具：

- [ggplot2](https://ggplot2.tidyverse.org/)：以图层组合图形，区分数据映射、标度、几何、坐标与主题。
- [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/)：借鉴 query/subject 表格字段命名，便于追踪比对端点与序列标识。
- [Circos](https://circos.ca/) 和 [circlize](https://github.com/jokergoo/circlize)：参考环形信息层级、连接带表达与多轨道的阅读秩序。
- [clinker](https://github.com/gamcil/clinker)、[gggenes](https://github.com/wilkox/gggenes) 和 [genoPlotR](https://genoplotr.r-forge.r-project.org/)：参考基因箭头、基因簇与线性共线性比较的可读表达。
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) 和 [gggenomes](https://thackl.github.io/gggenomes/)：参考注释布局、多表输入、特征堆叠及局部比较工作流。
- [ggview](https://github.com/idmn/ggview)：参考按实际导出尺寸预览的流程；ggchord 用 `view_ggchord()` 预览，用标准 `ggsave()` 保存。

这些项目为设计与比较提供了重要参考。具体可用功能、参数及运行依赖以 ggchord 的函数文档和 `DESCRIPTION` 为准。

## 继续阅读

- [开始使用](https://dangjem.github.io/ggchord/zh/installation.html)
- [章节索引](https://dangjem.github.io/ggchord/zh/about.html)
- [函数文档](https://dangjem.github.io/ggchord/zh/reference/index.html)
- [示例图集](https://dangjem.github.io/ggchord/zh/gallery.html)
- [其他信息](https://dangjem.github.io/ggchord/zh/about.html)

[Issues](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
