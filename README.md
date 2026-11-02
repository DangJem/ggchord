<p align="center">
  <a href="https://dangjem.github.io/ggchord/en/index.html"><img src="man/figures/logo.png" alt="ggchord" width="320"></a>
</p>

[English](README.md) · [现代汉语](README-Hans.md)

# ggchord: chord, circular and collinear sequence plots

An R package built on ggplot2. Combine independent layers for sequence backbones, interval links, genes, general features, labels and quantitative tracks.

This README documents **ggchord 0.16.5**. Read the migration guide when updating older code; formal release assets are listed on GitHub Releases.

## Install

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages(c("remotes", "dplyr"))
remotes::install_github("DangJem/ggchord@v0.16.5")
```

For a checkout of this version, use `pkgload::load_all(".")` from its root. A CRAN release is not claimed. These examples prepare tables with dplyr and the native `|>` pipe; dplyr is an example dependency.

## Three views

Select only the required fields from the bundled source-backed data before plotting; keep original metadata in the source objects.

<details>
<summary>Run first: shared data and appearance for chord and collinear views</summary>

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

### Chord: compare overall sequence relationships

Show relationships between aligned intervals.

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

![Chord: compare overall sequence relationships](man/figures/chord.png)

[Read the full guide →](https://dangjem.github.io/ggchord/en/real-data.html)

### Circular: read annotations on one sequence

Prepare the pUC19 annotations, then combine backbone and feature layers.

<details>
<summary>Run first: circular-map data and colours</summary>

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

![Circular: read annotations on one sequence](man/figures/circular.png)

[Read the full guide →](https://dangjem.github.io/ggchord/en/plasmids.html)

### Collinear: align and compare positions

Reuse the first section’s data and appearance, aligning at explicitly selected gene anchors.

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

![Collinear: align and compare positions](man/figures/collinear.png)

[Read the full guide →](https://dangjem.github.io/ggchord/en/modes.html)

`aes()` maps data; scales define mapping and guides; geoms define geometry and fixed styles; Positions place local geometry; coords select the global view; themes style non-data elements.

## Preview and save

RStudio Plots can rearrange labels when its panel changes. Preview and save at the same explicit dimensions:

```r
view_ggchord(p_collinear, width = 8, height = 5, units = "in")
ggsave("comparison.png", p_collinear, width = 8, height = 5, units = "in", dpi = 144)
```

## Continue

- [Get started](https://dangjem.github.io/ggchord/en/installation.html)
- [Task guides](https://dangjem.github.io/ggchord/en/articles.html)
- [Function reference](https://dangjem.github.io/ggchord/en/reference/index.html)
- [Gallery](https://dangjem.github.io/ggchord/en/gallery.html)
- [Migration guide](https://dangjem.github.io/ggchord/en/migration.html)
- [Changelog](https://dangjem.github.io/ggchord/en/news.html)

[Issues](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
