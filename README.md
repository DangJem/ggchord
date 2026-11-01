[English](README.md) · [现代汉语](README-Hans.md)

<img src="man/figures/logo.png" alt="ggchord" width="200" align="right">

# ggchord: chord, circular and collinear sequence plots

An R package built on ggplot2. Combine independent layers for sequence backbones, interval links, genes, general features, labels and quantitative tracks.

This README documents **ggchord 0.16.5**. Read the migration guide when updating older code; formal release assets are listed on GitHub Releases.

## Install

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord@v0.16.5")
```

For a checkout of this version, use `pkgload::load_all(".")` from its root. A CRAN release is not claimed.

![Collinear case with real sequences and all annotations](man/figures/reallinear.png)

Original repository FASTA, GFF3 and BLAST inputs; see the [real-data workflow](vignettes/ggchord_vignette.Rmd#real-data) for complete import, filtering code and interpretation.


## Three views

Select only the required fields from the bundled source-backed data before plotting; keep original metadata in the source objects.

```r
library(ggchord)
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)
seqs <- seq_data_example[c("accver", "length")]
links <- ribbon_data_example[c("qaccver", "saccver", "qstart", "qend", "sstart", "send")]
# The fixture preserves biological strands separately from its old display strands.
genes <- transform(gene_data_example, strand = source_strand)[c("accver", "start", "end", "strand")]
```

```r
# Optional appearance choices; none of these columns or scales is required input.
library(ggchord)
seq_palette <- setNames(c("#3867D6", "#E58A28", "#219B83", "#9B59B6"),
  c("MT108731.1", "MT118296.1", "OQ646790.1", "OR222515.1"))
gene_palette <- c("+" = "#3B6FD8", "-" = "#E1812C", "." = "#8291A6")
```

```r
library(ggchord)
feature_palette <- c(CDS = "#3973CB", promoter = "#8B5EC4", rep_origin = "#E8A42D",
  primer_bind = "#168F86", misc_feature = "#A0ADBD", protein_bind = "#CC6485")
```

```r
library(ggchord)
doc_theme <- theme_ggchord_publication(axis.text = element_text(size = 11),
  axis.text.check.overlap = TRUE) + theme(legend.position = "bottom")
doc_axis <- function(coordinates = "relative") scale_seq_position_continuous(
  coordinates = coordinates, n.breaks = 3, minor_breaks = NULL,
  labels = function(x) format(x, scientific = FALSE, big.mark = ",", trim = TRUE))
```

```r
p_chord <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .7) +
  scale_ribbon_fill_manual(values = seq_palette, guide = "none") + geom_seq(show.legend = FALSE) + geom_gene(position = position_strand(.1), width = .04) +
  geom_seq_label(size = 4, seq_label_radius = .84, seq_label_hjust = .5) + coord_chord() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) + doc_theme
p_chord
```

```r
data(plasmid_example_pUC19)
plasmid_seq <- plasmid_example_pUC19[c("accver", "length")]
feature_source <- find_common_features(plasmid_example_pUC19)
features <- feature_source[c("accver", "start", "end", "strand", "type", "anno")]
```

```r
circle_features <- features[features$type %in% c("CDS", "rep_origin"), ]
p_circular <- ggchord(seq_data = plasmid_seq) + geom_seq() +
  geom_feature(aes(feature_fill = type), data = circle_features, position = position_plasmid()) +
  geom_feature_label_repel(data = circle_features, position = position_plasmid(), size = 4) +
  geom_seq_center_label() + scale_feature_fill_manual(values = feature_palette, breaks = unique(circle_features$type)) +
  coord_circular() + theme_ggchord_publication(axis.text = element_text(size = 11))
p_circular
```

```r
anchor_rows <- unlist(lapply(split(seq_len(nrow(genes)), genes$accver), function(i)
  i[which.max(genes$end[i] - genes$start[i])]))
anchors <- genes[anchor_rows, c("accver", "start", "end")]
p_collinear <- ggchord(seq_data = seqs, ribbon_data = links, gene_data = genes) +
  geom_link_ribbon(aes(ribbon_fill = qaccver), alpha = .7) +
  scale_ribbon_fill_manual(values = seq_palette, guide = "none") + geom_seq() +
  geom_gene(position = position_strand(.1)) + geom_seq_label(size = 4) +
  coord_collinear(anchor = anchors, align = "center", orientation = "auto") + doc_axis() +
  scale_seq_colour_manual(values = seq_palette) +
  scale_gene_fill_manual(values = gene_palette) + doc_theme
p_collinear
```

| Chord / 弦图 | Circular / 圆图 | Collinear / 共线性图 |
| --- | --- | --- |
| ![Chord](man/figures/chord.png) | ![Circular](man/figures/circular.png) | ![Collinear](man/figures/collinear.png) |

`aes()` maps data; scales define mapping and guides; geoms define geometry and fixed styles; Positions place local geometry; coords select the global view; themes style non-data elements.

## Preview and save

RStudio Plots can rearrange labels when its panel changes. Preview and save at the same explicit dimensions:

```r
view_ggchord(p_collinear, width = 8, height = 5, units = "in")
ggsave("comparison.png", p_collinear, width = 8, height = 5, units = "in", dpi = 144)
```

## Continue

- [User guide](vignettes/ggchord_vignette.Rmd): real data, plasmids, formats, links, labels, tracks, projection, homology, styling, output and FAQ.
- [Migration guide](MIGRATION.md) and [NEWS](NEWS.md).
- [Build and preview the local website](website/README.md), including reference categories, terminology, gallery and deployment design.

[Issues](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
