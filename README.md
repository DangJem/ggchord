[English](README.md) · [现代汉语](README-Hans.md)

<img src="man/figures/logo.png" alt="ggchord" width="180" align="right">

# ggchord: chord, circular and collinear sequence plots

An R package built on ggplot2. Combine independent layers for sequence backbones, interval links, genes, general features, labels and quantitative tracks.

This README describes the local **0.16.5 candidate API**. The latest formal release remains 0.15.0. The new website source has not been deployed.

## Install

R ≥ 4.1.0 · ggplot2 ≥ 4.0.0

```r
install.packages("remotes")
remotes::install_github("DangJem/ggchord")
```

For a checkout of this version, use `pkgload::load_all(".")` from its root. A CRAN release is not claimed.

## Three views

These are synthetic teaching inputs, not biological observations.

```r
library(ggchord)
seqs <- data.frame(accver = c("A", "B", "C"), length = c(1000, 800, 1200))
links <- data.frame(qaccver = c("A", "B"), saccver = c("B", "C"),
  qstart = c(100, 200), qend = c(350, 500),
  sstart = c(150, 400), send = c(400, 700), pident = c(96, 89))
genes <- data.frame(accver = c("A", "B", "C"),
  start = c(150, 200, 400), end = c(300, 350, 650),
  strand = c("+", "+", "-"), anno = c("target", "target", "target"))
```

```r
p_chord <- ggchord(seqs, links, genes) +
  geom_seq() + geom_link_ribbon() + geom_gene() + geom_seq_label() +
  coord_chord()
```

```r
features <- data.frame(accver = "A", start = c(80, 380, 700),
  end = c(280, 610, 850), strand = c("+", "-", "."),
  type = c("CDS", "CDS", "region"), label = c("alpha", "beta", "origin"))
p_circular <- ggchord(seqs[1, ]) + geom_seq() +
  geom_feature(data = features) +
  geom_feature_label_repel(data = features) + coord_circular() +
  theme_ggchord(axis = element_blank())
```

```r
p_collinear <- ggchord(seqs, links, genes) +
  geom_link_ribbon(alpha = .3) + geom_seq() +
  geom_gene(position = position_strand(.1)) + geom_seq_label() +
  coord_collinear(anchor = genes[c("accver", "start", "end")],
    align = "center", orientation = "auto")
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

- [Full tutorial](vignettes/ggchord_vignette.Rmd): install → data → first plot → adjustments → advanced topics → FAQ.
- [Migration guide](MIGRATION.md) and [NEWS](NEWS.md).
- [Build and preview the local website](website/README.md), including reference categories, terminology, gallery and deployment design.

[Issues](https://github.com/DangJem/ggchord/issues) · [MIT license](LICENSE)
