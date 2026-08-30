🌐 Language: [简体中文](README-Hans.md) | [English](README.md)

# ggchord

`ggchord` is an R package for drawing layered multi-sequence chord diagrams
with `ggplot2`. It visualises sequence arcs, alignment ribbons, genomic
features, labels and position axes from tidy data frames.

## Installation

```r
install.packages("ggchord")
```

Install the development version from GitHub when needed:

```r
devtools::install_github("DangJem/ggchord")
```

## Quick start

```r
library(ggchord)

data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

ggchord(seq_data_example, ribbon_data_example, gene_data_example) +
  geom_seq() +
  geom_ribbon() +
  geom_gene() +
  geom_axis()
```

The required columns are:

| Data | Required columns |
| --- | --- |
| Sequences | `seq_id`, `length` |
| Alignments | `qaccver`, `saccver`, `length`, `pident`, `qstart`, `qend`, `sstart`, `send` |
| Features | `seq_id`, `start`, `end`, `strand`, `anno` |

## Documentation

- [English tutorial](https://dangjem.github.io/ggchord/articles/ggchord_vignette.html)
- [中文教程](https://dangjem.github.io/ggchord/articles/ggchord_guide_hans.html)
- [Function reference](https://dangjem.github.io/ggchord/reference/index.html)
- [Version history](NEWS.md)

Bug reports and feature requests are welcome through
[GitHub Issues](https://github.com/DangJem/ggchord/issues).

## License

MIT. See [LICENSE](LICENSE).
