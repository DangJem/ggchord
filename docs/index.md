# ggchord: Layered Multi-Sequence Chord Diagrams for ggplot2

## Overview

`ggchord` is an R package built on `ggplot2` for drawing **chord
diagrams of multi-sequence data** using the layered grammar of graphics.
You provide three tidy data frames (`seq_data`, `ribbon_data`,
`gene_data`) and stack `geom_*` layers like any other ggplot2 plot.

## Key Features

- True `ggplot2` layered style
- Multi-sequence, multi-ribbon and multi-gene support
- Per-sequence and per-strand parameters
- Themes, scales, `ggsave()`, `ggplot_build()` and
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)
  integration

## Installation

`ggchord` requires R (≥ 4.1.0) and `ggplot2` (≥ 4.0.0).

``` r

install.packages("ggchord")            # CRAN
devtools::install_github("DangJem/ggchord")   # development version
```

## Get started

| Need | Resource |
|----|----|
| Full tutorial | [ggchord tutorial](https://dangjem.github.io/ggchord/articles/ggchord_vignette.html) |
| 汉语完整指南 | [ggchord 汉语指南](https://dangjem.github.io/ggchord/articles/ggchord_guide_hans.html) |
| Function reference | [Function reference](https://dangjem.github.io/ggchord/reference/index.html) |
| Version history | [NEWS.md](https://dangjem.github.io/ggchord/NEWS.md) |

## Contributions & Feedback

Bug reports and feature requests are welcome via [GitHub
issues](https://github.com/DangJem/ggchord/issues). Pull requests are
also appreciated.

## License

This project is licensed under the MIT License — see the
[LICENSE](https://dangjem.github.io/ggchord/LICENSE) file for details.
