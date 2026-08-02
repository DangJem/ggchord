# ggchord: layered multi-sequence alignment chord diagrams for ggplot2

ggchord visualizes multi-sequence alignment results using ggplot2's
layered grammar. The `ggchord()` constructor handles data validation and
global settings; the `geom_*` layers are stacked as needed, each
responsible for its own layout parameters and visual rendering. The
layout is computed lazily when the plot is built (e.g. via
[`print()`](https://rdrr.io/r/base/print.html),
[`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html), or
[`ggplot_build()`](https://ggplot2.tidyverse.org/reference/ggplot_build.html)).

## Usage

``` r
ggchord(
  seq_data,
  ribbon_data = NULL,
  gene_data = NULL,
  title = NULL,
  rotation = 45,
  panel_margin = 0,
  show_legend = TRUE,
  debug = FALSE
)
```

## Arguments

- seq_data:

  data.frame/tibble, required. Basic sequence information

- ribbon_data:

  data.frame/tibble, optional. Alignment results

- gene_data:

  data.frame/tibble, optional. Gene annotation data

- title:

  Character. Main title of the plot, default NULL

- rotation:

  Numeric. Global rotation angle (degrees), default 45

- panel_margin:

  Optional numeric/list. Panel margin, default 0

- show_legend:

  Logical. Whether to show legends, default TRUE

- debug:

  Logical. Whether to output debug information, default FALSE

## Value

A ggchord object (inherits from ggplot) to which geom\_\* layers can be
added with +

## Examples

``` r
library(ggchord)
data(seq_data_example)
data(ribbon_data_example)
data(gene_data_example)

p <- ggchord(
  seq_data = seq_data_example,
  ribbon_data = ribbon_data_example,
  gene_data = gene_data_example
) +
  geom_seq() +
  geom_ribbon() +
  geom_gene() +
  geom_axis()
print(p)

```
