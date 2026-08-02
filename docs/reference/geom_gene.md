# Add a gene arrow layer

Draws gene annotation arrows on the chord diagram. Gene layout
parameters (offset, width, color scheme, etc.) are specified here. The
gene fill scale is kept independent from the ribbon's fill scale via a
separate internal aesthetic used by the ribbon layer.

## Usage

``` r
geom_gene(
  mapping = NULL,
  data = NULL,
  gene_offset = NULL,
  gene_width = NULL,
  gene_color_scheme = NULL,
  gene_colors = NULL,
  gene_order = NULL,
  gene_label_show = NULL,
  gene_label_size = NULL,
  gene_label_rotation = NULL,
  gene_label_radial_offset = NULL,
  gene_label_circum_offset = NULL,
  gene_label_circum_limit = NULL,
  show_legend = TRUE,
  show_label = NULL,
  label_size = NULL,
  legend_position = NULL,
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Default NULL (retrieved automatically from the layout)

- gene_offset:

  Optional numeric/vector/list. Radial offset of gene arrows, default
  0.1

- gene_width:

  Optional numeric/vector/list. Width of gene arrows, default 0.05

- gene_color_scheme:

  Character. "strand" or "manual", default "strand"

- gene_colors:

  Optional color vector. Fill color of gene arrows

- gene_order:

  Optional character vector. Display order of genes in the legend

- gene_label_show:

  Logical. Whether to show gene labels, default FALSE

- gene_label_size:

  Numeric. Label font size, default 2.5

- gene_label_rotation:

  Optional numeric/vector/list. Label rotation angle, default 0

- gene_label_radial_offset:

  Optional numeric/vector/list. Radial offset of labels, default 0

- gene_label_circum_offset:

  Optional numeric/vector/list. Circumferential offset of labels,
  default 0

- gene_label_circum_limit:

  Optional logical/vector/list. Whether to limit circumferential offset,
  default TRUE

- show_legend:

  Whether to show the legend, default TRUE

- show_label:

  Whether to show gene labels (overrides gene_label_show), default NULL

- label_size:

  Label font size (overrides gene_label_size), default NULL

- legend_position:

  Optional position of this layer's legend (the Strand or Gene
  Annotation legend): one of "left", "right", "top", "bottom" or
  "inside". When NULL the legend follows `theme(legend.position = ...)`
  together with the other legends. Can also be set with
  `theme(legend.position.gene = ...)`.

- ...:

  Additional arguments passed to
  [`geom_polygon()`](https://ggplot2.tidyverse.org/reference/geom_polygon.html)

## Value

A list of ggplot2 layers
