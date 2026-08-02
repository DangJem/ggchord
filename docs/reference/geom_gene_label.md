# Add a gene label layer

Draws the gene annotation labels on a chord diagram. This layer is
independent from
[`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md):
add it after
[`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md)
to annotate the gene arrows with their texts.

## Usage

``` r
geom_gene_label(
  mapping = NULL,
  data = NULL,
  gene_label_size = NULL,
  gene_label_rotation = NULL,
  gene_label_radial_offset = NULL,
  gene_label_circum_offset = NULL,
  gene_label_circum_limit = NULL,
  gene_label_wrap = NULL,
  show_legend = FALSE,
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Default NULL (retrieved automatically from the layout)

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

- gene_label_wrap:

  Numeric or NULL, default NULL. When set, long gene annotations are
  wrapped at this many characters (e.g. 15), which makes the labels
  narrower and less prone to overlap.

- show_legend:

  Whether to show the legend, default FALSE

- ...:

  Additional arguments passed to
  [`geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)

## Value

A list of ggplot2 layers. To let the labels avoid each other and the
genes (with leader lines), use
[`geom_gene_label_repel()`](https://dangjem.github.io/ggchord/reference/geom_gene_label_repel.md)
instead.

## Details

Long annotations can be wrapped with `gene_label_wrap`. For automatic
de-overlapping (with leader lines), use
[`geom_gene_label_repel()`](https://dangjem.github.io/ggchord/reference/geom_gene_label_repel.md)
instead.
