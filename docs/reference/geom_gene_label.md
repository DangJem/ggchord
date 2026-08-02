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
  gene_label_repel = FALSE,
  gene_label_wrap = NULL,
  gene_label_max_overlaps = Inf,
  gene_label_seed = 123,
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

- gene_label_repel:

  Logical, default FALSE. When TRUE, overlapping gene labels are
  automatically pushed apart (collision detection + automatic avoidance)
  so they do not cover each other.

- gene_label_wrap:

  Numeric or NULL, default NULL. When set, long gene annotations are
  wrapped at this many characters (e.g. 15), which makes the labels
  narrower and less prone to overlap.

- gene_label_max_overlaps:

  Numeric, default Inf. With `gene_label_repel = TRUE`, labels that
  still overlap more than this many other labels after de-overlapping
  are hidden (ggrepel-style). Use a finite value to declutter crowded
  plots.

- gene_label_seed:

  Numeric, default 123. Seed used by the de-overlap algorithm for
  reproducible results.

- show_legend:

  Whether to show the legend, default FALSE

- ...:

  Additional arguments passed to
  [`geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)

## Value

A list of ggplot2 layers

## Details

Long or crowded labels can be handled with `gene_label_repel` (automatic
de-overlapping), `gene_label_wrap` (wrapping long annotations) and
`gene_label_max_overlaps` (hiding labels that still overlap too many
others).
