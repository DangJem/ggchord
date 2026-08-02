# Add a repelled gene label layer (ggrepel-style)

Like
[`geom_gene_label()`](https://dangjem.github.io/ggchord/reference/geom_gene_label.md),
but the labels are placed with a force-based simulation that pushes them
away from the genes and from each other (similar to
[`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)).
Labels that move far enough from their anchor are connected to it with a
leader line, and labels that still overlap too many others can be
hidden.

## Usage

``` r
geom_gene_label_repel(
  mapping = NULL,
  data = NULL,
  gene_label_size = NULL,
  gene_label_rotation = NULL,
  gene_label_radial_offset = NULL,
  gene_label_circum_offset = NULL,
  gene_label_circum_limit = NULL,
  gene_label_wrap = NULL,
  max_overlaps = Inf,
  box_padding = 0.25,
  point_padding = 0.1,
  min_segment_length = 0.5,
  force = 1,
  seed = 123,
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
  wrapped at this many characters (e.g. 15).

- max_overlaps:

  Numeric, default Inf. Hide labels that still overlap more than this
  many other labels after repulsion (ggrepel-style decluttering). Use a
  finite value to clean up crowded plots.

- box_padding:

  Numeric, default 0.25. Extra padding around each label box (data
  units).

- point_padding:

  Numeric, default 0.1. Extra padding around the anchor points (data
  units).

- min_segment_length:

  Numeric, default 0.5. Labels that moved less than this distance (data
  units) from their anchor do not draw a leader line.

- force:

  Numeric, default 1. Strength of the repulsive forces.

- seed:

  Numeric, default 123. Random seed for reproducibility.

- show_legend:

  Whether to show the legend, default FALSE

- ...:

  Additional arguments passed to
  [`geom_text()`](https://ggplot2.tidyverse.org/reference/geom_text.html)

## Value

A list of ggplot2 layers (a leader-line layer and a text layer).
