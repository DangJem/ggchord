# Add an automatically arranged gene label layer

Like
[`geom_gene_label()`](https://dangjem.github.io/ggchord/reference/geom_gene_label.md),
but labels are placed by one of three deterministic collision-avoiding
layouts. The default `"aligned"` layout uses orderly cardinal rails,
`"radial"` uses compact local offset tracks, and `"arc"` keeps text
close to and rotated with the sequence curve.

## Usage

``` r
geom_gene_label_repel(
  mapping = NULL,
  data = NULL,
  gene_label_layout = "aligned",
  gene_label_size = NULL,
  gene_label_wrap = NULL,
  gene_label_side = "outside",
  max_overlaps = Inf,
  gene_label_segment_linetype = "auto",
  show_legend = FALSE,
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Default NULL (retrieved automatically from the layout)

- gene_label_layout:

  Character, default `"aligned"`. Label layout: `"aligned"` uses
  horizontal labels on orderly top, bottom, left and right rails;
  `"radial"` uses horizontal labels on the nearest collision-free local
  offset track; `"arc"` rotates labels along the sequence tangent and
  keeps them close to their genes.

- gene_label_size:

  Numeric. Label font size, default 2.5

- gene_label_wrap:

  Numeric or NULL, default NULL. When set, long gene annotations are
  wrapped at this many characters (e.g. 15).

- gene_label_side:

  Character, default "outside". Which side of the arc the labels sit on.
  `"auto"` keeps the strand-based placement (same as before);
  `"outside"` moves labels that would be inside the chord (where they
  can overlap the ribbons) to the outside of their arc; `"inside"` does
  the opposite. Labels moved to the other side are connected with a
  dashed leader line (see `gene_label_segment_linetype`).

- max_overlaps:

  Numeric, default Inf. Hide labels that still overlap more than this
  many other labels after repulsion (ggrepel-style decluttering). Use a
  finite value to clean up crowded plots.

- gene_label_segment_linetype:

  Character or numeric, default "auto". Leader-line linetype. `"auto"`
  draws solid lines, except for labels that were moved to the other side
  of their arc, which are drawn dashed. Any other valid ggplot2 linetype
  (e.g. `"solid"`, `"dashed"`, `"dotted"`, or a numeric dash pattern) is
  used for all leader lines.

- show_legend:

  Whether to show the legend, default FALSE

- ...:

  Additional arguments passed to `geom_text()`

## Value

A list of ggplot2 layers (a leader-line layer and a text layer).

## Details

The local outside/inside label concepts are informed by SnapGene and
Geneious, but ggchord uses generic mode names and an independent
geometry implementation. See [SnapGene feature
labels](https://support.snapgene.com/hc/en-us/articles/10383722725524-Display-Feature-Labels-Below-or-Inside-a-Map)
and [Geneious label
options](https://manual.geneious.com/en/latest/Sequences.html).

Low-level force, padding, orientation and segment arguments used by
earlier releases have been removed and now produce an error. Use
`gene_label_layout` for automatic placement, or
[`geom_gene_label()`](https://dangjem.github.io/ggchord/reference/geom_gene_label.md)
for manual rotation and offsets.

## Examples

``` r
library(ggchord)
data(seq_data_example)
data(gene_data_example)
p <- ggchord(seq_data_example, gene_data = gene_data_example) +
  geom_seq() + geom_gene() + geom_gene_label_repel()
p
```
