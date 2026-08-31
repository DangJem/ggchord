# Draw generic genomic features

A general layer for CDS, tRNA, rRNA, repeat, CRISPR, promoter or
user-defined features. It supports directional arrows, blocks, notched
chevrons and lollipops, all generated against the sequence's real local
curve. Feature fill and geometry use independent role-specific scales.

## Usage

``` r
geom_feature(
  mapping = NULL,
  data = NULL,
  type = "type",
  category = NULL,
  label = "label",
  feature_shape = "arrow",
  feature_colors = NULL,
  feature_width = NULL,
  feature_offset = NULL,
  feature_order = NULL,
  show_legend = TRUE,
  legend_position = "right",
  ...
)
```

## Arguments

- mapping:

  Optional aesthetic mapping. Role aesthetics such as `seq_id`, `start`,
  `end` and `strand` may rename input columns; ordinary visual mappings
  are evaluated after geometry is generated.

- data:

  data.frame with `seq_id`, `start`, `end` and `strand`; optional
  `type`, `category` and `label`.

- type:

  Column name used as the feature type, default `"type"`.

- category:

  Optional column name used for colour grouping; defaults to `type`.

- label:

  Optional column name used for annotation text; defaults to `label`
  when present, otherwise `type`.

- feature_shape:

  Fixed feature geometry used when `feature_shape` is not mapped in
  `aes()`: `"arrow"`, `"block"`, `"chevron"`, or `"lollipop"`. The
  default is `"arrow"`. Use `aes(feature_shape = type)` together with
  [`scale_feature_shape_manual()`](https://dangjem.github.io/ggchord/reference/scale_feature_shape_manual.md)
  to map categories to geometry.

- feature_colors:

  Optional named color vector by feature value; unnamed vectors are
  recycled positionally.

- feature_width:

  Optional numeric or named vector controlling feature width; passed to
  `geom_gene(gene_width = ...)`.

- feature_offset:

  Optional numeric or named vector controlling feature offset; passed to
  `geom_gene(gene_offset = ...)`.

- feature_order:

  Optional feature order for the legend.

- show_legend:

  Logical. Show the feature legend, default TRUE.

- legend_position:

  Position of the feature legend: `"left"`, `"right"`, `"top"`,
  `"bottom"` or `"inside"`.

- ...:

  Additional arguments passed to
  [`geom_gene()`](https://dangjem.github.io/ggchord/reference/geom_gene.md).

## Value

A list of ggplot2 layers

## Examples

``` r
library(ggchord)
data(seq_data_example)
features <- data.frame(seq_id = "MT108731.1",
                       start = 1000, end = 4000,
                       strand = "+", type = "CDS")
p <- ggchord(seq_data_example) + geom_seq() + geom_feature(features)
p
```
