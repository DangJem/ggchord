# Add a sequence label layer

Places sequence labels at the midpoint of each sequence arc, radially
offset from the arc. Labels can be styled via their radial offset,
rotation and font size.

## Usage

``` r
geom_seq_label(
  mapping = NULL,
  data = NULL,
  seq_label_radius = NULL,
  seq_label_rotation = NULL,
  seq_label_size = NULL,
  seq_labels = NULL,
  show_legend = FALSE,
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Default NULL (retrieved automatically from the layout)

- seq_label_radius:

  Optional numeric/vector. Radial position of the labels as a multiplier
  of the sequence arc radius (e.g. 1.15 = 15 the arc), default NULL
  (1.15)

- seq_label_rotation:

  Optional numeric/vector. Additional label rotation (degrees), default
  NULL (0)

- seq_label_size:

  Optional numeric/vector. Label font size, default NULL (3)

- seq_labels:

  Optional character vector. Override the label texts (defaults to the
  sequence labels from
  [`geom_seq()`](https://dangjem.github.io/ggchord/reference/geom_seq.md)
  or the sequence IDs)

- show_legend:

  Whether to show the legend, default FALSE

- ...:

  Additional arguments passed to
  [`geom_text()`](https://rdrr.io/pkg/ggplot2/man/geom_text.html)

## Value

A list of ggplot2 layers
