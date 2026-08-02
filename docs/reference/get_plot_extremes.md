# Calculate plot extremes

Extracts x/y coordinate extremes from all plot elements (sequence arcs,
ribbons, gene arrows, etc.) for adjusting the plot range

## Usage

``` r
get_plot_extremes(
  allRibbon = NULL,
  seqArcs = NULL,
  axisLines = NULL,
  axisTicks = NULL,
  gene_arrows = NULL,
  gene_polys = NULL,
  seq_labels = NULL,
  show_axis = FALSE
)
```

## Arguments

- allRibbon:

  data.frame, ribbon data (with x, y columns), default NULL

- seqArcs:

  List, sequence arc data (each element is a data frame with x, y,
  seq_id), default NULL

- axisLines:

  data.frame, axis line data (with x, y, seq_id columns), default NULL

- axisTicks:

  data.frame, tick mark data (with x0, y0, x1, y1, label_x, label_y
  columns), default NULL

- gene_arrows:

  data.frame, gene label data (with text_x, text_y columns), default
  NULL

- gene_polys:

  data.frame, gene arrow polygon data (with x, y columns), default NULL

- show_axis:

  Logical, whether to include extreme value calculation for axis-related
  elements, default FALSE

## Value

List containing x_min (minimum x), x_max (maximum x), y_min (minimum y),
y_max (maximum y)
