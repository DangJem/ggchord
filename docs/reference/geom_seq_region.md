# Highlight regions along sequence arcs

Draws rectangular bands on sequence arcs for one or more coordinate
intervals. This is useful for marking loci, repeats, CRISPR arrays, or
other user-defined regions without turning them into gene arrows.

## Usage

``` r
geom_seq_region(
  mapping = NULL,
  data = NULL,
  regions = NULL,
  region_fill = "#F59E0B",
  region_color = "#B45309",
  region_alpha = 0.25,
  region_width = 0.08,
  region_offset = 0,
  region_side = c("inside", "outside", "auto"),
  show_legend = FALSE,
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Optional data.frame; an alias for `regions` when `regions` is NULL.

- regions:

  data.frame with at least `seq_id`, `start`, `end`; optional `label`,
  `category` and `color`.

- region_fill:

  Character. Default fill colour for regions, default `"#F59E0B"`.

- region_color:

  Character. Outline colour, default `"#B45309"`.

- region_alpha:

  Numeric (0-1). Region alpha, default 0.25.

- region_width:

  Numeric. Band width in chord radius units, default 0.08.

- region_offset:

  Numeric. Radial offset from the sequence arc, default 0.

- region_side:

  Character. `"inside"`, `"outside"` or `"auto"`; `"auto"` places
  regions inside the chord when possible.

- show_legend:

  Logical. Whether to show a legend for category colours, default FALSE.

- ...:

  Additional arguments passed to `geom_polygon()`

## Value

A list of ggplot2 layers

## Examples

``` r
library(ggchord)
data(seq_data_example)
regions <- data.frame(seq_id = "MT108731.1",
                      start = 1000, end = 4000,
                      color = "orange")
p <- ggchord(seq_data_example) + geom_seq() + geom_seq_region(regions)
p
```
