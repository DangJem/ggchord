# Get the chord layout from the package environment

Returns the most recently computed chord layout (after the plot was
built, e.g. via [`print()`](https://rdrr.io/r/base/print.html) or
`ggplot_build()`). This is useful for building custom layers or
annotations on top of the chord geometry.

## Usage

``` r
get_chord_layout()
```

## Value

A chord layout list containing the computed geometry (sequence arcs,
ribbon polygons, gene arrows, axis elements, extremes, colors, etc.)

## Examples

``` r
library(ggchord)
data(seq_data_example)
data(ribbon_data_example)
p <- ggchord(seq_data_example, ribbon_data_example) + geom_seq() + geom_ribbon()
invisible(ggplot2::ggplot_build(p))
names(get_chord_layout()$seq_arcs)
#> NULL
```
