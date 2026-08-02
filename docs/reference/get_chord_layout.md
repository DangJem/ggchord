# Get the chord layout from the package environment

Returns the most recently computed chord layout (after the plot was
built, e.g. via [`print()`](https://rdrr.io/r/base/print.html) or
[`ggplot_build()`](https://rdrr.io/pkg/ggplot2/man/ggplot_build.html)).
This is useful for building custom layers or annotations on top of the
chord geometry.

## Usage

``` r
get_chord_layout()
```

## Value

A chord layout list containing the computed geometry (sequence arcs,
ribbon polygons, gene arrows, axis elements, extremes, colors, etc.)
