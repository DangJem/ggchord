# Chord diagram coordinate system

A lightweight Coord. It creates placeholder coordinates in
[`ggchord()`](https://dangjem.github.io/ggchord/reference/ggchord.md),
which are replaced at print time with the actual extents computed from
the layout.

## Usage

``` r
coord_chord(layout = NULL)
```

## Arguments

- layout:

  Chord layout object (passed internally by ggchord(), may be NULL)

## Value

A Coord object for ggplot2 `+` composition

## Examples

``` r
library(ggchord)
data(seq_data_example)
p <- ggchord(seq_data_example) + coord_chord() + geom_seq()
#> Coordinate system already present.
#> ℹ Adding new coordinate system, which will replace the existing one.
p
```
