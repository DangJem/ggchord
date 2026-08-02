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
