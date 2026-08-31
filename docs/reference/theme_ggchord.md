# ggchord themes

Theme helpers remove Cartesian axes and provide consistent typography
and legend spacing for chord diagrams. \`theme_ggchord()\` uses a
restrained white publication canvas; the other variants provide a still
quieter canvas, a dark canvas, or tighter print typography.

## Usage

``` r
theme_ggchord(base_size = 11, base_family = "")

theme_ggchord_minimal(base_size = 11, base_family = "")

theme_ggchord_dark(base_size = 11, base_family = "")

theme_ggchord_publication(base_size = 9, base_family = "")
```

## Arguments

- base_size:

  Base font size in points.

- base_family:

  Base font family.

## Value

A ggplot2 theme object.
