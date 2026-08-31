# ggchord themes

Theme helpers remove Cartesian axes and provide consistent typography
and legend spacing for chord diagrams. \`theme_ggchord()\` preserves the
package's established default appearance; the other variants provide a
quieter canvas, a dark canvas, or a print-oriented white canvas.

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
