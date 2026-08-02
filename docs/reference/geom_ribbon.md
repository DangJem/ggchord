# Add an alignment ribbon layer

Draws colored ribbons corresponding to alignment results. Color scheme
and spacing parameters are specified here.

## Usage

``` r
geom_ribbon(
  mapping = NULL,
  data = NULL,
  ribbon_color_scheme = NULL,
  ribbon_colors = NULL,
  ribbon_alpha = NULL,
  ribbon_ctrl_point = NULL,
  ribbon_gap = NULL,
  alpha = NULL,
  ribbon_outline_color = "black",
  ribbon_outline_width = 0.05,
  ribbon_outline_linetype = 1,
  show_legend = TRUE,
  legend_position = NULL,
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Default NULL (retrieved automatically from the layout)

- ribbon_color_scheme:

  Character. Color scheme "pident", "query" or "single", default
  "pident"

- ribbon_colors:

  Optional color vector. Ribbon color parameters

- ribbon_alpha:

  Numeric (0-1). Ribbon transparency, default 0.35

- ribbon_ctrl_point:

  Optional vector/list. Bezier control points, default c(0,0)

- ribbon_gap:

  Optional numeric/vector. Spacing between sequences and ribbons,
  default 0.15

- alpha:

  Ribbon transparency (overrides ribbon_alpha), defaults to the value
  used in the layout

- ribbon_outline_color:

  Character. Color of the ribbon outline (border), default "black"

- ribbon_outline_width:

  Numeric. Line width of the ribbon outline, default 0.05

- ribbon_outline_linetype:

  Numeric or character. Line type of the ribbon outline, default 1
  (solid); see `linetype` in ggplot2 for options

- show_legend:

  Whether to show the legend, default TRUE

- legend_position:

  Optional position of this layer's legend (the Identity( "inside". When
  NULL the legend follows `theme(legend.position = ...)` together with
  the other legends. Can also be set with
  `theme(legend.position.ribbon = ...)`.

- ...:

  Additional arguments passed to
  [`geom_polygon()`](https://ggplot2.tidyverse.org/reference/geom_polygon.html)

## Value

A list of ggplot2 layers
