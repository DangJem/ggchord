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
  ribbon_color_by = NULL,
  ribbon_color_limits = NULL,
  ribbon_color_breaks = NULL,
  ribbon_color_name = NULL,
  ribbon_alpha = NULL,
  ribbon_alpha_by = NULL,
  ribbon_alpha_range = c(0.15, 0.9),
  ribbon_ctrl_point = NULL,
  ribbon_gap = NULL,
  alpha = NULL,
  ribbon_outline_color = "black",
  ribbon_outline_width = 0.05,
  ribbon_outline_linetype = 1,
  ribbon_outline_by = NULL,
  ribbon_outline_colors = NULL,
  ribbon_linetype_by = NULL,
  ribbon_linetypes = NULL,
  ribbon_direction = c("none", "alpha", "outline", "linetype"),
  ribbon_direction_colors = c(same = "black", reverse = "grey50"),
  ribbon_direction_linetypes = c(same = "solid", reverse = "dashed"),
  ribbon_direction_alpha = c(same = 1, reverse = 0.45),
  show_legend = TRUE,
  legend_position = "left",
  legend_key_width = NULL,
  legend_key_height = NULL,
  ...
)
```

## Arguments

- mapping:

  Default NULL (uses pre-computed data)

- data:

  Default NULL (retrieved automatically from the layout)

- ribbon_color_scheme:

  Character. Color scheme `"pident"`, `"query"`, `"subject"` or
  `"single"`, default `"pident"`

- ribbon_colors:

  Optional color vector. Ribbon color parameters

- ribbon_color_by:

  Optional character column name. When set, ribbon fill is mapped to a
  continuous colourbar for that numeric column instead of `pident` (e.g.
  `"bitscore"`).

- ribbon_color_limits:

  Optional numeric length-2 limits for `ribbon_color_by`.

- ribbon_color_breaks:

  Optional numeric breaks for the `ribbon_color_by` colourbar.

- ribbon_color_name:

  Optional legend title for the `ribbon_color_by` colourbar (defaults to
  the column name).

- ribbon_alpha:

  Numeric (0-1). Ribbon transparency, default 0.35

- ribbon_alpha_by:

  Optional character column name. When set, alpha is scaled continuously
  from that numeric column.

- ribbon_alpha_range:

  Numeric length-2. Alpha range used by `ribbon_alpha_by`, default
  `c(0.15, 0.9)`.

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

- ribbon_outline_by:

  Optional discrete column name. When set, outline colour is mapped by
  that column and `ribbon_outline_colors` controls the palette.

- ribbon_outline_colors:

  Optional named color vector for `ribbon_outline_by`; unnamed vectors
  are recycled positionally.

- ribbon_linetype_by:

  Optional discrete column name. When set, outline linetype is mapped by
  that column and `ribbon_linetypes` controls the values.

- ribbon_linetypes:

  Optional named linetype vector for `ribbon_linetype_by`.

- ribbon_direction:

  Character. How to visually distinguish same- vs reverse-orientation
  alignments: `"none"`, `"alpha"`, `"outline"` or `"linetype"`.

- ribbon_direction_colors:

  Named color vector with `same` and `reverse` entries, used when
  `ribbon_direction = "outline"`.

- ribbon_direction_linetypes:

  Named linetype vector with `same` and `reverse` entries, used when
  `ribbon_direction = "linetype"`.

- ribbon_direction_alpha:

  Named numeric vector with `same` and `reverse` entries, used when
  `ribbon_direction = "alpha"`.

- show_legend:

  Whether to show the legend, default TRUE

- legend_position:

  Position of this layer's legend (the Identity( colourbar): one of
  "left", "right", "top", "bottom" or "inside", default "left". Pass
  NULL to let the legend follow `theme(legend.position = ...)` together
  with the other legends.

- legend_key_width:

  Optional width of the Identity( Accepts a grid unit, e.g.
  `unit(1, "cm")`, or a number interpreted as centimetres. Default NULL
  uses the ggplot2 default width.

- legend_key_height:

  Optional height of the Identity( Accepts a grid unit, e.g.
  `unit(5, "cm")`, or a number interpreted as centimetres. Default NULL
  lets the vertical bar fill the available height (and uses a fixed
  height for horizontal legends).

- ...:

  Additional arguments passed to `geom_polygon()`

## Value

A list of ggplot2 layers

## Examples

``` r
library(ggchord)
data(seq_data_example)
data(ribbon_data_example)
p <- ggchord(seq_data_example, ribbon_data_example) +
  geom_seq() + geom_ribbon()
p
```
