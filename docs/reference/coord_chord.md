# Chord diagram coordinate system

Controls global rotation, aspect ratio, coordinate limits, clipping and
the strategy used to fit chord geometry and labels.

## Usage

``` r
coord_chord(
  rotation = 45,
  ratio = 1,
  xlim = NULL,
  ylim = NULL,
  expand = FALSE,
  clip = "off",
  fit = c("labels", "geometry", "manual")
)
```

## Arguments

- rotation:

  Global clockwise layout rotation in degrees, default 45.

- ratio:

  Fixed y/x aspect ratio, default 1.

- xlim, ylim:

  Optional user limits. Explicit limits take priority over automatically
  fitted limits.

- expand:

  Logical. Expand coordinate limits, default FALSE. Automatic fitting
  already includes a small safety margin; set TRUE to request the
  additional ggplot2 coordinate expansion.

- clip:

  Whether drawing is clipped to the panel, default `"off"`.

- fit:

  Fitting strategy: `"labels"` includes measured label boxes,
  `"geometry"` fits geometric elements only, and `"manual"` requires
  explicit `xlim` and `ylim`.

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
