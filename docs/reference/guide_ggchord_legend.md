# Guides for ggchord role aesthetics

Thin wrappers around \[ggplot2::guide_legend()\] and
\[ggplot2::guide_colourbar()\] with compact defaults suitable for chord
diagrams. All ordinary ggplot2 guide arguments remain available.

## Usage

``` r
guide_ggchord_legend(
  title = waiver(),
  theme = NULL,
  position = NULL,
  direction = NULL,
  override.aes = list(),
  nrow = NULL,
  ncol = NULL,
  reverse = FALSE,
  order = 0,
  ...
)

guide_ggchord_colourbar(
  title = waiver(),
  theme = NULL,
  nbin = NULL,
  display = "raster",
  alpha = NA,
  draw.ulim = TRUE,
  draw.llim = TRUE,
  angle = NULL,
  position = NULL,
  direction = NULL,
  reverse = FALSE,
  order = 0,
  available_aes = c("colour", "color", "fill", "ribbon_fill", "gene_fill",
    "feature_fill", "region_fill"),
  ...
)
```

## Arguments

- title:

  Guide title.

- theme:

  Optional guide-specific theme.

- position, direction:

  Guide position and direction.

- override.aes:

  A list of legend-key aesthetic overrides.

- nrow, ncol:

  Legend key layout.

- reverse:

  Whether to reverse the key order.

- order:

  Guide order.

- ...:

  Additional arguments passed to the corresponding ggplot2 guide.

- nbin:

  Number of colourbar bins.

- display:

  Colourbar display method.

- alpha:

  Colourbar alpha override.

- draw.ulim, draw.llim:

  Whether to draw upper/lower limit ticks.

- angle:

  Tick-label angle.

- available_aes:

  Aesthetics accepted by the guide.

## Value

A ggplot2 guide object.
