# Ribbon fill scales

Ribbon fill scales

## Usage

``` r
scale_ribbon_fill_stepsn(
  ...,
  colours,
  values = NULL,
  colors,
  guide = ggplot2::guide_coloursteps(available_aes = "ribbon_fill")
)

scale_ribbon_fill_gradientn(
  ...,
  colours,
  values = NULL,
  colors,
  guide = ggplot2::guide_colourbar(available_aes = "ribbon_fill")
)

scale_ribbon_fill_manual(..., values)

scale_ribbon_fill_identity(...)
```

## Arguments

- ...:

  Arguments passed to the corresponding ggplot2 scale.

- colours, colors:

  Gradient colours.

- values:

  Manual values, or positions for gradient colours where supported by
  the corresponding ggplot2 scale.

- guide:

  A guide object; defaults accept the \`ribbon_fill\` aesthetic.

## Value

A ggplot2 scale.
