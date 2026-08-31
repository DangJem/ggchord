# Feature shape scale

Maps feature categories to the four geometry types understood by
\[geom_feature()\]: \`"arrow"\`, \`"block"\`, \`"chevron"\`, and
\`"lollipop"\`. Shape values affect the actual feature geometry as well
as its legend key.

## Usage

``` r
scale_feature_shape_manual(
  ...,
  values,
  name = "Feature",
  limits = NULL,
  guide = "none"
)
```

## Arguments

- ...:

  Arguments passed to \[ggplot2::discrete_scale()\].

- values:

  A named or unnamed character vector of feature geometry names.

- name:

  Scale and guide title, default \`"Feature"\`.

- limits:

  Optional category order. Named \`values\` use their name order by
  default so shape and fill guides can merge cleanly.

- guide:

  Guide specification. The default is \`"none"\` because
  \[geom_feature()\] combines shape silhouettes into its fill legend
  when both aesthetics describe the same feature categories. Supply
  \[guide_ggchord_legend()\] for a separate shape guide.

## Value

A ggplot2 discrete scale for the \`feature_shape\` aesthetic.
