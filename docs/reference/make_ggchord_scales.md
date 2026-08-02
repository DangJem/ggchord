# Build the list of scales for a computed layout

Build the list of scales for a computed layout

## Usage

``` r
make_ggchord_scales(
  layout,
  has_seq = FALSE,
  has_gene = FALSE,
  legend_position = NULL,
  legend_box = NULL,
  positions = list(),
  legend_key_length = NULL
)
```

## Arguments

- legend_position:

  The plot theme's \`legend.position\` (character).

- legend_box:

  The plot theme's \`legend.box\` setting. When the legend is at the
  top/bottom or the legend box is laid out horizontally
  (\`"horizontal"\`), a \`unit(1, "null")\` colorbar key height
  collapses to zero height in ggplot2 (the Identity( used in those cases
  so the colorbar stays visible; otherwise the colorbar fills the
  available height.

- positions:

  Named list with per-legend position overrides (\`seq\`, \`ribbon\`,
  \`gene\`), each \`NULL\` or one of "left", "right", "top", "bottom",
  "inside". Overrides make that legend sit in its own legend box at the
  given position instead of following the theme's \`legend.position\`.
