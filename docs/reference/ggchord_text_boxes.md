# Estimate the axis-aligned text boxes for a set of text labels.

\`text_x\`/\`text_y\` are the points selected by \`hjust\`/\`vjust\`,
not the rendered text centre. This helper returns both the anchor
(\`x\`/\`y\`) and the centre (\`cx\`/\`cy\`) together with the full
axis-aligned width/height of the text box (\`bw\`/\`bh\`). The same
projection is used by the repulsion solver, the obstacle boxes and the
adaptive coordinate limits so all three agree.

## Usage

``` r
ggchord_text_boxes(
  df,
  x_col = "text_x",
  y_col = "text_y",
  text_col = "text",
  angle_col = "text_angle",
  size_col = "size",
  hjust_col = "hjust",
  vjust_col = "vjust",
  units_per_inch = 0.35,
  box_padding = 0
)
```
