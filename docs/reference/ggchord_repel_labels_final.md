# Final deterministic de-overlap pass for repelled gene labels.

Runs after the horizontal/justification and arc-side adjustments so the
solver uses the exact rendered text boxes. It also treats sequence,
group and axis labels as hard rectangular obstacles.

## Usage

``` r
ggchord_repel_labels_final(
  gl,
  units_per_inch = 0.35,
  box_padding = 0.25,
  repel_boxes = NULL,
  max_iter = 500
)
```
