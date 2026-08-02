# De-overlap gene labels

Detects overlapping gene label boxes (estimated from the text size) and
pushes the labels apart until they no longer collide. Optionally hides
labels that still overlap more than \`max_overlaps\` other labels
(ggrepel-style decluttering).

## Usage

``` r
ggchord_label_deoverlap(
  gl,
  units_per_inch = 0.35,
  seed = 123,
  max_overlaps = Inf
)
```
