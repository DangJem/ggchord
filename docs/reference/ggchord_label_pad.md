# Estimate the coordinate margin (in data units) needed so that the text labels rendered by the gene/sequence label layers stay inside the figure.

Kept for backwards compatibility. Plot limits are now computed
adaptively by
[`ggchord_adaptive_limits()`](https://dangjem.github.io/ggchord/reference/ggchord_adaptive_limits.md),
which fits the actual text boxes rather than adding a single
conservative margin on every side.

## Usage

``` r
ggchord_label_pad(layout)
```
