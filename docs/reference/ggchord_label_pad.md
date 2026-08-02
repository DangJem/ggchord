# Estimate the coordinate margin (in data units) needed so that the text labels rendered by the gene/sequence label layers stay inside the figure.

The measured text width (inches) is converted to data units with a
calibration factor calibrated for square figures (the standard ggchord
layout); the value is intentionally slightly conservative.

## Usage

``` r
ggchord_label_pad(layout)
```
