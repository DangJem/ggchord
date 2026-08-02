# Hide text labels that overlap the plot content or each other

Estimates each label box (using the measured text size) and sets `label`
to NA when the box overlaps the given content points or another label
box.

## Usage

``` r
ggchord_hide_text_overlaps(df, content_pts, units_per_inch = 0.35)
```
