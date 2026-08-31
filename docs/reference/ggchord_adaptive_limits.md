# Compute coordinate limits that fit the rendered text boxes

Instead of adding one global text-width pad on every side, this helper
measures the actual gene, sequence and axis label boxes and expands only
the sides that need it. x and y are fitted independently:
\`coord_fixed()\` preserves equal physical units without requiring a
square data range. This lets wide or tall rendered content use the
available panel more efficiently.

## Usage

``` r
ggchord_adaptive_limits(layout)
```
