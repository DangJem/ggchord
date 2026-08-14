# Compute coordinate limits that fit the rendered text boxes

The chord geometry is placed in a square, fixed-aspect panel. Instead of
adding one global text-width pad on every side, this helper measures the
actual gene/sequence/group/axis label boxes and expands only the sides
that need it. The result is a tighter plot that uses the available panel
area.

## Usage

``` r
ggchord_adaptive_limits(layout)
```
