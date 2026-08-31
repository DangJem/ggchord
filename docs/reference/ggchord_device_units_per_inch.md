# Convert physical text dimensions to the current fixed-aspect plot scale.

Text is rendered in millimetres, whereas chord geometry is expressed in
data units. A fixed data-units-per-inch constant therefore cannot
describe the same label on both a small and a large output device. This
helper uses the current device dimensions and the undecorated chord
span; importantly, it does not feed already-expanded label limits back
into the estimate. The latter used to make leader-line clipping grow
with the labels themselves and produced conspicuously large,
output-size-dependent gaps.

## Usage

``` r
ggchord_device_units_per_inch(x, y, fallback_inches = 6, margin_inches = 1.25)
```
