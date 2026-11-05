Use the same four sequences and change one parameter at a time. Compare the position of the backbones and the space between them with the default plot.

## Default plot {#prepare}

{{code:frontbackbone}}

{{figure:frontbackbone}}

## Rotate the plot {#orientation}

{{code:frontrotation}}

{{figure:frontrotation}}

`rotation` sets the counterclockwise angle and defaults to 0°. With `rotation = 90`, the backbones and labels turn 90° from the default plot. Compare the position of the gaps with the default plot. Rotation can free space on one side for layout while preserving coordinate order within each sequence.

## Increase the space between backbones {#parameters}

{{code:frontgap}}

{{figure:frontgap}}

`seq_gap` controls the gap proportion between adjacent backbones. The default is `0.03`. This example sets only the gap after MT108731.1 to `0.12`. Other gaps retain `0.03`. The last sequence’s gap closes the ring. Compare the space between adjacent arcs, then check that the longer sequence still has the longer arc.

## When do you need further adjustments? {#clearance}

If ribbons overlap arrows in a comparison, see [Ribbon clearance](chord-routing.html#clearance). To adjust the radius or bend of one backbone, see [Per-sequence parameters](parameters.html).

<span id="guides"></span>Next, [adjust radius, curvature and display direction per sequence](chord-geometry.html).
