{{code:courseblockb}}

{{prepare:frontbackbone}}

## Rotate the plot {#orientation}

{{code:frontrotation}}

{{figure:frontrotation}}

`rotation` sets the counterclockwise angle and defaults to 0°. With `rotation = 90`, the backbones and labels turn 90° from the default plot. Compare the position of the gaps with the default plot. Rotation can free space on one side for layout while preserving coordinate order within each sequence.

## Increase the space between backbones {#parameters}

{{code:frontgap}}

{{figure:frontgap}}

`seq_gap` controls the gap proportion between adjacent backbones. The default is `0.03`. This example sets only the gap after MT108731.1 to `0.12`. Other gaps retain `0.03`. The last sequence’s gap closes the ring. Compare the space between adjacent arcs, then check that the longer sequence still has the longer arc.

## Radius: move one arc inward {#shape}

{{code:frontradius}}

{{figure:frontradius}}

The default radius is `1`. Three arcs keep it and `MT118296.1` uses `0.72`. Compare with the [default backbones](multi-sequence-start.html#draw). That arc sits closer to the centre and still ends at about 32,090 bp. Radius is a layout unit. Use the ticks to read sequence length.

## Curvature: draw a straight backbone {#curvature}

{{code:frontcurvature}}

{{figure:frontcurvature}}

`seq_curvature = 1` draws a standard arc, `0` draws a straight backbone, and negative values bend in the opposite direction. The other three backbones keep the standard arc. To combine radius and curvature, try each setting separately before using both.
