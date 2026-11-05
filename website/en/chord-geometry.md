The previous section adjusted overall rotation and backbone gaps. Use the same four sequences and change only `MT118296.1`. Each example starts from defaults so you can see the effect of one parameter.

## Radius: move one arc inward {#shape}

{{code:frontradius}}

{{figure:frontradius}}

The default radius is `1`. Three arcs keep it and `MT118296.1` uses `0.72`. Compare with the [default backbones](chord-start.html#draw). That arc sits closer to the centre and still ends at about 32,090 bp. Radius is a layout unit. Use the ticks to read sequence length.

## Curvature: draw a straight backbone {#curvature}

{{code:frontcurvature}}

{{figure:frontcurvature}}

`seq_curvature = 1` draws a standard arc, `0` draws a straight backbone, and negative values bend in the opposite direction. The other three backbones keep the standard arc. To combine radius and curvature, try each setting separately before using both.

## Orientation: start the coordinates at the other end {#direction}

{{code:frontorientation}}

{{figure:frontorientation}}

`seq_orientation = -1` reverses the sequence's display direction. Compared with the default, this backbone’s zero tick and arrow swap ends. Input coordinates and biological strands still come from the tables. The three vectors are named by sequence ID, so they still match if you reorder the sequences with `seq_order`.

## Order: arrange the backbones {#order}

{{code:frontorder}}

{{figure:frontorder}}

`seq_order` takes the sequence IDs in display order. Here `OQ646790.1` moves ahead of `MT118296.1`. Link rows still match the same IDs after reordering.

Next, [change the path and endpoint clearance of a ribbon](chord-routing.html).
