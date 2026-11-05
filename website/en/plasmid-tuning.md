Compare the default circular map with a single backbone and a rotated map.

## Default double backbone {#baseline}

{{code:firstcircle}}

{{figure:firstcircle}}

## Switch to a single backbone {#style}

{{code:frontcirclestyle}}

{{figure:frontcirclestyle}}

`seq_style = "single"` replaces the double backbone with one line. The centre label and bp ticks remain visible.

## Rotate the map {#rotation}

{{code:firstfeatures,frontcirclerotation}}

{{figure:frontcirclerotation}}

The default `rotation = 0` places the sequence origin at 12 o’clock, with coordinates increasing clockwise. Here `rotation = 90` moves it to 9 o’clock, rotating the ticks and replication origin together. The annotation still covers the same source bp interval. Next, [prepare your own sequence and annotations](plasmid-data.html).
